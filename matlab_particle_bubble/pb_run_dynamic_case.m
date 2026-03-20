function dynamic = pb_run_dynamic_case(par, doPlot)
% PB_RUN_DYNAMIC_CASE Quasi-static gap evolution under driving and repulsion.
Nt = ceil(par.dynamicTend / par.dynamicDt);

time    = zeros(Nt,1);
gap     = zeros(Nt,1);
urel    = zeros(Nt,1);
Fhydro  = zeros(Nt,1);
Ftheory = zeros(Nt,1);
Fdlvo   = zeros(Nt,1);
nCells  = zeros(Nt,1);
minDr   = zeros(Nt,1);

currentGap = par.initialGap;
sampleGrid = struct();
last = Nt;

for n = 1:Nt
    grid = pb_build_adaptive_grid(currentGap, par);
    solUnit = pb_solve_reynolds_pressure(grid, currentGap, 1.0, par);
    Frep = pb_disjoining_force(currentGap, par);

    resistance = max(solUnit.force, eps);
    currentUrel = max((par.Fdrive - Frep) / resistance, 0.0);
    currentFhydro = resistance * currentUrel;
    nextGap = max(par.hFloor, currentGap - par.dynamicDt * currentUrel);

    time(n)    = (n-1) * par.dynamicDt;
    gap(n)     = currentGap;
    urel(n)    = currentUrel;
    Fhydro(n)  = currentFhydro;
    Ftheory(n) = pb_taylor_lubrication_force(currentGap, currentUrel, par);
    Fdlvo(n)   = Frep;
    nCells(n)  = numel(grid.centers);
    minDr(n)   = min(grid.widths);

    if isempty(fieldnames(sampleGrid)) && currentGap <= par.showGridSampleAtGap
        sampleGrid = grid;
        sampleGrid.gap = currentGap;
        sampleGrid.pressure = pb_solve_reynolds_pressure(grid, currentGap, currentUrel, par).pressure;
    end

    currentGap = nextGap;
    if currentGap <= par.stopGap
        last = n;
        break;
    end
end

dynamic.time    = time(1:last);
dynamic.gap     = gap(1:last);
dynamic.urel    = urel(1:last);
dynamic.Fhydro  = Fhydro(1:last);
dynamic.Ftheory = Ftheory(1:last);
dynamic.Fdlvo   = Fdlvo(1:last);
dynamic.nCells  = nCells(1:last);
dynamic.minDr   = minDr(1:last);
dynamic.sampleGrid = sampleGrid;
dynamic.terminalGap = dynamic.gap(end);
dynamic.finalApproachSpeed = dynamic.urel(end);

if doPlot
    pb_plot_dynamic_results(dynamic, par);
end

fprintf('Dynamic run final gap          = %.3e m\n', dynamic.terminalGap);
fprintf('Dynamic run final approach vel = %.3e m/s\n', dynamic.finalApproachSpeed);
end
