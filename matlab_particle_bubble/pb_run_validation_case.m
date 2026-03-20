function validation = pb_run_validation_case(par, doPlot)
% PB_RUN_VALIDATION_CASE Compare numerical drainage force with Taylor theory.
gapList = par.validationGapList(:);
nCases  = numel(gapList);

Fnum       = zeros(nCases,1);
Ftheory    = zeros(nCases,1);
relError   = zeros(nCases,1);
nCells     = zeros(nCases,1);
minDr      = zeros(nCases,1);
sampleGrid = struct();

for k = 1:nCases
    gap  = gapList(k);
    grid = pb_build_adaptive_grid(gap, par);
    sol  = pb_solve_reynolds_pressure(grid, gap, par.validationUrel, par);

    Fnum(k)      = sol.force;
    Ftheory(k)   = pb_taylor_lubrication_force(gap, par.validationUrel, par);
    relError(k)  = abs(Fnum(k) - Ftheory(k)) / max(abs(Ftheory(k)), eps);
    nCells(k)    = numel(grid.centers);
    minDr(k)     = min(grid.widths);

    if isempty(fieldnames(sampleGrid)) && gap <= par.showGridSampleAtGap
        sampleGrid = grid;
        sampleGrid.gap = gap;
        sampleGrid.pressure = sol.pressure;
    end
end

validation.gap = gapList;
validation.forceNumerical   = Fnum;
validation.forceTheory      = Ftheory;
validation.relativeError    = relError;
validation.nCells           = nCells;
validation.minDr            = minDr;
validation.sampleGrid       = sampleGrid;
validation.meanRelativeError = mean(relError);
validation.maxRelativeError  = max(relError);

if doPlot
    pb_plot_validation_results(validation);
end

fprintf('Validation mean relative error  = %.3f %%\n', 100.0*validation.meanRelativeError);
fprintf('Validation max relative error   = %.3f %%\n', 100.0*validation.maxRelativeError);
end
