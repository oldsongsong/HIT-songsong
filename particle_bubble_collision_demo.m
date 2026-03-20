function results = particle_bubble_collision_demo(mode)
% PARTICLE_BUBBLE_COLLISION_DEMO
% -------------------------------------------------------------------------
% Research-oriented MATLAB prototype for particle-bubble thin-film drainage
% with a dynamic geometry-adaptive Cartesian grid.
%
% This code is intentionally focused on the lubrication-drainage stage that
% dominates particle-bubble 'collision' at very small gaps. It borrows the
% adaptive Cartesian-grid spirit of geometry-adaptive IB-LBM methods, but it
% solves the axisymmetric Reynolds lubrication equation in the thin film to
% keep the prototype compact and editable in MATLAB.
%
% Key ingredients
%   1) Axisymmetric thin-film drainage between a rigid particle and a bubble.
%   2) Dynamic non-uniform Cartesian grid in the radial direction.
%   3) Refinement indicator driven primarily by projected Lagrangian-marker
%      density, augmented by film-thickness information.
%   4) Prescribed-speed validation against the classical Taylor lubrication
%      force.
%   5) Coupled drainage dynamics with gravity/buoyancy driving force.
%
% Usage
%   results = particle_bubble_collision_demo();
%   results = particle_bubble_collision_demo('validation');
%   results = particle_bubble_collision_demo('dynamic');
%
% Output
%   A struct containing validation and/or dynamic drainage results.
% -------------------------------------------------------------------------

if nargin < 1
    mode = 'all';
end

par = defaultParameters();

switch lower(mode)
    case 'validation'
        results = runValidationCase(par, true);
    case 'dynamic'
        results = runDynamicCase(par, true);
    otherwise
        results.validation = runValidationCase(par, true);
        results.dynamic    = runDynamicCase(par, true);
end

end

function par = defaultParameters()
% Physical properties.
par.mu   = 1.0e-3;        % liquid viscosity [Pa*s]
par.rhoL = 1000.0;        % liquid density [kg/m^3]
par.rhoP = 2500.0;        % particle density [kg/m^3]
par.rhoB = 1.2;           % bubble gas density [kg/m^3]
par.g    = 9.81;          % gravity [m/s^2]

% Geometry.
par.Rp   = 0.50e-3;       % particle radius [m]
par.Rb   = 0.80e-3;       % bubble radius [m]
par.Reff = (par.Rp * par.Rb) / (par.Rp + par.Rb);
par.Rmax = 0.45e-3;       % radial truncation of film domain [m]

% Lagrangian-marker settings used by the adaptive-grid indicator.
par.nMarkers      = 240;
par.thetaMax      = 22.0 * pi / 180.0;
par.markerSigma   = 0.80;
par.lagWeightMix  = 0.60;
par.filmWeightMix = 0.30;
par.bandWeightMix = 0.10;
par.hRefine       = 3.0e-6;
par.bandFactor    = 3.0;

% Adaptive Cartesian grid parameters.
par.nBase            = 20;
par.maxLevel         = 5;
par.refineThresholds = [0.48, 0.42, 0.36, 0.31, 0.27];
par.targetCellsOnGap = 8.0;
par.hFloor           = 0.05e-6;

% Validation settings.
par.validationGapList = logspace(log10(0.25e-6), log10(5.0e-6), 9);
par.validationUrel    = 1.0e-3;     % prescribed approach speed [m/s]

% Dynamic approach/collision settings.
Vp = 4.0/3.0*pi*par.Rp^3;
Vb = 4.0/3.0*pi*par.Rb^3;
par.mp = par.rhoP * Vp;
par.mb = par.rhoB * Vb;
par.Fdrive = (par.rhoP - par.rhoL)*Vp*par.g + (par.rhoL - par.rhoB)*Vb*par.g;
par.dynamicDt         = 2.0e-6;
par.dynamicTend       = 6.0e-3;
par.initialGap        = 20.0e-6;
par.initialUrel       = 0.0;
par.hamaker           = 2.0e-20;    % weak repulsion regularization [J]
par.stopGap           = 0.10e-6;

% Plotting.
par.showGridSampleAtGap = 1.0e-6;
end

function validation = runValidationCase(par, doPlot)
gapList = par.validationGapList(:);
nCases  = numel(gapList);

Fnum      = zeros(nCases,1);
Ftheory   = zeros(nCases,1);
relError  = zeros(nCases,1);
nCells    = zeros(nCases,1);
minDr     = zeros(nCases,1);
sampleGrid = struct();

for k = 1:nCases
    gap  = gapList(k);
    grid = buildAdaptiveGrid(gap, par);
    sol  = solveLubricationPressure(grid, gap, par.validationUrel, par);

    Fnum(k)    = sol.force;
    Ftheory(k) = taylorLubricationForce(gap, par.validationUrel, par);
    relError(k) = abs(Fnum(k) - Ftheory(k)) / max(abs(Ftheory(k)), eps);
    nCells(k) = numel(grid.centers);
    minDr(k)  = min(grid.widths);

    if isempty(fieldnames(sampleGrid)) && gap <= par.showGridSampleAtGap
        sampleGrid = grid;
        sampleGrid.gap = gap;
        sampleGrid.pressure = sol.pressure;
    end
end

validation.gap = gapList;
validation.forceNumerical = Fnum;
validation.forceTheory    = Ftheory;
validation.relativeError  = relError;
validation.nCells         = nCells;
validation.minDr          = minDr;
validation.sampleGrid     = sampleGrid;
validation.meanRelativeError = mean(relError);
validation.maxRelativeError  = max(relError);

if doPlot
    plotValidationResults(validation);
end

fprintf('Validation mean relative error  = %.3f %%\n', 100.0*validation.meanRelativeError);
fprintf('Validation max relative error   = %.3f %%\n', 100.0*validation.maxRelativeError);

end

function dynamic = runDynamicCase(par, doPlot)
Nt = ceil(par.dynamicTend / par.dynamicDt);

time    = zeros(Nt,1);
gap     = zeros(Nt,1);
urel    = zeros(Nt,1);
Fhydro  = zeros(Nt,1);
Ftheory = zeros(Nt,1);
Fdlvo   = zeros(Nt,1);
nCells  = zeros(Nt,1);
minDr   = zeros(Nt,1);

currentGap  = par.initialGap;
currentUrel = par.initialUrel;

sampleGrid = struct();
last = Nt;

for n = 1:Nt
    grid = buildAdaptiveGrid(currentGap, par);
    solUnit = solveLubricationPressure(grid, currentGap, 1.0, par);
    Frep = disjoiningForce(currentGap, par);

    resistance = max(solUnit.force, eps);
    currentUrel = max((par.Fdrive - Frep) / resistance, 0.0);
    currentFhydro = resistance * currentUrel;
    nextGap  = max(par.hFloor, currentGap - par.dynamicDt * currentUrel);

    time(n)    = (n-1) * par.dynamicDt;
    gap(n)     = currentGap;
    urel(n)    = currentUrel;
    Fhydro(n)  = currentFhydro;
    Ftheory(n) = taylorLubricationForce(currentGap, currentUrel, par);
    Fdlvo(n)   = Frep;
    nCells(n)  = numel(grid.centers);
    minDr(n)   = min(grid.widths);

    if isempty(fieldnames(sampleGrid)) && currentGap <= par.showGridSampleAtGap
        sampleGrid = grid;
        sampleGrid.gap = currentGap;
        sampleGrid.pressure = solveLubricationPressure(grid, currentGap, currentUrel, par).pressure;
    end

    currentGap  = nextGap;

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
    plotDynamicResults(dynamic, par);
end

fprintf('Dynamic run final gap          = %.3e m\n', dynamic.terminalGap);
fprintf('Dynamic run final approach vel = %.3e m/s\n', dynamic.finalApproachSpeed);

end

function grid = buildAdaptiveGrid(gap, par)
faces = linspace(0.0, par.Rmax, par.nBase + 1);

for level = 1:par.maxLevel
    centers = 0.5 * (faces(1:end-1) + faces(2:end));
    widths  = diff(faces);
    indicator = refinementIndicator(centers, widths, gap, par);

    splitMask = indicator > par.refineThresholds(level);
    if ~any(splitMask)
        break;
    end

    newFaces = faces(1);
    for i = 1:numel(widths)
        if splitMask(i)
            mid = 0.5 * (faces(i) + faces(i+1));
            newFaces = [newFaces, mid, faces(i+1)]; %#ok<AGROW>
        else
            newFaces = [newFaces, faces(i+1)]; %#ok<AGROW>
        end
    end
    faces = unique(newFaces, 'stable');
end

centers = 0.5 * (faces(1:end-1) + faces(2:end));
widths  = diff(faces);
[rhoP, rhoB] = markerDensity(centers, gap, par);

grid.faces   = faces(:);
grid.centers = centers(:);
grid.widths  = widths(:);
grid.rhoP    = rhoP(:);
grid.rhoB    = rhoB(:);
grid.indicator = refinementIndicator(centers, widths, gap, par);
end

function indicator = refinementIndicator(r, dr, gap, par)
[rhoP, rhoB] = markerDensity(r, gap, par);
combined = rhoP .* rhoB;
combined = combined / max(max(combined), eps);

hLocal = filmThickness(r, gap, par);
filmMetric = exp(-hLocal / par.hRefine);
filmMetric = filmMetric / max(max(filmMetric), eps);

bandRadius = par.bandFactor * sqrt(2.0 * par.Reff * max(gap, par.hFloor));
bandMetric = double(r <= bandRadius);

resolutionMetric = min(1.0, dr ./ (par.targetCellsOnGap * max(hLocal, par.hFloor)));

indicator = par.lagWeightMix  * combined + ...
            par.filmWeightMix * filmMetric + ...
            par.bandWeightMix * bandMetric + ...
            0.25 * resolutionMetric;
end

function [rhoP, rhoB] = markerDensity(r, gap, par)
theta = linspace(0.0, par.thetaMax, par.nMarkers);

rp = par.Rp * sin(theta);
rb = par.Rb * sin(theta);

% Equal arclength weights on each surface meridian.
dsp = [diff(theta), theta(end)-theta(end-1)] * par.Rp;
dsb = [diff(theta), theta(end)-theta(end-1)] * par.Rb;

baseDx = par.Rmax / par.nBase;
sigma  = par.markerSigma * baseDx;

rhoP = zeros(size(r));
rhoB = zeros(size(r));
for k = 1:numel(theta)
    rhoP = rhoP + dsp(k) * exp(-0.5*((r-rp(k))/sigma).^2) / (sqrt(2.0*pi) * sigma);
    rhoB = rhoB + dsb(k) * exp(-0.5*((r-rb(k))/sigma).^2) / (sqrt(2.0*pi) * sigma);
end

% When the film gets thinner, overlapping projected densities should be more
% important in the central drainage region.
proximityBoost = exp(-gap / par.hRefine);
rhoP = rhoP * (1.0 + 0.5 * proximityBoost);
rhoB = rhoB * (1.0 + 0.5 * proximityBoost);
end

function h = filmThickness(r, gap, par)
h = max(gap, par.hFloor) + 0.5 * (r.^2) / par.Reff;
end

function sol = solveLubricationPressure(grid, gap, urel, par)
N  = numel(grid.centers);
r  = grid.centers;
dr = grid.widths;
rf = grid.faces;

A = spalloc(N, N, 3*N);
b = -12.0 * par.mu * urel * ones(N,1);

for i = 1:N
    rc = max(r(i), 0.5*dr(i));
    vol = rc * dr(i);

    aw = 0.0;
    ae = 0.0;

    if i > 1
        rw = rf(i);
        dxw = r(i) - r(i-1);
        hw  = filmThickness(rw, gap, par);
        aw  = rw * hw^3 / dxw;
        A(i,i-1) =  aw / vol;
    end

    if i < N
        re = rf(i+1);
        dxe = r(i+1) - r(i);
        he  = filmThickness(re, gap, par);
        ae  = re * he^3 / dxe;
        A(i,i+1) =  ae / vol;
    else
        re = rf(end);
        dxe = rf(end) - r(end);
        he  = filmThickness(re, gap, par);
        ae  = re * he^3 / dxe;
        % Dirichlet p(Rmax)=0, so boundary contribution only modifies diag.
    end

    A(i,i) = -(aw + ae) / vol;
    b(i) = b(i);
end

p = A \ b;
p = max(p, 0.0);
force = 2.0 * pi * sum(p .* r .* dr);

sol.pressure = p;
sol.force    = force;
end

function Ft = taylorLubricationForce(gap, urel, par)
Ft = 6.0 * pi * par.mu * par.Reff^2 * urel / max(gap, par.hFloor);
end

function Frep = disjoiningForce(gap, par)
heff = max(gap, par.hFloor);
Frep = par.hamaker * par.Reff / (6.0 * heff^2);
end

function plotValidationResults(validation)
figure('Color','w','Position',[100 100 1100 800]);

subplot(2,2,1);
loglog(validation.gap*1e6, validation.forceNumerical, 'o-', 'LineWidth', 1.5, 'MarkerSize', 6); hold on;
loglog(validation.gap*1e6, validation.forceTheory, 's--', 'LineWidth', 1.5, 'MarkerSize', 6);
xlabel('gap [\mum]'); ylabel('drainage force [N]');
legend('numerical','Taylor theory','Location','southwest');
title('Hydrodynamic drainage force'); grid on;

subplot(2,2,2);
semilogx(validation.gap*1e6, 100.0*validation.relativeError, 'o-', 'LineWidth', 1.5, 'MarkerSize', 6);
xlabel('gap [\mum]'); ylabel('relative error [%]');
title('Validation error'); grid on;

subplot(2,2,3);
semilogx(validation.gap*1e6, validation.nCells, 'o-', 'LineWidth', 1.5, 'MarkerSize', 6); hold on;
semilogx(validation.gap*1e6, validation.minDr*1e6, 's--', 'LineWidth', 1.5, 'MarkerSize', 6);
xlabel('gap [\mum]'); ylabel('cell count / min \Deltar [\mum]');
legend('cell count','min \Deltar','Location','best');
title('Adaptive-grid response'); grid on;

subplot(2,2,4);
if isfield(validation.sampleGrid, 'centers')
    stairs(validation.sampleGrid.faces*1e6, [validation.sampleGrid.pressure(1); validation.sampleGrid.pressure; 0], 'LineWidth', 1.2); hold on;
    yyaxis right;
    plot(validation.sampleGrid.centers*1e6, validation.sampleGrid.indicator, 'k-', 'LineWidth', 1.5);
    ylabel('refinement indicator [-]');
    yyaxis left;
    ylabel('pressure [Pa]');
    xlabel('r [\mum]');
    title(sprintf('Sample adaptive grid at gap = %.2f \\mum', validation.sampleGrid.gap*1e6));
    grid on;
else
    text(0.1, 0.5, 'No sample grid stored.', 'FontSize', 12);
    axis off;
end

sgtitle('Particle-bubble lubrication validation with adaptive Cartesian grid');
end

function plotDynamicResults(dynamic, par)
figure('Color','w','Position',[120 120 1100 800]);

subplot(2,2,1);
plot(dynamic.time*1e3, dynamic.gap*1e6, 'b-', 'LineWidth', 1.5);
xlabel('time [ms]'); ylabel('minimum gap [\mum]');
title('Gap evolution during drainage'); grid on;

subplot(2,2,2);
plot(dynamic.time*1e3, dynamic.urel*1e3, 'r-', 'LineWidth', 1.5);
xlabel('time [ms]'); ylabel('approach speed [mm/s]');
title('Relative approach speed'); grid on;

subplot(2,2,3);
plot(dynamic.time*1e3, dynamic.Fhydro*1e6, 'k-', 'LineWidth', 1.5); hold on;
plot(dynamic.time*1e3, dynamic.Ftheory*1e6, 'k--', 'LineWidth', 1.2);
plot(dynamic.time*1e3, dynamic.Fdlvo*1e6, 'm-.', 'LineWidth', 1.2);
yline(par.Fdrive*1e6, 'b:', 'LineWidth', 1.5);
xlabel('time [ms]'); ylabel('force [\muN]');
legend('numerical hydro','Taylor hydro','disjoining','driving force','Location','best');
title('Force balance'); grid on;

subplot(2,2,4);
plot(dynamic.time*1e3, dynamic.nCells, 'g-', 'LineWidth', 1.5); hold on;
plot(dynamic.time*1e3, dynamic.minDr*1e6, 'c--', 'LineWidth', 1.5);
xlabel('time [ms]'); ylabel('cell count / min \Deltar [\mum]');
legend('cell count','min \Deltar','Location','best');
title('Dynamic adaptive-grid response'); grid on;

sgtitle('Adaptive particle-bubble drainage simulation');
end
