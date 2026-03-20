function par = pb_default_parameters()
% PB_DEFAULT_PARAMETERS Default physical, numerical and adaptive-grid parameters.

% Physical properties.
par.mu   = 1.0e-3;
par.rhoL = 1000.0;
par.rhoP = 2500.0;
par.rhoB = 1.2;
par.g    = 9.81;

% Geometry.
par.Rp   = 0.50e-3;
par.Rb   = 0.80e-3;
par.Reff = (par.Rp * par.Rb) / (par.Rp + par.Rb);
par.Rmax = 0.45e-3;

% Marker-driven refinement parameters.
par.nMarkers      = 240;
par.thetaMax      = 22.0 * pi / 180.0;
par.markerSigma   = 0.80;
par.lagWeightMix  = 0.60;
par.filmWeightMix = 0.30;
par.bandWeightMix = 0.10;
par.hRefine       = 3.0e-6;
par.bandFactor    = 3.0;

% Adaptive-grid controls.
par.nBase            = 20;
par.maxLevel         = 5;
par.refineThresholds = [0.48, 0.42, 0.36, 0.31, 0.27];
par.targetCellsOnGap = 8.0;
par.hFloor           = 0.05e-6;

% Validation settings.
par.validationGapList = logspace(log10(0.25e-6), log10(5.0e-6), 9);
par.validationUrel    = 1.0e-3;

% Quasi-static approach settings.
Vp = 4.0/3.0*pi*par.Rp^3;
Vb = 4.0/3.0*pi*par.Rb^3;
par.mp = par.rhoP * Vp;
par.mb = par.rhoB * Vb;
par.Fdrive            = (par.rhoP - par.rhoL)*Vp*par.g + (par.rhoL - par.rhoB)*Vb*par.g;
par.dynamicDt         = 2.0e-6;
par.dynamicTend       = 6.0e-3;
par.initialGap        = 20.0e-6;
par.initialUrel       = 0.0;
par.hamaker           = 2.0e-20;
par.stopGap           = 0.10e-6;
par.showGridSampleAtGap = 1.0e-6;
end
