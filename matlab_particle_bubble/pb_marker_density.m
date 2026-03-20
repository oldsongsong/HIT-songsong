function [rhoP, rhoB] = pb_marker_density(r, gap, par)
% PB_MARKER_DENSITY Project particle/bubble Lagrangian markers onto radial grid.

theta = linspace(0.0, par.thetaMax, par.nMarkers);
rp = par.Rp * sin(theta);
rb = par.Rb * sin(theta);

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

proximityBoost = exp(-gap / par.hRefine);
rhoP = rhoP * (1.0 + 0.5 * proximityBoost);
rhoB = rhoB * (1.0 + 0.5 * proximityBoost);
end
