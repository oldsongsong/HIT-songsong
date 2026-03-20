function indicator = pb_refinement_indicator(r, dr, gap, par)
% PB_REFINEMENT_INDICATOR Combined geometric/physical refinement metric.
[rhoP, rhoB] = pb_marker_density(r, gap, par);
combined = rhoP .* rhoB;
combined = combined / max(max(combined), eps);

hLocal = pb_film_thickness(r, gap, par);
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
