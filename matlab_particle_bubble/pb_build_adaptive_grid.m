function grid = pb_build_adaptive_grid(gap, par)
% PB_BUILD_ADAPTIVE_GRID Binary refinement of the radial Cartesian grid.
faces = linspace(0.0, par.Rmax, par.nBase + 1);

for level = 1:par.maxLevel
    centers = 0.5 * (faces(1:end-1) + faces(2:end));
    widths  = diff(faces);
    indicator = pb_refinement_indicator(centers, widths, gap, par);

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
[rhoP, rhoB] = pb_marker_density(centers, gap, par);

grid.faces     = faces(:);
grid.centers   = centers(:);
grid.widths    = widths(:);
grid.rhoP      = rhoP(:);
grid.rhoB      = rhoB(:);
grid.indicator = pb_refinement_indicator(centers, widths, gap, par);
end
