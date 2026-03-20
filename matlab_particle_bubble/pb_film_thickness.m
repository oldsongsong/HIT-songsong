function h = pb_film_thickness(r, gap, par)
% PB_FILM_THICKNESS Axisymmetric paraboloid thin-film thickness model.
h = max(gap, par.hFloor) + 0.5 * (r.^2) / par.Reff;
end
