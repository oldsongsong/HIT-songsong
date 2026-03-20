function sol = pb_solve_reynolds_pressure(grid, gap, urel, par)
% PB_SOLVE_REYNOLDS_PRESSURE Solve axisymmetric Reynolds lubrication equation.
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
        hw  = pb_film_thickness(rw, gap, par);
        aw  = rw * hw^3 / dxw;
        A(i,i-1) =  aw / vol;
    end

    if i < N
        re = rf(i+1);
        dxe = r(i+1) - r(i);
        he  = pb_film_thickness(re, gap, par);
        ae  = re * he^3 / dxe;
        A(i,i+1) =  ae / vol;
    else
        re = rf(end);
        dxe = rf(end) - r(end);
        he  = pb_film_thickness(re, gap, par);
        ae  = re * he^3 / dxe;
    end

    A(i,i) = -(aw + ae) / vol;
end

p = A \ b;
p = max(p, 0.0);
force = 2.0 * pi * sum(p .* r .* dr);

sol.pressure = p;
sol.force    = force;
end
