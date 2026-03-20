function Ft = pb_taylor_lubrication_force(gap, urel, par)
% PB_TAYLOR_LUBRICATION_FORCE Small-gap Taylor lubrication reference force.
Ft = 6.0 * pi * par.mu * par.Reff^2 * urel / max(gap, par.hFloor);
end
