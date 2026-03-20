function Frep = pb_disjoining_force(gap, par)
% PB_DISJOINING_FORCE Weak repulsive regularization for near-contact stability.
heff = max(gap, par.hFloor);
Frep = par.hamaker * par.Reff / (6.0 * heff^2);
end
