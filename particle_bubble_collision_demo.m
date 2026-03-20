function results = particle_bubble_collision_demo(mode)
% PARTICLE_BUBBLE_COLLISION_DEMO
% -------------------------------------------------------------------------
% Compatibility wrapper for the modular MATLAB particle-bubble drainage
% framework stored in ./matlab_particle_bubble.
%
% This wrapper keeps the old entry name while redirecting execution to the
% reorganized, easier-to-extend research framework.
% -------------------------------------------------------------------------

if nargin < 1
    mode = 'all';
end

thisFile = mfilename('fullpath');
thisDir  = fileparts(thisFile);
frameworkDir = fullfile(thisDir, 'matlab_particle_bubble');
if exist(frameworkDir, 'dir') == 7
    addpath(frameworkDir);
else
    error('Framework directory not found: %s', frameworkDir);
end

results = pb_run(mode);
end
