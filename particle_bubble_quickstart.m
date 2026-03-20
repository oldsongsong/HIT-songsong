% PARTICLE_BUBBLE_QUICKSTART
% -------------------------------------------------------------------------
% Minimal calling examples for the modular particle-bubble MATLAB framework.
%
% Recommended usage:
%   1) Run this file section-by-section in MATLAB.
%   2) Start from the wrapper call examples.
%   3) Only after that, modify parameters through the direct-module examples.
% -------------------------------------------------------------------------

%% Example 1: simplest call (run everything)
results_all = particle_bubble_collision_demo();

%% Example 2: run validation only
results_validation = particle_bubble_collision_demo('validation');

%% Example 3: run dynamic drainage only
results_dynamic = particle_bubble_collision_demo('dynamic');

%% Example 4: direct modular usage with custom parameters
rootDir = fileparts(mfilename('fullpath'));
frameworkDir = fullfile(rootDir, 'matlab_particle_bubble');
addpath(frameworkDir);

par = pb_default_parameters();
par.validationUrel = 2.0e-3;          % user-defined approach speed [m/s]
par.initialGap     = 15.0e-6;         % user-defined initial gap [m]
par.dynamicTend    = 4.0e-3;          % user-defined end time [s]

validation_custom = pb_run_validation_case(par, true);
dynamic_custom    = pb_run_dynamic_case(par, true);

%% Example 5: lowest-level call (grid + Reynolds solver only)
gap  = 1.0e-6;
urel = 1.0e-3;

grid = pb_build_adaptive_grid(gap, par);
sol  = pb_solve_reynolds_pressure(grid, gap, urel, par);

fprintf('Manual call: gap = %.3e m, force = %.3e N\n', gap, sol.force);
