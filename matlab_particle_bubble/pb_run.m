function results = pb_run(mode)
% PB_RUN Main entry for the modular particle-bubble drainage framework.
%
% Usage:
%   results = pb_run();
%   results = pb_run('validation');
%   results = pb_run('dynamic');

if nargin < 1
    mode = 'all';
end

par = pb_default_parameters();

switch lower(mode)
    case 'validation'
        results = pb_run_validation_case(par, true);
    case 'dynamic'
        results = pb_run_dynamic_case(par, true);
    otherwise
        results.validation = pb_run_validation_case(par, true);
        results.dynamic    = pb_run_dynamic_case(par, true);
end
end
