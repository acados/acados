%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.

%


example_dir = fileparts(which('acados_env_variables_windows'));

acados_dir = fullfile(example_dir, '..', '..');
casadi_dir = fullfile(acados_dir, 'external', 'casadi-matlab');
matlab_interface_dir = fullfile(acados_dir, 'interfaces', 'acados_matlab_octave');

addpath(matlab_interface_dir);
addpath(casadi_dir);

setenv('ACADOS_INSTALL_DIR', acados_dir);
setenv('ENV_RUN', 'true');