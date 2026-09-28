%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.

%

addpath(pwd)
addpath(fullfile(pwd, '..', 'test'))

%% check that environment variables are provided
try
    check_casadi_availibility();
    require_env_variable('LD_LIBRARY_PATH');
    require_env_variable('ACADOS_INSTALL_DIR');
    if is_octave()
        require_env_variable('OCTAVE_PATH');
    else
        require_env_variable('MATLABPATH');
    end
catch exception
    exit_with_error(exception);
end



%% ocp tests
try
    disc_dyn_example_ocp;
catch exception
    exit_with_error(exception);
end