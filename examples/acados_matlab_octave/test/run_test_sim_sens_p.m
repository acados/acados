%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.

%

%% check that environment variables are provided

addpath(pwd)

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

%% sim tests
try
    test_sens_forw_p;
catch exception
    exit_with_error(exception);
end
