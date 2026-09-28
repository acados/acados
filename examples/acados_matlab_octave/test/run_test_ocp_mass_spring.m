%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.

%

addpath(pwd)

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
    test_ocp_linear_mass_spring;
    test_ocp_linear_mass_spring_new;
catch exception
    exit_with_error(exception);
end

fprintf('\nrun_tests_ocp: success!\n\n');
