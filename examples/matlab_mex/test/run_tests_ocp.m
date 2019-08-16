%% check that environment variables are provided
require_env_variable('LD_LIBRARY_PATH');
require_env_variable('ACADOS_INSTALL_DIR');

if is_octave()
    require_env_variable('OCTAVE_PATH');
else
    require_env_variable('MATLABPATH');
end

%% ocp tests
try
    test_ocp;
catch error
    exit_with_error(error);
end

fprintf('\nrun_tests_ocp: success!\n\n');
