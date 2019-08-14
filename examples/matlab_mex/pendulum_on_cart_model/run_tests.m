%% check that environment variables are provided
require_env_variable('LD_LIBRARY_PATH');
require_env_variable('ACADOS_INSTALL_DIR');

if is_octave()
    require_env_variable('OCTAVE_PATH');
else
    require_env_variable('MATLABPATH');
end


%% sim tests
try
    test_sens_forw_ci;
%     test_sens_adj;
%     test_sens_hess;

%     c = 1;
%     b = [c zeros(2)];
catch error
    exit_with_error(error);
end

fprintf('\nRUN_TESTS: success!\n\n');
