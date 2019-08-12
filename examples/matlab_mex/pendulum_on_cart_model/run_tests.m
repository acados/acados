try
    % check that env.sh has been run
    env_run = getenv('ENV_RUN');
    if (~strcmp(env_run, 'true'))
        error('ERROR: env.sh has not been sourced! Before executing this example.');
    end
catch error
    exit_with_error(error);
end



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

% exit(0);
