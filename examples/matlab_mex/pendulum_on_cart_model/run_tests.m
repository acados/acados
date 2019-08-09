% 1/inf

try
    test_sens_forw;
%     test_sens_adj;
%     test_sens_hess;

%     c = 1;
%     b = [c zeros(2)];
catch
    fprintf('\nRUN_TEST: FAIL -> exit\n');
    quit(1);
end

fprintf('\nRUN_TESTS: success!\n\n');

exit()
