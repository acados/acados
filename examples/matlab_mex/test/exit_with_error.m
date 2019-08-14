function exit_with_error(error)
    fprintf(error.message);
    fprintf('\nRUN_TEST: FAIL -> exit\n');
    quit(1);
end