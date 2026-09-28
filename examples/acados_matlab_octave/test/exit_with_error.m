%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.

%

function exit_with_error(error)
    disp('in exit_with_error')
    disp(error)
    disp(error.message)
    fprintf(error.message);
    fprintf('\nRUN_TEST: FAIL -> exit\n');
    error('')
end