%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.

%

function check_casadi_availibility()
    try
        import casadi.*
        disp('casadi import success');
    catch error
        exit_with_error(error);
    end
end