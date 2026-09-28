%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.

%

function print_casadi_expression(f)
    for ii = 1:length(f)
        disp(f(ii,:));
    end
    disp(' ');
end