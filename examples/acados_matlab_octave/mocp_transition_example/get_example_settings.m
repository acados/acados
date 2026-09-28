%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.

function settings = get_example_settings()
    settings.X0 = [2.0, 0.0];
    settings.PENALTY_X = 1e0;
    settings.T_HORIZON = 1.0;
    settings.N_HORIZON = 25;

    settings.L2_COST_V = 1e-1;
    settings.L2_COST_P = 1e0;
    settings.L2_COST_A = 1e-3;
    settings.WITH_X_BOUNDS = true;
end
