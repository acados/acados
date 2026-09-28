%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.



function [p_global, m, l, coefficients, coefficient_vals, knots, p_global_values] = create_p_global(lut)

    import casadi.*
    m = MX.sym('m');
    l = MX.sym('l');
    p_global = {m, l};
    p_global_values = [0.1; 0.8];

    large_scale = false;
    if lut
        if large_scale
            % large scale lookup table
            knots = {0:2^8-1,0:2^8-1};
            coefficient_vals = 0.1*ones((2^8-3)^2, 1);
        else
            % small scale lookup table
            knots = {0:2^4,0:2^4};
            coefficient_vals = 0.1*ones((2^4-3)^2, 1);
        end

        coefficients = MX.sym('coefficient', numel(coefficient_vals), 1);
        p_global{end+1} = coefficients;
        p_global_values = [p_global_values; coefficient_vals(:)];
    else
        coefficient_vals = [];
        knots = [];
        coefficients = MX.sym('coefficient', 0, 1);
    end

    p_global = vertcat(p_global{:});
end
