%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.

%


function generate_c_code_explicit_ode_with_cost_state(context, model, model_dir)

    import casadi.*

    % check type
    if isa(model.x, 'casadi.SX')
        cost_state = SX.sym('cost_state');
    else
        cost_state = MX.sym('cost_state');
    end
    x_with_cost = [model.x; cost_state];

    add_explicit_ode_function_definitions(context, model, model_dir, x_with_cost, model.f_expl_expr_with_cost);
end

