%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.

%


function generate_c_code_nonlinear_least_squares(context, model, target_dir, stage_type)

    import casadi.*

    %% load model
    x = model.x;
    u = model.u;
    z = model.z;
    p = model.p;

    % check type
    if isa(x(1), 'casadi.SX')
        isSX = true;
    else
        isSX = false;
    end

    if strcmp(stage_type, 'initial')

        if ~isempty(model.cost_y_expr_0)
            fun = model.cost_y_expr_0;
        elseif isempty(model.cost_y_expr_0) && ~isempty(model.cost_y_expr)
            disp('path used')
            fun = model.cost_y_expr;
        else
            error('empty cost_y_expr is not allowed.\nPlease use SX.zeros(1) or %s',...
            'linear least squares formulation for a zero cost term.');
        end
        % generate jacobians
        jac_x = jacobian(fun, x);
        jac_u = jacobian(fun, u);
        % output symbolics
        ny_0 = length(fun);
        if isSX
            y_0 = SX.sym('y', ny_0, 1);
        else
            y_0 = MX.sym('y', ny_0, 1);
        end
        % generate hessian
        y_0_adj = jtimes(fun, [u; x], y_0, true);
        y_0_hess = jacobian(y_0_adj, [u; x], struct('symmetric', isSX));
        dy_dz = jacobian(fun, z);
        % add functions to context
        context.add_function_definition([model.name,'_cost_y_0_fun'], {x, u, z, p}, {fun}, target_dir, 'cost');
        context.add_function_definition([model.name,'_cost_y_0_fun_jac_ut_xt'], {x, u, z, p}, {fun, [jac_u'; jac_x'], dy_dz}, target_dir, 'cost');
        if context.opts.generate_hess
            context.add_function_definition([model.name,'_cost_y_0_hess'], {x, u, z, y_0, p}, {y_0_hess}, target_dir, 'cost');
        end
    elseif strcmp(stage_type, 'path')
        fun = model.cost_y_expr;
        if isempty(fun)
            error('empty cost_y_expr is not allowed.\nPlease use SX.zeros(1) or %s',...
                'linear least squares formulation for a zero cost term.');
        end
        % generate jacobians
        jac_x = jacobian(fun, x);
        jac_u = jacobian(fun, u);
        % output symbolics
        ny = length(fun);
        if isSX
            y = SX.sym('y', ny, 1);
        else
            y = MX.sym('y', ny, 1);
        end
        % generate hessian
        y_adj = jtimes(fun, [u; x], y, true);
        y_hess = jacobian(y_adj, [u; x], struct('symmetric', isSX));
        dy_dz = jacobian(fun, z);
        % add functions to context
        context.add_function_definition([model.name,'_cost_y_fun'], {x, u, z, p}, {fun}, target_dir, 'cost');
        context.add_function_definition([model.name,'_cost_y_fun_jac_ut_xt'], ...
                                {x, u, z, p}, {fun, [jac_u'; jac_x'], dy_dz}, target_dir, 'cost');

        if context.opts.generate_hess
            context.add_function_definition([model.name,'_cost_y_hess'], {x, u, z, y, p}, {y_hess}, target_dir, 'cost');
        end

    elseif strcmp(stage_type, 'terminal')
        fun = model.cost_y_expr_e;
        if isempty(fun)
            error('empty cost_y_expr_e is not allowed.\nPlease use SX.zeros(1) or %s',...
                'linear least squares formulation for a zero cost term.');
        end
        % generate jacobians
        jac_x = jacobian(fun, x);
        dy_dz = jacobian(fun, z);
        % output symbolics
        ny_e = length(fun);
        if any(which_depends(fun, model.u))
            error('terminal cost cannot depend on u.');
        end
        if any(which_depends(fun, model.z))
            error('terminal cost cannot depend on z.');
        end
        if isSX
            y_e = SX.sym('y', ny_e, 1);
            u = SX.sym('u', 0, 0);
        else
            y_e = MX.sym('y', ny_e, 1);
            u = MX.sym('u', 0, 0);
        end
        % generate hessian
        y_e_adj = jtimes(fun, x, y_e, true);
        y_e_hess = jacobian(y_e_adj, x, struct('symmetric', isSX));
        % add functions to context
        context.add_function_definition([model.name,'_cost_y_e_fun'], {x, u, z, p}, {fun}, target_dir, 'cost');
        context.add_function_definition([model.name,'_cost_y_e_fun_jac_ut_xt'], {x, u, z, p}, {fun, jac_x', dy_dz}, target_dir, 'cost');

        if context.opts.generate_hess
            context.add_function_definition([model.name,'_cost_y_e_hess'], {x, u, z, y_e, p}, {y_e_hess}, target_dir, 'cost');
        end
    end

end

