%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.




function generate_c_code_discrete_dynamics(context, model, model_dir)

    import casadi.*

    %% load model
    x = model.x;
    u = model.u;
    p = model.p;
    pi = model.pi;
    nx = length(x);

    if isempty(model.disc_dyn_expr)
        error('Field `disc_dyn_expr` is required for discrete dynamics.')
    end
    phi = model.disc_dyn_expr;
    nx1 = length(phi);

    % check type
    if isa(x(1), 'casadi.SX')
        isSX = true;
    else
        isSX = false;
    end

    ux = vertcat(u, x);

    % generate jacobians
    if isempty(model.disc_dyn_custom_jac_ux_expr)
        jac_ux = jacobian(phi, ux);
    else
        jac_ux = model.disc_dyn_custom_jac_ux_expr;
    end

    % generate adjoint
    adj_ux = jtimes(phi, ux, pi, true);
    % generate hessian
    if context.opts.generate_hess
        if isempty(model.disc_dyn_custom_hess_ux_expr)
            hess_ux = jacobian(adj_ux, ux, struct('symmetric', isSX));
        else
            hess_ux = model.disc_dyn_custom_hess_ux_expr;
        end
        % lower triangular is sufficient
        hess_ux = tril(hess_ux);
        context.add_function_definition([model.name,'_dyn_disc_phi_fun_jac_hess'], {x, u, pi, p}, {phi, jac_ux', hess_ux}, model_dir, 'dyn');
    end

    context.add_function_definition([model.name,'_dyn_disc_phi_fun'], {x, u, p}, {phi}, model_dir, 'dyn');
    context.add_function_definition([model.name,'_dyn_disc_phi_fun_jac'], {x, u, p}, {phi, jac_ux'}, model_dir, 'dyn');
end