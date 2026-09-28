%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.

%


function generate_c_code_implicit_ode(context, model, model_dir)

    import casadi.*

    %% load model
    x = model.x;
    u = model.u;
    p = model.p;
    xdot = model.xdot;
    z = model.z;

    nx = length(x);
    nz = length(z);

    % check type
    if isa(x(1), 'casadi.SX')
        isSX = true;
    else
        isSX = false;
    end
    if isempty(model.f_impl_expr)
        error("Field `f_impl_expr` is required for integrator type IRK.")
    end

    f_impl = model.f_impl_expr;

    %% generate jacobians
    jac_x = jacobian(f_impl, x);
    jac_xdot = jacobian(f_impl, xdot);
    jac_u = jacobian(f_impl, u);
    jac_z = jacobian(f_impl, z);

    if context.opts.sens_forw_p
        jac_p = jacobian(f_impl, p);
    end

    %% generate hessian
    if context.opts.generate_hess
        x_xdot_z_u = [x; xdot; z; u];
        if isSX
            multiplier  = SX.sym('multiplier', nx + nz);
        else
            multiplier  = MX.sym('multiplier', nx + nz);
        end
        % hessian computed as forward over adjoint !!!
        ADJ = jtimes(f_impl, x_xdot_z_u, multiplier, true);
        HESS = jacobian(ADJ, x_xdot_z_u, struct('symmetric', isSX));
    end

    fun_name = [model.name,'_impl_dae_fun'];
    context.add_function_definition(fun_name, {x, xdot, u, z, p}, {f_impl}, model_dir, 'dyn');

    fun_name = [model.name,'_impl_dae_fun_jac_x_xdot_z'];
    context.add_function_definition(fun_name, {x, xdot, u, z, p}, {f_impl, jac_x, jac_xdot, jac_z}, model_dir, 'dyn');

    fun_name = [model.name,'_impl_dae_jac_x_xdot_u_z'];
    context.add_function_definition(fun_name, {x, xdot, u, z, p}, {jac_x, jac_xdot, jac_u, jac_z}, model_dir, 'dyn');

    fun_name = [model.name,'_impl_dae_fun_jac_x_xdot_u'];
    context.add_function_definition(fun_name, {x, xdot, u, z, p}, {f_impl, jac_x, jac_xdot, jac_u}, model_dir, 'dyn');

    if context.opts.generate_hess
        fun_name = [model.name,'_impl_dae_hess'];
        context.add_function_definition(fun_name, {x, xdot, u, z, multiplier, p}, {HESS}, model_dir, 'dyn');
    end

    if context.opts.sens_forw_p
        fun_name = [model.name,'_impl_dae_jac_p'];
        context.add_function_definition(fun_name, {x, xdot, u, z, p}, {jac_p}, model_dir, 'dyn');
    end
end