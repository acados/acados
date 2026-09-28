%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.

%


% Shared implementation for the explicit ODE code generation: ERK and ERK_WITH_COST
function add_explicit_ode_function_definitions(context, model, model_dir, x, f_expl)

    import casadi.*

    %% load model
    u = model.u;
    p = model.p;
    nx = length(x);
    nu = length(u);
    np = length(p);

    % check type
    if isa(x, 'casadi.SX')
        isSX = true;
    else
        isSX = false;
    end

    % setup expressions
    if isSX
        Sx = SX.sym('Sx', nx, nx);
        Su = SX.sym('Su', nx, nu);
        lambdaX = SX.sym('lambdaX', nx, 1);
        vdeX = SX.zeros(nx, nx);
        vdeU = SX.zeros(nx, nu) + jacobian(f_expl, u);
        if context.opts.sens_forw_p
            Sp   = SX.sym('Sp', nx, np);
            vdeP = SX.zeros(nx, np) + jacobian(f_expl, p);  % f_p
        end
    else
        Sx = MX.sym('Sx', nx, nx);
        Su = MX.sym('Su', nx, nu);
        lambdaX = MX.sym('lambdaX', nx, 1);
        vdeX = MX.zeros(nx, nx);
        vdeU = MX.zeros(nx, nu) + jacobian(f_expl, u);
        if context.opts.sens_forw_p
            Sp   = MX.sym('Sp', nx, np);
            vdeP = MX.zeros(nx, np) + jacobian(f_expl, p);  % f_p
        end
    end

    vdeX = vdeX + jtimes(f_expl, x, Sx);
    vdeU = vdeU + jtimes(f_expl, x, Su);

    if context.opts.sens_forw_p
        vdeP = vdeP + jtimes(f_expl, x, Sp);   % A*Sp + f_p
    end

    % 'true' at the end tells to transpose the jacobian before multiplication => reverse mode
    adj = jtimes(f_expl, [x;u], lambdaX, true);

    if context.opts.generate_hess
        S_forw = vertcat(horzcat(Sx, Su), horzcat(zeros(nu,nx), eye(nu)));
        hess = S_forw.'*jtimes(adj, [x;u], S_forw);
        % vectorized lower triangular Hessian
        hess_vec = hess(hess.sparsity().makeDense().get_lower()+1);
    end

    fun_name = [model.name,'_expl_ode_fun'];
    context.add_function_definition(fun_name, {x, u, p}, {f_expl}, model_dir, 'dyn');

    fun_name = [model.name,'_expl_vde_forw'];
    context.add_function_definition(fun_name, {x, Sx, Su, u, p}, {f_expl, vdeX, vdeU}, model_dir, 'dyn');

    fun_name = [model.name,'_expl_vde_adj'];
    context.add_function_definition(fun_name, {x, lambdaX, u, p}, {adj}, model_dir, 'dyn');

    if context.opts.generate_hess
        fun_name = [model.name,'_expl_ode_hess'];
        context.add_function_definition(fun_name, {x, Sx, Su, lambdaX, u, p}, {adj, hess_vec}, model_dir, 'dyn');
    end

    % param-direction forward VDE
    if context.opts.sens_forw_p
        fun_name = [model.name,'_expl_vde_forw_p'];
        context.add_function_definition(fun_name, {x, Sp, u, p}, {vdeP}, model_dir, 'dyn');
    end
end
