%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% The 2-Clause BSD License
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.;

%

function generate_c_code_conl_cost(context, model, target_dir, stage_type)

    import casadi.*

    %% load model
    x = model.x;
    u = model.u;
    z = model.z;
    p = model.p;
    t = model.t;

    % check type
    if isa(x(1), 'casadi.SX')
        isSX = true;
    else
        isSX = false;
    end

    if strcmp(stage_type, 'initial')
        yref = model.cost_r_in_psi_expr_0;
        y_expr = model.cost_y_expr_0;
        outer_expr = model.cost_psi_expr_0;
        res_expr = model.cost_r_in_psi_expr_0;
        custom_hess = model.cost_conl_custom_outer_hess_0;

        suffix_name_fun = '_conl_cost_0_fun';
        suffix_name_fun_jac_hess = '_conl_cost_0_fun_jac_hess';

    elseif strcmp(stage_type, 'path')
        yref = model.cost_r_in_psi_expr;
        y_expr = model.cost_y_expr;
        outer_expr = model.cost_psi_expr;
        res_expr = model.cost_r_in_psi_expr;
        custom_hess = model.cost_conl_custom_outer_hess;

        suffix_name_fun = '_conl_cost_fun';
        suffix_name_fun_jac_hess = '_conl_cost_fun_jac_hess';

    elseif strcmp(stage_type, 'terminal')
        yref = model.cost_r_in_psi_expr_e;
        y_expr = model.cost_y_expr_e;
        outer_expr = model.cost_psi_expr_e;
        res_expr = model.cost_r_in_psi_expr_e;
        custom_hess = model.cost_conl_custom_outer_hess_e;

        suffix_name_fun = '_conl_cost_e_fun';
        suffix_name_fun_jac_hess = '_conl_cost_e_fun_jac_hess';

        % create dummy u, z for terminal stage
        if isSX
            u = SX.sym('u', 0, 0);
            z = SX.sym('z', 0, 0);
        else
            u = MX.sym('u', 0, 0);
            z = MX.sym('z', 0, 0);
        end
    end

    % Check if required expressions are defined
    if isempty(y_expr)
        error(['cost_y_expr for stage ' stage_type ' is empty. Required for CONVEX_OVER_NONLINEAR cost.']);
    end
    if isempty(outer_expr)
        error(['cost_psi_expr for stage ' stage_type ' is empty. Required for CONVEX_OVER_NONLINEAR cost.']);
    end
    if isempty(res_expr)
        error(['cost_r_in_psi_expr for stage ' stage_type ' is empty. Required for CONVEX_OVER_NONLINEAR cost.']);
    end

    % Set up function names
    fun_name_cost_fun = [model.name suffix_name_fun];
    fun_name_cost_fun_jac_hess = [model.name suffix_name_fun_jac_hess];

    % Set up functions to be exported
    % inner_expr = y_expr - yref
    inner_expr = y_expr - yref;

    % outer_loss_fun: psi(residual, t, p) = outer_expr
    outer_loss_fun = Function('psi', {res_expr, t, p}, {outer_expr});

    % cost_expr = outer_loss_fun(inner_expr, t, p)
    cost_expr = outer_loss_fun(inner_expr, t, p);

    % outer_loss_grad_fun: gradient of outer loss w.r.t. residual
    outer_grad = jacobian(outer_expr, res_expr).';
    outer_loss_grad_fun = Function('outer_loss_grad', {res_expr, t, p}, {outer_grad});

    % Compute or use custom hessian
    if isempty(custom_hess)
        % Compute hessian of outer loss w.r.t. residual
        [hess, ~] = hessian(outer_loss_fun(res_expr, t, p), res_expr);
    else
        hess = custom_hess;
    end

    outer_hess_fun = Function('outer_hess', {res_expr, t, p}, {hess});
    outer_hess_expr = outer_hess_fun(inner_expr, t, p);

    % Check if hessian is diagonal
    outer_hess_is_diag = outer_hess_expr.sparsity().is_diag();

    % if residual dimension <= 4, do not exploit diagonal structure
    ny = length(res_expr);
    if ny <= 4
        outer_hess_is_diag = 0;
    end

    % Jacobians of inner expression w.r.t. u, x, z
    Jt_ux_expr = jacobian(inner_expr, vertcat(u, x)).';
    Jt_z_expr = jacobian(inner_expr, z).';

    % Add functions to context
    context.add_function_definition(fun_name_cost_fun, ...
        {x, u, z, yref, t, p}, {cost_expr}, target_dir, 'cost');

    context.add_function_definition(fun_name_cost_fun_jac_hess, ...
        {x, u, z, yref, t, p}, ...
        {cost_expr, outer_loss_grad_fun(inner_expr, t, p), Jt_ux_expr, Jt_z_expr, outer_hess_expr, outer_hess_is_diag}, ...
        target_dir, 'cost');

end

