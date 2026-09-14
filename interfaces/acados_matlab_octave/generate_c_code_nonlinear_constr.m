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


function generate_c_code_nonlinear_constr(context, model, target_dir, stage_type)

    import casadi.*

    %% load model
    x = model.x;
    u = model.u;
    p = model.p;
    z = model.z;

    if isa(x(1), 'casadi.SX')
        isSX = true;
    else
        isSX = false;
    end

    if strcmp(stage_type, 'initial')
        h = model.con_h_expr_0;
        suffix_name = '_constr_h_0';
        function_inputs = {x, u, z, p};
        ux = [u; x];
        is_terminal = false;
    elseif strcmp(stage_type, 'path')
        h = model.con_h_expr;
        suffix_name = '_constr_h';
        function_inputs = {x, u, z, p};
        ux = [u; x];
        is_terminal = false;
    elseif strcmp(stage_type, 'terminal')
        % NOTE: terminal node has no u, z
        h = model.con_h_expr_e;
        suffix_name = '_constr_h_e';
        function_inputs = {x, p};
        ux = x;
        is_terminal = true;

        if any(which_depends(h, model.u))
            error('terminal constraints cannot depend on u.');
        end
        if any(which_depends(h, model.z))
            error('terminal constraints cannot depend on z.');
        end
    else
        error("Unknown stage type.")
    end

    % multipliers for hessian
    nh = length(h);
    if isSX
        lam_h = SX.sym('lam_h', nh, 1);
    else
        lam_h = MX.sym('lam_h', nh, 1);
    end

    % generate jacobians and Hessian with respect to u, x
    jac_ux = jacobian(h, ux);
    adj_ux = jtimes(h, ux, lam_h, true);
    % see https://github.com/casadi/casadi/issues/3703
    hess_ux = jacobian(adj_ux, ux, struct('symmetric', isSX));
    % lower triangular is sufficient
    hess_ux = ca.tril(hess_ux)

    context.add_function_definition([model.name suffix_name '_fun'], ...
        function_inputs, {h}, target_dir, 'constr');
    context.add_function_definition([model.name suffix_name '_fun_jac_uxt_zt'], ...
        function_inputs, {h, jac_ux'}, target_dir, 'constr');

    if context.opts.generate_hess
        if is_terminal
            context.add_function_definition([model.name suffix_name '_fun_jac_uxt_zt_hess'], ...
                {x, lam_h, p}, {h, jac_ux', hess_ux}, target_dir, 'constr');
        else
            jac_z = jacobian(h, z);
            adj_z = jtimes(h, z, lam_h, true);
            hess_z = jacobian(adj_z, z, struct('symmetric', isSX));
            context.add_function_definition([model.name suffix_name '_fun_jac_uxt_zt_hess'], ...
                {x, u, lam_h, z, p}, {h, jac_ux', hess_ux, jac_z', hess_z}, target_dir, 'constr');
        end
    end
end
