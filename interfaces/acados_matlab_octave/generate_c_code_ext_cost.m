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



function generate_c_code_ext_cost(context, model, target_dir, stage_type)

    import casadi.*

    %% load model
    x = model.x;
    u = model.u;
    z = model.z;
    p = model.p;

    if strcmp(stage_type, "initial")
        if isempty(model.cost_expr_ext_cost_0)
            error('Field `cost_expr_ext_cost_0` is required for cost_type_0 == EXTERNAL.')
        end

        ext_cost_0 = model.cost_expr_ext_cost_0;
        % generate jacobian, hessian
        [full_hess, grad] = hessian(ext_cost_0, vertcat(u, x, z));
        % add functions to context
        context.add_function_definition([model.name,'_cost_ext_cost_0_fun'], {x, u, z, p}, {ext_cost_0}, target_dir, 'cost');
        context.add_function_definition([model.name,'_cost_ext_cost_0_fun_jac'], {x, u, z, p}, {ext_cost_0, grad}, target_dir, 'cost');
        if ~isempty(model.cost_expr_ext_cost_custom_hess_0)
            context.add_function_definition([model.name,'_cost_ext_cost_0_fun_jac_hess'], {x, u, z, p},...
                                        {ext_cost_0, grad, model.cost_expr_ext_cost_custom_hess_0}, target_dir, 'cost');
        else
            context.add_function_definition([model.name,'_cost_ext_cost_0_fun_jac_hess'], {x, u, z, p}, {ext_cost_0, grad, full_hess}, target_dir, 'cost');
        end

    elseif strcmp(stage_type, "path")
        if isempty(model.cost_expr_ext_cost)
            error('Field `cost_expr_ext_cost` is required for cost_type == EXTERNAL.')
        end
        ext_cost = model.cost_expr_ext_cost;
        % generate jacobian, hessian
        [full_hess, grad] = hessian(ext_cost, vertcat(u, x, z));
        % add functions to context
        context.add_function_definition([model.name,'_cost_ext_cost_fun'], {x, u, z, p}, {ext_cost}, target_dir, 'cost');
        context.add_function_definition([model.name,'_cost_ext_cost_fun_jac'], {x, u, z, p}, {ext_cost, grad}, target_dir, 'cost');
        if ~isempty(model.cost_expr_ext_cost_custom_hess)
            context.add_function_definition([model.name,'_cost_ext_cost_fun_jac_hess'], {x, u, z, p}, ...
                                        {ext_cost, grad, model.cost_expr_ext_cost_custom_hess}, target_dir, 'cost');
        else
            context.add_function_definition([model.name,'_cost_ext_cost_fun_jac_hess'], {x, u, z, p}, ...
                                        {ext_cost, grad, full_hess}, target_dir, 'cost');
        end

    elseif strcmp(stage_type, "terminal")
        if isempty(model.cost_expr_ext_cost_e)
            error('Field `cost_expr_ext_cost_e` is required for cost_type_e == EXTERNAL.')
        end
        ext_cost_e = model.cost_expr_ext_cost_e;
        if any(which_depends(ext_cost_e, model.u))
            error('terminal cost cannot depend on u.');
        end
        if any(which_depends(ext_cost_e, model.z))
            error('terminal cost cannot depend on z.');
        end
        % generate jacobians
        jac_x_e = jacobian(ext_cost_e, x);
        % generate hessians
        hes_xx_e = jacobian(jac_x_e', x);
        % add functions to context
        context.add_function_definition([model.name,'_cost_ext_cost_e_fun'], {x, p}, {ext_cost_e}, target_dir, 'cost');
        context.add_function_definition([model.name,'_cost_ext_cost_e_fun_jac'], {x, p}, {ext_cost_e, jac_x_e'}, target_dir, 'cost');
        if ~isempty(model.cost_expr_ext_cost_custom_hess_e)
            context.add_function_definition([model.name,'_cost_ext_cost_e_fun_jac_hess'], {x, p},...
                                        {ext_cost_e, jac_x_e', model.cost_expr_ext_cost_custom_hess_e}, target_dir, 'cost');
        else
            context.add_function_definition([model.name, '_cost_ext_cost_e_fun_jac_hess'], {x, p}, {ext_cost_e, jac_x_e', hes_xx_e}, target_dir, 'cost');
        end
    else
        error("Unknown stage type.")
    end

end

