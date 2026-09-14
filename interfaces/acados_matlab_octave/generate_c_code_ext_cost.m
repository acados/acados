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
        ext_cost = model.cost_expr_ext_cost_0;
        custom_hess = model.cost_expr_ext_cost_custom_hess_0;
        suffix_name = '_cost_ext_cost_0';
        diff_vars = vertcat(u, x, z);
        function_inputs = {x, u, z, p};
        error_message = 'Field `cost_expr_ext_cost_0` is required for cost_type_0 == EXTERNAL.';
    elseif strcmp(stage_type, "path")
        ext_cost = model.cost_expr_ext_cost;
        custom_hess = model.cost_expr_ext_cost_custom_hess;
        suffix_name = '_cost_ext_cost';
        diff_vars = vertcat(u, x, z);
        function_inputs = {x, u, z, p};
        error_message = 'Field `cost_expr_ext_cost` is required for cost_type == EXTERNAL.';
    elseif strcmp(stage_type, "terminal")
        ext_cost = model.cost_expr_ext_cost_e;
        custom_hess = model.cost_expr_ext_cost_custom_hess_e;
        suffix_name = '_cost_ext_cost_e';
        diff_vars = x;
        function_inputs = {x, p};
        error_message = 'Field `cost_expr_ext_cost_e` is required for cost_type_e == EXTERNAL.';
    else
        error("Unknown stage type.")
    end

    if isempty(ext_cost)
        error(error_message)
    end

    if strcmp(stage_type, "terminal")
        if any(which_depends(ext_cost, model.u))
            error('terminal cost cannot depend on u.');
        end
        if any(which_depends(ext_cost, model.z))
            error('terminal cost cannot depend on z.');
        end
    end

    [full_hess, grad] = hessian(ext_cost, diff_vars);

    context.add_function_definition([model.name suffix_name '_fun'], ...
        function_inputs, {ext_cost}, target_dir, 'cost');
    context.add_function_definition([model.name suffix_name '_fun_jac'], ...
        function_inputs, {ext_cost, grad}, target_dir, 'cost');

    if isempty(custom_hess)
        custom_hess = full_hess;
    end

    % lower triangular is sufficient
    custom_hess = ca.tril(custom_hess)

    context.add_function_definition([model.name suffix_name '_fun_jac_hess'], ...
        function_inputs, {ext_cost, grad, custom_hess}, target_dir, 'cost');

end

