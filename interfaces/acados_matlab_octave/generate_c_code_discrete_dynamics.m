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
    end

    context.add_function_definition([model.name,'_dyn_disc_phi_fun'], {x, u, p}, {phi}, model_dir, 'dyn');
    context.add_function_definition([model.name,'_dyn_disc_phi_fun_jac'], {x, u, p}, {phi, jac_ux'}, model_dir, 'dyn');
    if context.opts.generate_hess
        context.add_function_definition([model.name,'_dyn_disc_phi_fun_jac_hess'], {x, u, pi, p}, {phi, jac_ux', hess_ux}, model_dir, 'dyn');
    end

end