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


function generate_c_code_explicit_ode(context, model, model_dir)

    import casadi.*

    %% load model
    x = model.x;
    u = model.u;
    p = model.p;
    nx = length(x);
    nu = length(u);

    % check type
    if isa(x(1), 'casadi.SX')
        isSX = true;
    else
        isSX = false;
    end

    if isempty(model.f_expl_expr)
        error("Field `f_expl_expr` is required for integrator type ERK.")
    end

    f_expl = model.f_expl_expr;

    % setup expressions
    if isSX
        Sx = SX.sym('Sx', nx, nx);
        Su = SX.sym('Su', nx, nu);
        lambdaX = SX.sym('lambdaX', nx, 1);
        vdeX = SX.zeros(nx, nx);
        vdeU = SX.zeros(nx, nu) + jacobian(f_expl, u);
    else
        Sx = MX.sym('Sx', nx, nx);
        Su = MX.sym('Su', nx, nu);
        lambdaX = MX.sym('lambdaX', nx, 1);
        vdeX = MX.zeros(nx, nx);
        vdeU = MX.zeros(nx, nu) + jacobian(f_expl, u);
    end

    vdeX = vdeX + jtimes(f_expl, x, Sx);
    vdeU = vdeU + jtimes(f_expl, x, Su);

    % 'true' at the end tells to transpose the jacobian before multiplication => reverse mode
    adj = jtimes(f_expl, [x;u], lambdaX, true);

    if context.opts.generate_hess
        S_forw = vertcat(horzcat(Sx, Su), horzcat(zeros(nu,nx), eye(nu)));
        hess = S_forw.'*jtimes(adj, [x;u], S_forw);
        % TODO uncompress it ?????
        hess2 = [];
        for j = 1:nx+nu
            for i = j:nx+nu
                hess2 = [hess2; hess(i,j)];
            end
        end
    end

    fun_name = [model.name,'_expl_ode_fun'];
    context.add_function_definition(fun_name, {x, u, p}, {f_expl}, model_dir, 'dyn');

    fun_name = [model.name,'_expl_vde_forw'];
    context.add_function_definition(fun_name, {x, Sx, Su, u, p}, {f_expl, vdeX, vdeU}, model_dir, 'dyn');

    fun_name = [model.name,'_expl_vde_adj'];
    context.add_function_definition(fun_name, {x, lambdaX, u, p}, {adj}, model_dir, 'dyn');

    if context.opts.generate_hess
        fun_name = [model.name,'_expl_ode_hess'];
        context.add_function_definition(fun_name, {x, Sx, Su, lambdaX, u, p}, {adj, hess2}, model_dir, 'dyn');
    end

end
