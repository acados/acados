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


function generate_c_code_explicit_ode( model, opts )

%% import casadi
import casadi.*

casadi_opts = struct('mex', false, 'casadi_int', 'int', 'casadi_real', 'double');
check_casadi_version();

if nargin > 1
    if isfield(opts, 'sens_hess')
        generate_hess = opts.sens_hess;
    else
        generate_hess = 'false';
%        if opts.print_info
%        disp('generate_hess option was not set - default is false')
%        end
    end
else
    generate_hess = 'false';
end
generate_hess = 'true'; % TODO remove when not needed any more !!!


%% load model
% x
is_template = false;
if isa(model, 'acados_template_mex.acados_model_json')
    is_template = true;
    % names without sym
    x = model.x;
    nx = length(x);
    % check type
    if isa(x(1), 'casadi.SX')
        isSX = true;
    else
        isSX = false;
    end
    % u
    u = model.u;
    nu = length(u);
    % p
    p = model.p;
    np = length(p);


else
    x = model.sym_x;
    nx = length(x);
    % check type
    if isa(x(1), 'casadi.SX')
        isSX = true;
    else
        isSX = false;
    end
    % u
    if isfield(model, 'sym_u')
        u = model.sym_u;
        nu = length(u);
    else
        if isSX
            u = SX.sym('u',0, 0);
        else
            u = MX.sym('u',0, 0);
        end
        nu = 0;
    end
    % p
    if isfield(model, 'sym_p')
        p = model.sym_p;
        np = length(p);
    else
        if isSX
            p = SX.sym('p',0, 0);
        else
            p = MX.sym('p',0, 0);
        end
        np = 0;
    end
end

model_name = model.name;

if isfield(model, 'dyn_expr_f')
    f_expl = model.dyn_expr_f;
    model_name = [model_name, '_dyn'];
elseif isfield(model, 'expr_f')
    f_expl = model.expr_f;
else
    f_expl = model.f_expl_expr;
end



%% set up functions to be exported
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
expl_ode_fun = Function([model_name,'_expl_ode_fun'], {x, u, p}, {f_expl});

vdeX = vdeX + jtimes(f_expl, x, Sx);

vdeU = vdeU + jtimes(f_expl, x, Su);

expl_vde_for = Function([model_name,'_expl_vde_forw'], {x, Sx, Su, u, p}, {f_expl, vdeX, vdeU});

% 'true' at the end tells to transpose the jacobian before multiplication => reverse mode
adj = jtimes(f_expl, [x;u], lambdaX, true);

expl_vde_adj = Function([model_name,'_expl_vde_adj'], {x, lambdaX, u, p}, {adj});

S_forw = vertcat(horzcat(Sx, Su), horzcat(zeros(nu,nx), eye(nu)));
hess = S_forw.'*jtimes(adj, [x;u], S_forw);
% TODO uncompress it ?????
hess2 = [];
for j = 1:nx+nu
    for i = j:nx+nu
        hess2 = [hess2; hess(i,j)];
    end
end

if is_template
    return_dir = pwd;
    if ~exist( fullfile(pwd,'c_generated_code'), 'dir')
        mkdir('c_generated_code');
    end
    cd 'c_generated_code'
    model_dir = [model_name, '_model'];
    if ~exist(fullfile(pwd, model_dir), 'dir')
        mkdir(model_dir);
    end
    cd(model_dir)
end

expl_ode_hes = Function([model_name,'_expl_ode_hess'], {x, Sx, Su, lambdaX, u, p}, {adj, hess2});

%% generate C code
expl_ode_fun.generate([model_name,'_expl_ode_fun'], casadi_opts);
expl_vde_for.generate([model_name,'_expl_vde_forw'], casadi_opts);
expl_vde_adj.generate([model_name,'_expl_vde_adj'], casadi_opts);
if strcmp(generate_hess, 'true')
    expl_ode_hes.generate([model_name,'_expl_ode_hess'], casadi_opts);
end

if is_template
    cd(return_dir);
end

end
