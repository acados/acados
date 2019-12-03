%
% Copyright 2019 Gianluca Frison, Dimitris Kouzoupis, Robin Verschueren,
% Andrea Zanelli, Niels van Duijkeren, Jonathan Frey, Tommaso Sartor,
% Branimir Novoselnik, Rien Quirynen, Rezart Qelibari, Dang Doan,
% Jonas Koenemann, Yutao Chen, Tobias Schöls, Jonas Schlagenhauf, Moritz Diehl
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


function generate_c_code_ext_cost( model, opts )

%% import casadi
import casadi.*

casadi_version = CasadiMeta.version();
if ( strcmp(casadi_version(1:3),'3.4') || strcmp(casadi_version(1:3),'3.5')) % require casadi 3.4.x
    casadi_opts = struct('mex', false, 'casadi_int', 'int', 'casadi_real', 'double');
else % old casadi versions
    error('Please provide CasADi version 3.4 or 3.5 to ensure compatibility with acados')
end

%% load model
% x
x = model.sym_x;
nx = length(x);
% check type
if class(x(1)) == 'casadi.SX'
    isSX = true;
else
    isSX = false;
end
% u
u = model.sym_u;
nu = length(u);
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

model_name = model.name;

if isfield(model, 'cost_expr_ext_cost')
    ext_cost = model.cost_expr_ext_cost;
    % generate jacobians
    jac_x = jacobian(ext_cost, x);
    jac_u = jacobian(ext_cost, u);
    % generate hessians
    hes_uu = jacobian(jac_u', u);
    hes_xu = jacobian(jac_u', x);
    hes_ux = jacobian(jac_x', u);
    hes_xx = jacobian(jac_x', x);
    % Set up functions
    ext_cost_fun = Function([model_name,'_cost_ext_cost_fun'], {x, u, p}, {ext_cost});
    ext_cost_fun_jac_hes = Function([model_name,'_cost_ext_cost_fun_jac_hess'], {x, u, p}, {ext_cost, [jac_u'; jac_x'], [hes_uu, hes_xu; hes_ux, hes_xx]});
    % generate C code
    ext_cost_fun.generate([model_name,'_cost_ext_cost_fun'], casadi_opts);
    ext_cost_fun_jac_hes.generate([model_name,'_cost_ext_cost_fun_jac_hess'], casadi_opts);
end

if isfield(model, 'cost_expr_ext_cost_e')
    ext_cost_e = model.cost_expr_ext_cost_e;
    % generate jacobians
    jac_x_e = jacobian(ext_cost_e, x);
    % generate hessians
    hes_xx_e = jacobian(jac_x', x);
    % Set up functions
    ext_cost_e_fun = Function([model_name,'_cost_ext_cost_e_fun'], {x, p}, {ext_cost_e});
    ext_cost_e_fun_jac_hes = Function([model_name,'_cost_ext_cost_e_fun_jac_hess'], {x, p}, {ext_cost_e, jac_x_e', hes_xx_e});
    % generate C code
    ext_cost_e_fun.generate([model_name,'_cost_ext_cost_e_fun'], casadi_opts);
    ext_cost_e_fun_jac_hes.generate([model_name,'_cost_ext_cost_e_fun_jac_hess'], casadi_opts);
end


