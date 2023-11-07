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


function generate_c_code_ext_cost( model, opts, target_dir )

%% import casadi
import casadi.*

casadi_opts = struct('mex', false, 'casadi_int', 'int', 'casadi_real', 'double');
check_casadi_version();

% cd to target folder
if nargin > 2
    original_dir = pwd;
    if ~exist(target_dir, 'dir')
        mkdir(target_dir);
    end
    chdir(target_dir)
end

%% load model
% x
x = model.sym_x;
% check type
if isa(x(1), 'casadi.SX')
    isSX = true;
else
    isSX = false;
end
% u
u = model.sym_u;
% z
z = model.sym_z;
% p
if isfield(model, 'sym_p')
    p = model.sym_p;
else
    if isSX
        p = SX.sym('p',0, 0);
    else
        p = MX.sym('p',0, 0);
    end
end

model_name = model.name;

if isfield(model, 'cost_expr_ext_cost') && strcmp(model.cost_ext_fun_type, 'casadi') && strcmp(model.cost_type, 'ext_cost')
    ext_cost = model.cost_expr_ext_cost;
    % generate jacobian, hessian
    [full_hess, grad] = hessian(ext_cost, vertcat(u, x, z));
    % Set up functions
    ext_cost_fun = Function([model_name,'_cost_ext_cost_fun'], {x, u, z, p}, {ext_cost});
    ext_cost_fun_jac = Function([model_name,'_cost_ext_cost_fun_jac'], {x, u, z, p}, {ext_cost, grad});
    if isfield(model, 'cost_expr_ext_cost_custom_hess')
        ext_cost_fun_jac_hess = Function([model_name,'_cost_ext_cost_fun_jac_hess'], {x, u, z, p},...
                                     {ext_cost, grad, model.cost_expr_ext_cost_custom_hess});
    else
        ext_cost_fun_jac_hess = Function([model_name,'_cost_ext_cost_fun_jac_hess'], {x, u, z, p},...
                                     {ext_cost, grad, full_hess});
    end
    % generate C code
    ext_cost_fun.generate([model_name,'_cost_ext_cost_fun'], casadi_opts);
    ext_cost_fun_jac_hess.generate([model_name,'_cost_ext_cost_fun_jac_hess'], casadi_opts);
    ext_cost_fun_jac.generate([model_name,'_cost_ext_cost_fun_jac'], casadi_opts);
end

if isfield(model, 'cost_expr_ext_cost_0') && strcmp(model.cost_ext_fun_type_0, 'casadi') && strcmp(model.cost_type_0, 'ext_cost')
    ext_cost_0 = model.cost_expr_ext_cost_0;
    % generate jacobian, hessian
    [full_hess, grad] = hessian(ext_cost_0, vertcat(u, x, z));
    % Set up functions
    ext_cost_0_fun = Function([model_name,'_cost_ext_cost_0_fun'], {x, u, z, p}, {ext_cost_0});
    ext_cost_0_fun_jac = Function([model_name,'_cost_ext_cost_0_fun_jac'], {x, u, z, p}, {ext_cost_0, grad});
    if isfield(model, 'cost_expr_ext_cost_custom_hess_0')
        ext_cost_0_fun_jac_hess = Function([model_name,'_cost_ext_cost_0_fun_jac_hess'], {x, u, z, p},...
                                     {ext_cost_0, grad, model.cost_expr_ext_cost_custom_hess_0});
    else
        ext_cost_0_fun_jac_hess = Function([model_name,'_cost_ext_cost_0_fun_jac_hess'], {x, u, z, p}, {ext_cost_0, grad, full_hess});
    end
    % generate C code
    ext_cost_0_fun.generate([model_name,'_cost_ext_cost_0_fun'], casadi_opts);
    ext_cost_0_fun_jac.generate([model_name,'_cost_ext_cost_0_fun_jac'], casadi_opts);
    ext_cost_0_fun_jac_hess.generate([model_name,'_cost_ext_cost_0_fun_jac_hess'], casadi_opts);
end

if isfield(model, 'cost_expr_ext_cost_e') && strcmp(model.cost_ext_fun_type_e, 'casadi') && strcmp(model.cost_type_e, 'ext_cost')
    ext_cost_e = model.cost_expr_ext_cost_e;
    % generate jacobians
    jac_x_e = jacobian(ext_cost_e, x);
    % generate hessians
    hes_xx_e = jacobian(jac_x_e', x);
    % Set up functions
    ext_cost_e_fun = Function([model_name,'_cost_ext_cost_e_fun'], {x, p}, {ext_cost_e});
    ext_cost_e_fun_jac = Function([model_name,'_cost_ext_cost_e_fun_jac'], {x, p}, {ext_cost_e, jac_x_e'});
    if isfield(model, 'cost_expr_ext_cost_custom_hess_e')
        ext_cost_e_fun_jac_hess = Function([model_name,'_cost_ext_cost_e_fun_jac_hess'], {x, p},...
                                     {ext_cost_e, jac_x_e', model.cost_expr_ext_cost_custom_hess_e});
    else
        ext_cost_e_fun_jac_hess = Function([model_name,'_cost_ext_cost_e_fun_jac_hess'], {x, p}, {ext_cost_e, jac_x_e', hes_xx_e});
    end
    % generate C code
    ext_cost_e_fun.generate([model_name,'_cost_ext_cost_e_fun'], casadi_opts);
    ext_cost_e_fun_jac.generate([model_name,'_cost_ext_cost_e_fun_jac'], casadi_opts);
    ext_cost_e_fun_jac_hess.generate([model_name,'_cost_ext_cost_e_fun_jac_hess'], casadi_opts);
end

if nargin > 2
    chdir(original_dir)
end

end

