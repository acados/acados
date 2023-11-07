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


function generate_c_code_nonlinear_least_squares( model, opts, target_dir )

%% import casadi
import casadi.*

casadi_opts = struct('mex', false, 'casadi_int', 'int', 'casadi_real', 'double');
check_casadi_version();

%% load model
% x
x = model.sym_x;
nx = length(x);
% check type
if isa(x(1), 'casadi.SX')
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
z = model.sym_z;

model_name = model.name;

% cd to target folder
if nargin > 2
    original_dir = pwd;
    if ~exist(target_dir, 'dir')
        mkdir(target_dir);
    end
    chdir(target_dir)
end

if isfield(model, 'cost_expr_y_0')
    fun = model.cost_expr_y_0;
    if all(size(fun) == 0)
        error('empty cost_expr_y_0 is not allowed.\nPlease use SX.zeros(1) or %s',...
              'linear least squares formulation for a zero cost term.');
    end
    % generate jacobians
    jac_x = jacobian(fun, x);
    jac_u = jacobian(fun, u);
    % output symbolics
    ny_0 = length(fun);
    if isSX
        y_0 = SX.sym('y', ny_0, 1);
    else
        y_0 = MX.sym('y', ny_0, 1);
    end
    % generate hessian
    y_0_adj = jtimes(fun, x, y_0, true);
    y_0_hess = jacobian(y_0_adj, [u; x]);
    dy_dz = jacobian(fun, z)
    % Set up functions
    y_0_fun = Function([model_name,'_cost_y_0_fun'], {x, u, z, p}, {fun});
    y_0_fun_jac_ut_xt = Function([model_name,'_cost_y_0_fun_jac_ut_xt'], {x, u, z, p}, {fun, [jac_u'; jac_x'], dy_dz});
    y_0_hess = Function([model_name,'_cost_y_0_hess'], {x, u, z, y_0, p}, {y_0_hess});
    % generate C code
    y_0_fun.generate([model_name,'_cost_y_0_fun'], casadi_opts);
    y_0_fun_jac_ut_xt.generate([model_name,'_cost_y_0_fun_jac_ut_xt'], casadi_opts);
    y_0_hess.generate([model_name,'_cost_y_0_hess'], casadi_opts);
end

if isfield(model, 'cost_expr_y')
    fun = model.cost_expr_y;
    if all(size(fun) == 0)
        error('empty cost_expr_y is not allowed.\nPlease use SX.zeros(1) or %s',...
              'linear least squares formulation for a zero cost term.');
    end
    % generate jacobians
    jac_x = jacobian(fun, x);
    jac_u = jacobian(fun, u);
    % output symbolics
    ny = length(fun);
    if isSX
        y = SX.sym('y', ny, 1);
    else
        y = MX.sym('y', ny, 1);
    end
    % generate hessian
    y_adj = jtimes(fun, [u; x], y, true);
    y_hess = jacobian(y_adj, [u; x]);
    dy_dz = jacobian(fun, z);
    % Set up functions
    y_fun = Function([model_name,'_cost_y_fun'], {x, u, z, p}, {fun});
    y_fun_jac_ut_xt = Function([model_name,'_cost_y_fun_jac_ut_xt'],...
                              {x, u, z, p}, {fun, [jac_u'; jac_x'], dy_dz});
    y_hess = Function([model_name,'_cost_y_hess'], {x, u, z, y, p}, {y_hess});
    % generate C code
    y_fun.generate([model_name,'_cost_y_fun'], casadi_opts);
    y_fun_jac_ut_xt.generate([model_name,'_cost_y_fun_jac_ut_xt'], casadi_opts);
    y_hess.generate([model_name,'_cost_y_hess'], casadi_opts);
end

if isfield(model, 'cost_expr_y_e')
    fun = model.cost_expr_y_e;
    if all(size(fun) == 0)
        error('empty cost_expr_y_e is not allowed.\nPlease use SX.zeros(1) or %s',...
              'linear least squares formulation for a zero cost term.');
    end
    % generate jacobians
    jac_x = jacobian(fun, x);
    dy_dz = jacobian(fun, z);
    % output symbolics
    ny_e = length(fun);
    if isSX
        y_e = SX.sym('y', ny_e, 1);
        u = SX.sym('u', 0, 0)
    else
        y_e = MX.sym('y', ny_e, 1);
        u = MX.sym('u', 0, 0)
    end
    % generate hessian
    y_e_adj = jtimes(fun, x, y_e, true);
    y_e_hess = jacobian(y_e_adj, x);
    % Set up functions
    y_e_fun = Function([model_name,'_cost_y_e_fun'], {x, u, z, p}, {fun});
    y_e_fun_jac_ut_xt = Function([model_name,'_cost_y_e_fun_jac_ut_xt'], {x, u, z, p}, {fun, jac_x', dy_dz});
    y_e_hess = Function([model_name,'_cost_y_e_hess'], {x, u, z, y_e, p}, {y_e_hess});
    % generate C code
    y_e_fun.generate([model_name,'_cost_y_e_fun'], casadi_opts);
    y_e_fun_jac_ut_xt.generate([model_name,'_cost_y_e_fun_jac_ut_xt'], casadi_opts);
    y_e_hess.generate([model_name,'_cost_y_e_hess'], casadi_opts);
end

if nargin > 2
    chdir(original_dir)
end

end

