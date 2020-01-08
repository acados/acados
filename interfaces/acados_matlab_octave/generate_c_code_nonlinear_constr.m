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


function generate_c_code_nonlinear_constr( model, opts, target_dir )

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
if isa(x(1), 'casadi.SX')
    isSX = true;
else
    isSX = false;
end
% u
u = model.sym_u;
nu = length(u);
% z
if isfield(model, 'sym_z')
    z = model.sym_z;
    nz = length(z);
else
    if isSX
        z = SX.sym('z',0, 0);
    else
        z = MX.sym('z',0, 0);
    end
    nz = 0;
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

model_name = model.name;

% cd to target folder
if nargin > 2
    original_dir = pwd;
    if ~exist(target_dir, 'dir')
        mkdir(target_dir);
    end
    chdir(target_dir)
end

if isfield(model, 'constr_expr_h')
    h = model.constr_expr_h;
    % multipliers for hessian
    nh = length(h);
    if isSX
        lam_h = SX.sym('lam_h', nh, 1);
    else
        lam_h = MX.sym('lam_h', nh, 1);
    end
    % generate jacobians
    jac_ux = jacobian(h, [u; x]);
    jac_z  = jacobian(h, [z]);
    % generate adjoint
    adj_ux = jtimes(h, [u; x], lam_h, true);
    % generate hessian
    hess_ux = jacobian(adj_ux, [u; x]);
    % Set up functions
    h_fun = Function([model_name,'_constr_h_fun'], {x, u, z, p}, {h});
    h_fun_jac_uxt_zt = Function([model_name,'_constr_h_fun_jac_uxt_zt'], {x, u, z, p}, {h, jac_ux', jac_z'});
    h_fun_jac_ut_xt_hess = Function([model_name,'_constr_h_fun_jac_ut_xt_hess'], {x, u, lam_h, z, p}, {h, jac_ux', hess_ux});
    % generate C code
    h_fun.generate([model_name,'_constr_h_fun'], casadi_opts);
    h_fun_jac_uxt_zt.generate([model_name,'_constr_h_fun_jac_uxt_zt'], casadi_opts);
    h_fun_jac_ut_xt_hess.generate([model_name,'_constr_h_fun_jac_ut_xt_hess'], casadi_opts);
end

if isfield(model, 'constr_expr_h_e')
    h_e = model.constr_expr_h_e;
    % multipliers for hessian
    nh_e = length(h_e);
    if isSX
        lam_h_e = SX.sym('lam_h', nh_e, 1);
    else
        lam_h_e = MX.sym('lam_h', nh_e, 1);
    end
    % generate jacobians
    jac_x_e     = jacobian(h_e, x);
    % generate adjoint (TODO output also adjoint when hessian is computed ?????)
    adj_ux_e = jtimes(h_e, x, lam_h_e, true);
    % generate hessian
    hess_ux_e = jacobian(adj_ux_e, x);
    % Set up functions
    h_e_fun = Function([model_name,'_constr_h_e_fun'], {x, p}, {h_e});
    h_e_fun_jac_uxt_zt = Function([model_name,'_constr_h_e_fun_jac_uxt_zt'], {x, p}, {h_e, jac_x_e'});
    h_e_fun_jac_uxt_zt_hess = Function([model_name,'_constr_h_e_fun_jac_uxt_zt_hess'], {x, lam_h_e, p}, {h_e, jac_x_e', hess_ux_e});
    % generate C code
    h_e_fun.generate([model_name,'_constr_h_e_fun'], casadi_opts);
    h_e_fun_jac_uxt_zt.generate([model_name,'_constr_h_e_fun_jac_uxt_zt'], casadi_opts);
    h_e_fun_jac_uxt_zt_hess.generate([model_name,'_constr_h_e_fun_jac_uxt_zt_hess'], casadi_opts);
end

if nargin > 2
    chdir(original_dir)
end


end