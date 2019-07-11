%   This file is part of acados.
%
%   acados is free software; you can redistribute it and/or
%   modify it under the terms of the GNU Lesser General Public
%   License as published by the Free Software Foundation; either
%   version 3 of the License, or (at your option) any later version.
%
%   acados is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%   Lesser General Public License for more details.
%
%   You should have received a copy of the GNU Lesser General Public
%   License along with acados; if not, write to the Free Software Foundation,
%   Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
%
%   Author: Jonathan Frey: jonathanpaulfrey(at)gmail.com

function generate_c_code_gnsf(model)

%% import casadi
import casadi.*

casadi_version = CasadiMeta.version();
if strcmp(casadi_version(1:3),'3.4') % require casadi 3.4.x
	casadi_opts = struct('mex', false, 'casadi_int', 'int', 'casadi_real', 'double');
else % old casadi versions
	error('Please download and install CasADi version 3.4.x to ensure compatibility with acados')
end

%% import models
% model matrices
A  = model.dyn_gnsf_A;
B  = model.dyn_gnsf_B;
C  = model.dyn_gnsf_C;
E  = model.dyn_gnsf_E;
c  = model.dyn_gnsf_c;

L_x    = model.dyn_gnsf_L_x;
L_z    = model.dyn_gnsf_L_z;
L_xdot = model.dyn_gnsf_L_xdot;
L_u    = model.dyn_gnsf_L_u;

A_LO = model.dyn_gnsf_A_LO;
E_LO = model.dyn_gnsf_E_LO;
B_LO = model.dyn_gnsf_B_LO;

% state permutation vector: x_gnsf = dvecpe(x, ipiv)
ipiv_x = model.dyn_gnsf_ipiv_x;
idx_perm_x = model.dyn_gnsf_idx_perm_x;
ipiv_z = model.dyn_gnsf_ipiv_z;
idx_perm_z = model.dyn_gnsf_idx_perm_z;
ipiv_f = model.dyn_gnsf_ipiv_f;
idx_perm_f = model.dyn_gnsf_idx_perm_f;

% CasADi variables and expressions
% x
x = model.sym_x;
x1 = x(model.dyn_gnsf_idx_perm_x(1:model.dim_gnsf_nx1));
% check type
if class(x(1)) == 'casadi.SX'
    isSX = true;
else
    isSX = false;
end
% xdot
xdot = model.sym_xdot;
x1dot = xdot(model.dyn_gnsf_idx_perm_x(1:model.dim_gnsf_nx1));
% u
if isfield(model, 'sym_u')
    u = model.sym_u;
else
    if isSX
        u = SX.sym('u',0, 0);
    else
        u = MX.sym('u',0, 0);
    end
end
% z
if isfield(model, 'sym_z')
    z = model.sym_z;
	z1 = model.sym_z(model.dyn_gnsf_idx_perm_z(1:model.dim_gnsf_nz1));
else
    if isSX
        z = SX.sym('z',0, 0);
        z1 = SX.sym('z1',0, 0);
    else
        z = MX.sym('z',0, 0);
        z1 = MX.sym('z1',0, 0);
    end
end
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
% y
y = model.sym_gnsf_y;
% uhat
uhat = model.sym_gnsf_uhat;

% expressions
phi = model.dyn_gnsf_expr_phi;
f_lo = model.dyn_gnsf_expr_f_lo;

nontrivial_f_LO = model.dyn_gnsf_nontrivial_f_LO;
purely_linear = model.dyn_gnsf_purely_linear;

% name
model_name = model.name;
model_name = [model_name, '_dyn'];

% generate functions
jac_phi_y = jacobian(phi,y);
jac_phi_uhat = jacobian(phi,uhat);

if (strcmp(model.dyn_param_f, 'true'))
    phi_fun = Function([model_name,'_gnsf_phi_fun'], {y, uhat, p}, {phi});
    phi_fun_jac_y = Function([model_name,'_gnsf_phi_fun_jac_y'], {y, uhat, p}, {phi, jac_phi_y});
    phi_jac_y_uhat = Function([model_name,'_gnsf_phi_jac_y_uhat'], {y, uhat, p}, {jac_phi_y, jac_phi_uhat});

    f_lo_fun_jac_x1k1uz = Function([model_name,'_gnsf_f_lo_fun_jac_x1k1uz'], {x1, x1dot, z1, u, p}, ...
        {f_lo, [jacobian(f_lo,x1), jacobian(f_lo,x1dot), jacobian(f_lo,u), jacobian(f_lo,z1)]});
else
    phi_fun = Function([model_name,'_gnsf_phi_fun'], {y, uhat}, {phi});
    phi_fun_jac_y = Function([model_name,'_gnsf_phi_fun_jac_y'], {y, uhat}, {phi, jac_phi_y});
    phi_jac_y_uhat = Function([model_name,'_gnsf_phi_jac_y_uhat'], {y, uhat}, {jac_phi_y, jac_phi_uhat});

    f_lo_fun_jac_x1k1uz = Function([model_name,'_gnsf_f_lo_fun_jac_x1k1uz'], {x1, x1dot, z1, u}, ...
        {f_lo, [jacobian(f_lo,x1), jacobian(f_lo,x1dot), jacobian(f_lo,u), jacobian(f_lo,z1)]});
end

% get_matrices function
dummy = x(1);

get_matrices_fun = Function([model_name,'_gnsf_get_matrices_fun'], {dummy},...
     {A, B, C, E, L_x, L_xdot, L_z, L_u, A_LO, c, E_LO, B_LO, nontrivial_f_LO, purely_linear, ipiv_x, ipiv_z});


%% generate functions
f_lo_fun_jac_x1k1uz.generate([model_name,'_gnsf_f_lo_fun_jac_x1k1uz'], casadi_opts);
phi_fun.generate([model_name,'_gnsf_phi_fun'], casadi_opts);
phi_fun_jac_y.generate([model_name,'_gnsf_phi_fun_jac_y'], casadi_opts);
phi_jac_y_uhat.generate([model_name,'_gnsf_phi_jac_y_uhat'], casadi_opts);
get_matrices_fun.generate([model_name,'_gnsf_get_matrices_fun'], casadi_opts);

end
