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
function generate_c_code_gnsf(gnsf)
%% import casadi
import casadi.*

if CasadiMeta.version()=='3.4.0'
	% casadi 3.4
	casadi_opts = struct('mex', false, 'casadi_int', 'int', 'casadi_real', 'double');
else
	% old casadi versions
	error('Please download and install Casadi 3.4.0 to ensure compatibility with acados')
end

% nx  = gnsf.nx;
% nu  = gnsf.nu;
% nz  = gnsf.nz;
% 
% nx1 = gnsf.nx1;
% nx2 = gnsf.nx2;
% n_out = gnsf.n_out;
% ny = gnsf.ny;
% nuhat = gnsf.nuhat;

A  = gnsf.A;
B  = gnsf.B;
C  = gnsf.C;
E  = gnsf.E;
c  = gnsf.c;

L_x    = gnsf.L_x;
L_z    = gnsf.L_z;
L_xdot = gnsf.L_xdot;
L_u    = gnsf.L_u;


A_LO = gnsf.A_LO;
f_lo = gnsf.f_lo_expr;

% functions
phi = gnsf.phi_expr;
y = gnsf.y;
uhat = gnsf.uhat;

jac_phi_y = jacobian(phi,y);
jac_phi_uhat = jacobian(phi,uhat);

model_name_prefix = gnsf.name;

phi_fun = Function([model_name_prefix,'phi_fun'], {y,uhat}, {phi});
phi_fun_jac_y = Function([model_name_prefix,'phi_fun_jac_y'], {y,uhat}, {phi, jac_phi_y});
phi_jac_y_uhat = Function([model_name_prefix,'phi_jac_y_uhat'], {y,uhat}, {jac_phi_y, jac_phi_uhat});

phi_jac_y = Function([model_name_prefix,'phi_jac_y_uhat'], {y,uhat}, {jac_phi_y});

% Linear output
x1 = gnsf.x(1:gnsf.nx1);
x1dot = gnsf.xdot(1:gnsf.nx1);
u = gnsf.u;
z = gnsf.z;

% x2 = gnsf.x(1+gnsf.nx1 : gnsf.nx);

jac_f_x1 = jacobian(f_lo,x1);
jac_f_u  = jacobian(f_lo,u);
jac_f_z  = jacobian(f_lo,z);
jac_f_k1 = jacobian(f_lo,x1dot);

f_lo_fun_jac_x1k1uz = Function([model_name_prefix,'f_lo_fun_jac_x1k1uz'], {x1, x1dot, z, u}, ...
    {f_lo, [jac_f_x1, jac_f_k1, jac_f_z, jac_f_u]});

% f_lo_fun = Function([model_name_prefix,'f_lo_fun_jac_x1k1uz'], {x1, x1dot, z, u}, {f_lo});


%% generate functions
dummy = SX.sym('dummy');

model_matrices = SX.zeros(size([A(:); B(:); C(:); E(:); L_x(:); L_xdot(:); L_z(:); L_u(:); A_LO(:); c(:)])) + ...
    [A(:); B(:); C(:); E(:); L_x(:); L_xdot(:); L_z(:); L_u(:); A_LO(:); c(:)];
get_matrices_fun = Function([model_name_prefix,'get_matrices_fun'], {dummy}, {model_matrices(:)});
get_matrices_fun.generate([model_name_prefix,'get_matrices_fun'], casadi_opts);

% generate phi, f_LO
f_lo_fun_jac_x1k1uz.generate([model_name_prefix,'f_lo_fun_jac_x1k1uz'], casadi_opts);
phi_fun.generate([model_name_prefix,'phi_fun'], casadi_opts);
phi_fun_jac_y.generate([model_name_prefix,'phi_fun_jac_y'], casadi_opts);
phi_jac_y_uhat.generate([model_name_prefix,'phi_jac_y_uhat'], casadi_opts);

end
