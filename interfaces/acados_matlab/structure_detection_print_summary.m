function structure_detection_print_summary(gnsf, model)
%
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

%% Description
% this function prints the most important info after determining a GNSF
% reformulation of the implicit model "initial_model" into "gnsf", which is
% equivalent to the "reordered_model".

% % GNSF
% get dimensions
nx  = gnsf.nx;
nu  = gnsf.nu;
nz  = gnsf.nz;

nx1 = gnsf.nx1;
nx2 = gnsf.nx2;

nz1 = gnsf.nz1;
nz2 = gnsf.nz2;

np = gnsf.np;
n_out = gnsf.n_out;
ny = gnsf.ny;
nuhat = gnsf.nuhat;

%
f_impl_expr = model.dyn_expr_f;
n_nodes_initial = model.dyn_expr_f.n_nodes();
%x_old = model.x;
%f_impl_old = model.f_impl_expr;

x = gnsf.x;
z = gnsf.z;

phi_current = gnsf.phi_expr;

%% PRINT SUMMARY -- STRUCHTRE DETECTION
disp(' ');
disp('*********************************************************************************************');
disp(' ');
disp('******************        SUCCESS: GNSF STRUCTURE DETECTION COMPLETE !!!      ***************');
disp(' ');
disp('*********************************************************************************************');
disp(' ');
disp(['========================= STRUCTURE DETECTION SUMMARY ====================================']);
disp(' ');
disp('-------- Nonlinear Static Feedback type system --------');
disp(' ');
disp(' successfully transcribed dynamic system model into GNSF structure ');
disp(' ');
disp(['reduced dimension of nonlinearity phi from        ', sprintf('%6s', num2str(nx+nz)),      ' to ', sprintf('%6s', num2str(gnsf.n_out))]);
disp(' ');
disp(['reduced input dimension of nonlinearity phi from  ', sprintf('%6s', num2str(2*nx+nz+nu)), ' to ', sprintf('%6s', num2str(gnsf.ny + gnsf.nuhat))]);
disp(' ');
disp(['reduced number of nodes in CasADi expression of']);
disp(['nonlinearity phi from                             '...
    , sprintf('%6s', num2str(n_nodes_initial)), ' to ', sprintf('%6s', num2str(phi_current.n_nodes()))]);
disp(' ');
disp('----------- Linear Output System (LOS) ---------------');
if gnsf.nx2 + gnsf.nz2 >0
    disp(' ');
    disp(['introduced Linear Output System of size           ', sprintf('%6s', num2str(gnsf.nx2 + gnsf.nz2)),'']);
    disp(' ');
end
if gnsf.nx2 >0
    disp('consisting of the states:');
    disp(' ');
    disp(x(gnsf.nx1+1:gnsf.nx));
    disp(' ');
end
if gnsf.nz2 >0
    disp('and algebraic variables:');
    disp(' ');
    disp(z(gnsf.nz1+1:gnsf.nz));
    disp(' ');
end
if gnsf.purely_linear == 1
    disp(' ');
    disp('Model is fully linear!');
	disp(' ');
end

%compare_x = (x_old == x);
%if ~compare_x.is_constant()
%disp(' ');
%disp('--------------------------------------------------------------------------------------------------');
%disp('NOTE: permuted state vector x, such that the implicit model, can take it in the same order as GNSF');
%disp(' ');
%disp(' OLD / initial state vector read as: ');
%disp(x_old);
%disp(' ');
%disp(' whereas NEW / permuted state vector reads as: ');
%disp(x);
%end
% TODO permutation information

%compare_f = (f_impl_old == f_impl_expr);
%if ~compare_f.is_constant()
%disp(' ');
%disp('--------------------------------------------------------------------------------------------------');
%disp('NOTE: permuted implicit function, such that the implicit & GNSF model, have the same order');
%disp(' ');
%disp(' OLD / initial state f_impl read as: ');
%disp(' ');
%print_casadi_expression(f_impl_old);
%disp(' ');
%disp(' whereas NEW / permuted f_impl reads as: ');
%disp(' ');
%print_casadi_expression(f_impl_expr);
%disp(' ');
%end
if ~isequal(gnsf.idx_perm_x, [1:nx])
	disp(' ');
	disp('--------------------------------------------------------------------------------------------------');
	disp('NOTE: permuted differential state vector x, such that x_gnsf = x(idx_perm_x) with idx_perm_x =');
	gnsf.idx_perm_x
	disp(' ');
end

if nz~= 0 && ~isequal(gnsf.idx_perm_z, [1:nz])
	disp(' ');
	disp('--------------------------------------------------------------------------------------------------');
	disp('NOTE: permuted algebraic state vector z, such that z_gnsf = z(idx_perm_z) with idx_perm_z =');
	gnsf.idx_perm_z
	disp(' ');
end

if ~isequal(gnsf.idx_perm_f, [1:nx+nz])
	disp(' ');
	disp('--------------------------------------------------------------------------------------------------');
	disp('NOTE: permuted rhs expression vector f, such that f_gnsf = f(idx_perm_f) with idx_perm_f =');
	gnsf.idx_perm_f
	disp(' ');
end

%% print GNSF dimensions
format short
disp('--------------------------------------------------------------------------------------------------------');
disp(' ');
disp('The dimensions of the GNSF reformulated model read as:');
disp(' ');
%T_dim = table(nx, nu, nz, np, nx1, nz1, n_out, ny, nuhat);
%disp( T_dim )
nx
nu
nz
np
nz1
nz1
n_out
ny
nuhat
format short e

end
