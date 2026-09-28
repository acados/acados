%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.





% NOTE: `acados` currently supports both an old MATLAB/Octave interface (< v0.4.0)
% as well as a new interface (>= v0.4.0).

% THIS EXAMPLE still uses the OLD interface. If you are new to `acados` please start
% with the examples that have been ported to the new interface already.
% see https://github.com/acados/acados/issues/1196#issuecomment-2311822122)


function model = linear_mass_spring_model()

import casadi.*

%% dims
num_mass = 4;

nx = 2*num_mass;
nu = num_mass-1;

%% symbolic variables
if 1
	sym_x = SX.sym('x', nx, 1); % states
	sym_u = SX.sym('u', nu, 1); % controls
	sym_xdot = SX.sym('xdot',size(sym_x)); %state derivatives
else
	sym_x = MX.sym('x', nx, 1); % states
	sym_u = MX.sym('u', nu, 1); % controls
	sym_xdot = MX.sym('xdot',size(sym_x)); %state derivatives
end

%% dynamics
% continuous time
Ac = zeros(nx, nx);
for ii=1:num_mass
	Ac(ii,num_mass+ii) = 1.0;
	Ac(num_mass+ii,ii) = -2.0;
end
for ii=1:num_mass-1
	Ac(num_mass+ii,ii+1) = 1.0;
	Ac(num_mass+ii+1,ii) = 1.0;
end

Bc = zeros(nx, nu);
for ii=1:nu
	Bc(num_mass+ii, ii) = 1.0;
end

c_const = zeros(nx, 1);
% just to test gnsf with nontrivial c_LO
% for ii=1:nx
%     c_const(ii) = (-1)^ii * 1e-8;
% end

% discrete time
Ts = 0.5; % sampling time
M = expm([Ts*Ac, Ts*Bc; zeros(nu, 2*nx/2+nu)]);
A = M(1:nx,1:nx);
B = M(1:nx,nx+1:end);

dyn_expr_f_expl = Ac*sym_x + Bc*sym_u + c_const;
dyn_expr_f_impl = dyn_expr_f_expl - sym_xdot;
dyn_expr_phi = A*sym_x + B*sym_u;

%% constraints
constr_expr_h_0 = sym_u;
constr_expr_h = [sym_u; sym_x];
constr_expr_h_e = sym_x;

%% nonlinear least squares
cost_expr_y_0 = sym_u;
cost_expr_y = [sym_u; sym_x];
cost_expr_y_e = sym_x;

%% external cost
yr_u = zeros(nu, 1);
yr_x = zeros(nx, 1);
dWu = 2*ones(nu, 1);
dWx = ones(nx, 1);

ymyr_0 = sym_u - yr_u;
ymyr = [sym_u; sym_x] - [yr_u; yr_x];
ymyr_e = sym_x - yr_x;

cost_expr_ext_cost_0 = 0.5 * ymyr_0' * (dWu .* ymyr_0);
cost_expr_ext_cost = 0.5 * ymyr' * ([dWu; dWx] .* ymyr);
cost_expr_ext_cost_e = 0.5 * ymyr_e' * (dWx .* ymyr_e);

%% populate structure
model.nx = nx;
model.nu = nu;
model.sym_x = sym_x;
model.sym_xdot = sym_xdot;
model.sym_u = sym_u;
model.dyn_expr_f_expl = dyn_expr_f_expl;
model.dyn_expr_f_impl = dyn_expr_f_impl;
model.dyn_expr_phi = dyn_expr_phi;
model.constr_expr_h_0 = constr_expr_h_0;
model.constr_expr_h = constr_expr_h;
model.constr_expr_h_e = constr_expr_h_e;
model.cost_expr_y_0 = cost_expr_y_0;
model.cost_expr_y = cost_expr_y;
model.cost_expr_y_e = cost_expr_y_e;
model.cost_expr_ext_cost_0 = cost_expr_ext_cost_0;
model.cost_expr_ext_cost = cost_expr_ext_cost;
model.cost_expr_ext_cost_e = cost_expr_ext_cost_e;
end