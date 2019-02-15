function model = linear_mass_spring_model()

import casadi.*

%% dims
num_mass = 4;

nx = 2*num_mass;
nu = num_mass-1;

%% symbolic variables
sym_x = MX.sym('x', nx, 1); % states
sym_u = MX.sym('u', nu, 1); % controls
sym_xdot = MX.sym('xdot',size(sym_x)); %state derivatives

%% dynamics
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

expr_f_expl = Ac*sym_x + Bc*sym_u;
expr_f_impl = expr_f_expl - sym_xdot;

%% constraints
expr_h = [sym_u; sym_x];
expr_h_e = [sym_x];

%% nonlnear least squares
expr_y = [sym_u; sym_x];
expr_y_e = [sym_x];

%% external cost
yr_u = zeros(nu, 1);
yr_x = zeros(nx, 1);
dWu = 2*ones(nu, 1);
dWx = ones(nx, 1);

ymyr = [sym_u; sym_x] - [yr_u; yr_x];
ymyr_e = sym_x - yr_x;

expr_ext_cost = 0.5 * ymyr' * ([dWu; dWx] .* ymyr);
expr_ext_cost_e = 0.5 * ymyr_e' * (dWx .* ymyr_e);

%% populate structure
model.nx = nx;
model.nu = nu;
model.sym_x = sym_x;
model.sym_xdot = sym_xdot;
model.sym_u = sym_u;
model.expr_f_expl = expr_f_expl;
model.expr_f_impl = expr_f_impl;
model.expr_h = expr_h;
model.expr_h_e = expr_h_e;
model.expr_y = expr_y;
model.expr_y_e = expr_y_e;
model.expr_ext_cost = expr_ext_cost;
model.expr_ext_cost_e = expr_ext_cost_e;

