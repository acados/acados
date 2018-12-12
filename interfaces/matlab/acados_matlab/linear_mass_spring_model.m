function model = linear_mass_spring_model()

import casadi.*

num_mass = 4;

nx = 2*num_mass;
nu = num_mass-1;

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

sym_x = MX.sym('x', nx, 1); % states
sym_u = MX.sym('u', nu, 1); % controls
sym_xdot = MX.sym('xdot',size(sym_x)); %state derivatives

expr_expl = Ac*sym_x + Bc*sym_u;
expr_impl = expr_expl - sym_xdot;

model.nx = nx;
model.nu = nu;
model.sym_x = sym_x;
model.sym_xdot = sym_xdot;
model.sym_u = sym_u;
model.expr_expl = expr_expl;
model.expr_impl = expr_impl;

