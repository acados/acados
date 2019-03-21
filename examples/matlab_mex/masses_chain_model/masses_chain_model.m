function model = masses_chain_model(nfm)

import casadi.*

% number of masses
nm = nfm+1; %3;

% Environment
g = 9.81;     % [N/kg]
L = 0.033;
D = 1.0;
m = 0.03;
x0 = zeros(3,1);

% dims
nx = (nm-1)*2*3;
nu = 3;



%% symbolic variables
sym_u = SX.sym('u', nu, 1); % controls
sym_x = [];
str_x = [];
for ii=1:nm-1
	p = SX.sym(['p' num2str(ii)], 3);
	v = SX.sym(['v' num2str(ii)], 3);
	tmp_x = struct('p', p, 'v', v);
	str_x = [str_x; tmp_x];
	sym_x = [sym_x; casadi_struct2vec(tmp_x)];
end
sym_xdot = SX.sym('xdot',size(sym_x)); %state derivatives



%% dynamics
% compute forces
F = {};
for ii=1:nm-1
	if ii==1
		dist = str_x(1).p - x0;
	else
		dist = str_x(ii).p - str_x(ii-1).p;
	end
	tmp = D * (1 - L/sqrt(dist.'*dist));
	F = {F{:}, tmp*dist};
end
% setup ode
expr_f_expl = [];
for ii=1:nm-2
	f = 1/m * (F{ii+1} - F{ii}) - [0; 0; g];
	expr_f_expl = [expr_f_expl; casadi_vec(tmp_x, 'p', str_x(ii).v, 'v', f)];
end
expr_f_expl = [expr_f_expl; casadi_vec(tmp_x, 'p', str_x(end).v, 'v', sym_u)];

expr_f_impl = expr_f_expl - sym_xdot;

% function to compute rest positions
f_fun_jac = Function('f_fun_jac', {sym_x, sym_u}, {expr_f_expl, jacobian(expr_f_expl, [sym_x; sym_u])});

% compute end rest position => reference
x_ref = [linspace(0, 1, nm); zeros(5, nm)];
x_ref = x_ref(:,2:end);
x_ref = x_ref(:);
u_ref = zeros(nu, 1);

[fun_tmp, jac_tmp] = f_fun_jac(x_ref, u_ref);
while norm(full(fun_tmp)) > 1e-10
	fun = full(fun_tmp);
	jac = full(jac_tmp);
	delta = - jac \ fun;
	x_ref = x_ref + delta(1:nx);
	[fun_tmp, jac_tmp] = f_fun_jac(x_ref, u_ref);
	norm(full(fun_tmp));
end
%x_ref

% compute start rest position => x0
x0 = [zeros(1, nm); linspace(0, 1.5, nm); linspace(0, 0.5, nm); zeros(3, nm)];
x0 = x0(:,2:end);
x0 = x0(:);
u0 = zeros(nu, 1);

[fun_tmp, jac_tmp] = f_fun_jac(x0, u0);
while norm(full(fun_tmp)) > 1e-10
	fun = full(fun_tmp);
	jac = full(jac_tmp);
	delta = - jac \ fun;
	x0 = x0 + delta(1:nx);
	[fun_tmp, jac_tmp] = f_fun_jac(x0, u0);
	norm(full(fun_tmp));
end
%x0



%% cost
expr_y = [sym_x; sym_u];
expr_y_e = [sym_x];



%% constraints
expr_h = SX.zeros(0);
expr_h_e = SX.zeros(0);



%% populate structure
mode = struct;
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
model.x0 = x0;
model.x_ref = x_ref;

