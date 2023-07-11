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

% discrete dynamics: casadi RK integrator

% Fixed step Runge-Kutta 4 integrator

M   = 2; % RK4 steps per interval
FUN = Function('f', {sym_x, sym_u}, {expr_f_expl});
DT  = SX.sym('DT', 1);
H   = DT/M;
X0  = sym_x; %SX.sym('X0', nx);
U   = sym_u; %SX.sym('U', nu);
X   = X0;
for j=1:M
   k1 = FUN(X, U);
   k2 = FUN(X + H/2 * k1, U);
   k3 = FUN(X + H/2 * k2, U);
   k4 = FUN(X + H * k3, U);
   X  = X+H/6*(k1 +2*k2 +2*k3 +k4);
end

np = 1;
sym_p = DT;
expr_phi = X;



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
model.np = np;
model.sym_x = sym_x;
model.sym_xdot = sym_xdot;
model.sym_u = sym_u;
model.sym_p = sym_p;
model.expr_f_expl = expr_f_expl;
model.expr_f_impl = expr_f_impl;
model.expr_phi = expr_phi;
model.expr_h = expr_h;
model.expr_h_e = expr_h_e;
model.expr_y = expr_y;
model.expr_y_e = expr_y_e;
model.x0 = x0;
model.x_ref = x_ref;

