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

expr_f_expl = Ac*sym_x + Bc*sym_u + c_const;
expr_f_impl = expr_f_expl - sym_xdot;
expr_phi = A*sym_x + B*sym_u;

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
model.expr_phi = expr_phi;
model.expr_h = expr_h;
model.expr_h_e = expr_h_e;
model.expr_y = expr_y;
model.expr_y_e = expr_y_e;
model.expr_ext_cost = expr_ext_cost;
model.expr_ext_cost_e = expr_ext_cost_e;

