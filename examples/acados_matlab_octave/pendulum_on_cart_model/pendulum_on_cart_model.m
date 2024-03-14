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

function model = pendulum_on_cart_model()

import casadi.*

%% system dimensions
nx = 4;
nu = 1;

%% system parameters
M = 1;    % mass of the cart [kg]
m = 0.1;  % mass of the ball [kg]
l = 0.8;  % length of the rod [m]
g = 9.81; % gravity constant [m/s^2]

%% named symbolic variables
p = SX.sym('p');         % horizontal displacement of cart [m]
theta = SX.sym('theta'); % angle of rod with the vertical [rad]
v = SX.sym('v');         % horizontal velocity of cart [m/s]
dtheta = SX.sym('dtheta'); % angular velocity of rod [rad/s]
F = SX.sym('F');         % horizontal force acting on cart [N]

%% (unnamed) symbolic variables
sym_x = vertcat(p, theta, v, dtheta);
sym_xdot = SX.sym('xdot', nx, 1);
sym_u = F;

%% dynamics
%expr_f_expl = vertcat(v, ...
%                      dtheta, ...
%                      (- l*m*sin(theta)*dtheta.^2 + F + g*m*cos(theta)*sin(theta))/(M + m - m*cos(theta).^2), ...
%                      (- l*m*cos(theta)*sin(theta)*dtheta.^2 + F*cos(theta) + g*m*sin(theta) + M*g*sin(theta))/(l*(M + m - m*cos(theta).^2)));
sin_theta = sin(theta);
cos_theta = cos(theta);
denominator = M + m - m*cos_theta.^2;
dyn_expr_f_expl = vertcat(v, ...
						 dtheta, ...
                         (- l*m*sin_theta*dtheta.^2 + F + g*m*cos_theta*sin_theta)/denominator, ...
                         (- l*m*cos_theta*sin_theta*dtheta.^2 + F*cos_theta + g*m*sin_theta + M*g*sin_theta)/(l*denominator));
dyn_expr_f_impl = dyn_expr_f_expl - sym_xdot;

%% constraints
constr_expr_h_0 = sym_u;
constr_expr_h = sym_u;

%% cost
% generic cost
W_x = diag([1e3, 1e3, 1e-2, 1e-2]);
W_u = 1e-2;
cost_expr_ext_cost_e = 0.5 * sym_x'* W_x * sym_x;
cost_expr_ext_cost = cost_expr_ext_cost_e + 0.5 * sym_u' * W_u * sym_u;
cost_expr_ext_cost_0 = 0.5 * sym_u' * W_u * sym_u;

% nonlinear least sqares
cost_expr_y_0 = sym_u;
cost_W_0 = W_u;
cost_expr_y = vertcat(sym_x, sym_u);
cost_W = blkdiag(W_x, W_u);
cost_expr_y_e = sym_x;
cost_W_e = W_x;

% linear least squares
ny_0 = nu; % number of outputs in initial cost term
cost_Vx_0 = zeros(ny_0,nx);
cost_Vu_0 = eye(nu);
cost_y_ref_0 = zeros(ny_0, 1);

ny = nx+nu; % number of outputs in lagrange term
cost_Vx = [eye(nx); zeros(nu,nx)]; % state-to-output matrix in lagrange term
cost_Vu = [zeros(nx, nu); eye(nu)]; % input-to-output matrix in lagrange term
cost_y_ref = zeros(ny, 1); % output reference in lagrange term
    
ny_e = nx; % number of outputs in terminal cost term
cost_Vx_e = eye(ny_e, nx);
cost_y_ref_e = zeros(ny_e, 1);


%% populate structure
model.nx = nx;
model.nu = nu;
model.sym_x = sym_x;
model.sym_xdot = sym_xdot;
model.sym_u = sym_u;

model.dyn_expr_f_expl = dyn_expr_f_expl;
model.dyn_expr_f_impl = dyn_expr_f_impl;

model.constr_expr_h_0 = constr_expr_h_0;
model.constr_expr_h = constr_expr_h;

model.cost_expr_ext_cost_0 = cost_expr_ext_cost_0;
model.cost_expr_ext_cost = cost_expr_ext_cost;
model.cost_expr_ext_cost_e = cost_expr_ext_cost_e;

model.cost_expr_y_0 = cost_expr_y_0;
model.cost_W_0 = cost_W_0;
model.cost_expr_y = cost_expr_y;
model.cost_W = cost_W;
model.cost_expr_y_e = cost_expr_y_e;
model.cost_W_e = cost_W_e;

model.cost_Vx_0 = cost_Vx_0;
model.cost_Vu_0 = cost_Vu_0;
model.cost_y_ref_0 = cost_y_ref_0;
model.cost_Vx = cost_Vx;
model.cost_Vu = cost_Vu;
model.cost_y_ref = cost_y_ref;
model.cost_Vx_e = cost_Vx_e;
model.cost_y_ref_e = cost_y_ref_e;

end
