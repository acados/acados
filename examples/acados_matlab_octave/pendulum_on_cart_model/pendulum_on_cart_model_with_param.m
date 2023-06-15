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

function model = pendulum_on_cart_model_with_param()

import casadi.*

%% system dimensions
nx = 4;
nu = 1;

%% system parameters
M = SX.sym('M');  % mass of the cart [kg]
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
expr_f_expl = vertcat(v, ...
                      dtheta, ...
                      (- l*m*sin_theta*dtheta.^2 + F + g*m*cos_theta*sin_theta)/denominator, ...
                      (- l*m*cos_theta*sin_theta*dtheta.^2 + F*cos_theta + g*m*sin_theta + M*g*sin_theta)/(l*denominator));
expr_f_impl = expr_f_expl - sym_xdot;

%% constraints
expr_h = sym_u;

%% cost
W_x = diag([1e3, 1e3, 1e-2, 1e-2]);
W_u = 1e-2;
expr_ext_cost_e = sym_x'* W_x * sym_x;
expr_ext_cost = expr_ext_cost_e + sym_u' * W_u * sym_u;
% nonlinear least sqares
cost_expr_y = vertcat(sym_x, sym_u);
W = blkdiag(W_x, W_u);
model.cost_expr_y_e = sym_x;
model.W_e = W_x;


sym_p = M;

%% populate structure
model.nx = nx;
model.nu = nu;
model.sym_x = sym_x;
model.sym_xdot = sym_xdot;
model.sym_u = sym_u;
model.sym_p = sym_p;
model.expr_f_expl = expr_f_expl;
model.expr_f_impl = expr_f_impl;
model.expr_h = expr_h;
model.expr_ext_cost = expr_ext_cost;
model.expr_ext_cost_e = expr_ext_cost_e;

model.cost_expr_y = cost_expr_y;
model.W = W;

end
