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



function model = pendulum_on_cart_model_dae()

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
omega = SX.sym('omega'); % angular velocity of rod [rad/s]
F = SX.sym('F');         % horizontal force acting on cart [N]

sym_z = SX.sym('z');

%% (unnamed) symbolic variables
sym_x = vertcat(p, theta, v, omega);
sym_xdot = SX.sym('xdot', nx, 1);
sym_u = F;

%% dynamics
expr_f_impl = vertcat(v, ...
                      omega, ...
                      (- l*m*sin(theta)*omega.^2 + F + g*m*cos(theta)*sin(theta))/(M + m - m*cos(theta).^2), ...
                      (- l*m*sym_z + F*cos(theta) + g*m*sin(theta) + M*g*sin(theta))/(l*(M + m - m*cos(theta).^2)), ...
                      cos(theta)*sin(theta)*omega.^2) ...
                  - [sym_xdot; sym_z];

%% constraints
expr_h = sym_u;

%% cost
W_x = diag([1e3, 1e3, 1e-2, 1e-2]);
W_u = 1e-2;
expr_ext_cost_e = sym_x'* W_x * sym_x;
expr_ext_cost = expr_ext_cost_e + sym_u' * W_u * sym_u;

%% populate structure
model.nx = nx;
model.nu = nu;
model.sym_x = sym_x;
model.sym_xdot = sym_xdot;
model.sym_z = sym_z;
model.sym_u = sym_u;
model.dyn_expr_f_impl = expr_f_impl;
model.constr_expr_h = expr_h;
model.cost_expr_ext_cost = expr_ext_cost;
model.cost_expr_ext_cost_e = expr_ext_cost_e;

end