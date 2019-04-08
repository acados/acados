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
omega = SX.sym('omega'); % angular velocity of rod [rad/s]
F = SX.sym('F');         % horizontal force acting on cart [N]

%% (unnamed) symbolic variables
sym_x = vertcat(p, theta, v, omega);
sym_xdot = SX.sym('xdot', nx, 1);
sym_u = F;

%% dynamics
expr_f_expl = vertcat(v, ...
                      omega, ...
                      (- l*m*sin(theta)*omega.^2 + F + g*m*cos(theta)*sin(theta))/(M + m - m*cos(theta).^2), ...
                      (- l*m*cos(theta)*sin(theta)*omega.^2 + F*cos(theta) + g*m*sin(theta) + M*g*sin(theta))/(l*(M + m - m*cos(theta).^2)));
expr_f_impl = expr_f_expl - sym_xdot;

%% populate structure
model.nx = nx;
model.nu = nu;
model.sym_x = sym_x;
model.sym_xdot = sym_xdot;
model.sym_u = sym_u;
model.expr_f_expl = expr_f_expl;
model.expr_f_impl = expr_f_impl;

