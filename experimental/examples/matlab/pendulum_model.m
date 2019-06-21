function [ode_model, nx, nu] = pendulum_model()

import casadi.*

M = 1;    % mass of the cart [kg]
m = 0.1;  % mass of the ball [kg]
g = 9.81; % gravity constant [m/s^2]
l = 0.8;  % length of the rod [m]

p = SX.sym('p');         % horizontal displacement [m]
theta = SX.sym('theta'); % angle with the vertical [rad]
v = SX.sym('v');         % horizontal velocity [m/s]
omega = SX.sym('omega'); % angular velocity [rad/s]
F = SX.sym('F');         % horizontal force [N]

ode_rhs = vertcat(v, ...
                  omega, ...
                  (- l*m*sin(theta)*omega.^2 + F + g*m*cos(theta)*sin(theta))/(M + m - m*cos(theta).^2), ...
                  (- l*m*cos(theta)*sin(theta)*omega.^2 + F*cos(theta) + g*m*sin(theta) + M*g*sin(theta))/(l*(M + m - m*cos(theta).^2)));

x = vertcat(p, theta, v, omega);
nx = length(x);
nu = 1;

ode_model = Function('pendulum', {x, F}, {ode_rhs}, {'x', 'u'}, {'rhs'});

end