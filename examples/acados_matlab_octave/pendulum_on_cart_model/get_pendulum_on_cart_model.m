%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.

function model = get_pendulum_on_cart_model()

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
    x = vertcat(p, theta, v, dtheta);
    xdot = SX.sym('xdot', nx, 1);
    u = F;

    sin_theta = sin(theta);
    cos_theta = cos(theta);
    denominator = M + m - m*cos_theta.^2;
    f_expl_expr = vertcat(v, ...
                             dtheta, ...
                             (- l*m*sin_theta*dtheta.^2 + F + g*m*cos_theta*sin_theta)/denominator, ...
                             (- l*m*cos_theta*sin_theta*dtheta.^2 + F*cos_theta + g*m*sin_theta + M*g*sin_theta)/(l*denominator));
    f_impl_expr = f_expl_expr - xdot;

    % populate
    model = AcadosModel();
    model.x = x;
    model.xdot = xdot;
    model.u = u;

    model.f_expl_expr = f_expl_expr;
    model.f_impl_expr = f_impl_expr;
    model.name = 'pendulum_on_cart';
end
