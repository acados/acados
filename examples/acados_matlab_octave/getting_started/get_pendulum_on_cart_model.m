%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.

function model = get_pendulum_on_cart_model(varargin)
    % GET_PENDULUM_ON_CART_MODEL returns AcadosModel describing the pendulum on cart dynamics.
    % MODEL = GET_PENDULUM_ON_CART_MODEL() returns AcadosModel with implicit and explicit continuous dynamics defined.
    % MODEL = GET_PENDULUM_ON_CART_MODEL(DELTA_T) returns AcadosModel with implicit and explicit continuous dynamics and discrete dynamics defined.
    % MODEL = GET_PENDULUM_ON_CART_MODEL(DELTA_T, WITH_PARAM) returns AcadosModel with implicit and explicit continuous dynamics and discrete dynamics defined. Also the mass M is modeled as a parameter.

    import casadi.*

    %% system dimensions
    nx = 4;
    nu = 1;

    %% system parameters
    if nargin > 1 && varargin{2} % parametric model
        M = SX.sym('M');  % mass of the cart [kg]
        param = M;
    else
        M = 1;    % mass of the cart [kg]
        param = [];
    end

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

    % discrete dynamics
    if nargin > 0
        delta_t = varargin{1};
        disc_dyn_expr = x + delta_t * f_expl_expr; % explicit Euler
    else
        disc_dyn_expr = [];
    end

    % populate
    model = AcadosModel();
    model.x = x;
    model.xdot = xdot;
    model.u = u;
    model.p = param;

    model.f_expl_expr = f_expl_expr;
    model.f_impl_expr = f_impl_expr;
    model.disc_dyn_expr = disc_dyn_expr;
    model.name = 'pendulum';
end
