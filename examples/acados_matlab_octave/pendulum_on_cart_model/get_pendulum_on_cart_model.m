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
