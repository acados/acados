%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.

%

function [ model ] = get_pendulum_dae_model()
    %% this function generates an implicit ODE / index-1 DAE model,
    % which consists of a CasADi expression f_impl_expr
    % that depends on the symbolic CasADi variables x, xdot, u, z,
    % and a model name, which will be used as a prefix for generated C
    % functions later on;
    
    %% CasADi
    import casadi.*

    model_name_prefix = 'pendulum_dae';
    
    %% Parameters (taken from Rien Quirynen's Master Thesis)
    % NOTE: removed torque from parameters and made it the control.
    m = 2;  % mass
    g = 9.81; % grav const
    I = 0.1; % moment of inertia

    %% Set up States & Controls
    xpos    = SX.sym('xpos');     % Differential States
    ypos    = SX.sym('ypos');
    alpha   = SX.sym('alpha');     
    vx      = SX.sym('vx');
    vy      = SX.sym('vy');
    valpha  = SX.sym('valpha');
    x = vertcat(xpos, ypos, alpha, vx, vy, valpha);
    
    ax      = SX.sym('ax');     % Algebraic states
    ay      = SX.sym('ay');
    aalpha  = SX.sym('aalpha');
    Fx      = SX.sym('Fx');
    Fy      = SX.sym('Fy');
    z = vertcat(ax, ay, aalpha, Fx, Fy);
    
    u       = SX.sym('u');  % Controls % applied torque
    
    %% xdot
    xpos_dot    = SX.sym('xpos_dot');     % Differential States
    ypos_dot    = SX.sym('ypos_dot');
    alpha_dot   = SX.sym('alpha_dot');     
    vx_dot      = SX.sym('vx_dot');
    vy_dot      = SX.sym('vy_dot');
    valpha_dot  = SX.sym('valpha_dot');
    
    xdot = [xpos_dot; ypos_dot; alpha_dot; vx_dot; vy_dot; valpha_dot];
    
    %% Dynamics: implicit DAE formulation (index-1)
    % x = vertcat(xpos, ypos, alpha, vx, vy, valpha);
    % z = vertcat(ax, ay, aalpha, Fx, Fy);
    f_impl = vertcat(xpos_dot - vx, ...
                     ypos_dot - vy, ...
                     alpha_dot - valpha, ...
                     vx_dot - ax, ...
                     vy_dot - ay, ...
                     valpha_dot - aalpha, ...
                     m * ax - Fx, ...
                     m * ay + m * g - Fy, ...
                     I * aalpha - u - Fx * ypos + Fy * xpos, ...
                     ax + vy * valpha + ypos * aalpha, ...
                     ay - vx * valpha - xpos * aalpha);
    %% constraint
    expr_h = ax^2 + ay^2;
                 
    %% initial value
%     x0 = [1; -5; 1; 0.1; -0.5; 0.1];
%     z0 = [-1.5; -0.3; -0.3; -3; 19];
%     u0 = 1;
    model = AcadosModel();
    model.f_impl_expr = f_impl;
    model.x = x;
    model.xdot = xdot;
    model.u = u;
    model.z = z;
    model.con_h_expr = expr_h;
    model.name = model_name_prefix;
end