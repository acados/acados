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

function [ model ] = pendulum_dae_model()
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

    model.expr_f_impl = f_impl;
    model.sym_x = x;
    model.sym_xdot = xdot;
    model.sym_u = u;
    model.sym_z = z;
    model.name = model_name_prefix;
    model.expr_h = expr_h;
    
end