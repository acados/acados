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

function [ model ] = export_pendulum_ode_model()
    %% this function generates an implicit ODE / index-1 DAE model,
    % which consists of a CasADi expression f_impl_expr, f_expl_expr
    % that depends on the symbolic CasADi variables x, xdot, u, z,
    % and a model name, which will be used as a prefix for generated C
    % functions later on;
%% this model is based on the explicit pendulum model
% but formulated implicitly to test implicit integrators with it.
    %% CasADi
    import casadi.*
    casadi_version = CasadiMeta.version();
    if ( strcmp(casadi_version(1:3),'3.4') || strcmp(casadi_version(1:3),'3.5')) % require casadi 3.4.x
        casadi_opts = struct('mex', false, 'casadi_int', 'int', 'casadi_real', 'double');
    else % old casadi versions
        error('Please provide CasADi version 3.4 or 3.5 to ensure compatibility with acados')
    end
    model_name = 'pendulum_ode';
    

    %% Constants
    M = 1;
    m = 0.1;
    g = 9.81;
    l = 0.8;

    %% Set up States & Controls
    x1      = SX.sym('x1');
    theta   = SX.sym('theta');
    v1      = SX.sym('v1');
    dtheta  = SX.sym('dtheta');
    
    x = [x1; v1; theta; dtheta];

    % Controls
    F = SX.sym('F');
    u = F;
    
    %% xdot
    x1_dot      = SX.sym('x1_dot');
    theta_dot   = SX.sym('theta_dot');
    v1_dot      = SX.sym('v1_dot');
    dtheta_dot  = SX.sym('dtheta_dot');

    xdot = [ x1_dot; theta_dot; v1_dot; dtheta_dot ];
    
    %% algebraic variables
    z = [];
    
    %% Dynamics     
    denominator = M + m - m*cos(theta)*cos(theta);
    f_expl = [  v1; ...
                (-m*l*sin(theta)*dtheta*dtheta + m*g*cos(theta)*sin(theta)+F)/denominator; ...
                dtheta; ...
                (-m*l*cos(theta)*sin(theta)*dtheta*dtheta + F*cos(theta)+(M+m)*g*sin(theta))/(l*denominator)];
    
    f_impl = xdot - f_expl;
    
    model.f_impl_expr = f_impl;
    model.f_expl_expr = f_expl;
    model.x = x;
    model.xdot = xdot;
    model.u = u;
    model.z = z;
    model.name = model_name;

end

