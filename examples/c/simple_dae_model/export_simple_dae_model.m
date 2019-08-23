%
% Copyright 2019 Gianluca Frison, Dimitris Kouzoupis, Robin Verschueren,
% Andrea Zanelli, Niels van Duijkeren, Jonathan Frey, Tommaso Sartor,
% Branimir Novoselnik, Rien Quirynen, Rezart Qelibari, Dang Doan,
% Jonas Koenemann, Yutao Chen, Tobias Schöls, Jonas Schlagenhauf, Moritz Diehl
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

function [ model ] = export_simple_dae_model()
    %% this function generates an implicit ODE / index-1 DAE model,
    % which consists of a CasADi expression f_impl_expr
    % that depends on the symbolic CasADi variables x, xdot, u, z,
    % and a model name, which will be used as a prefix for generated C
    % functions later on;
    
    %% CasADi
    import casadi.*
    casadi_version = CasadiMeta.version();
    if strcmp(casadi_version(1:3),'3.4') % require casadi 3.4.x
        casadi_opts = struct('mex', false, 'casadi_int', 'int', 'casadi_real', 'double');
    else % old casadi versions
        error('Please download and install CasADi version 3.4.x to ensure compatibility with acados')
    end
    model_name_prefix = 'simple_dae';
       
    %% parameters
    alpha = 1.0;
    
    %% Set up States & Controls
    x1    = SX.sym('x1');     % Differential States
    x2    = SX.sym('x2');
    x = vertcat(x1, x2);
    
    z1      = SX.sym('z1');     % Algebraic states
    z2      = SX.sym('z2');
    z = vertcat(z1, z2);
    
    u1      = SX.sym('u1');     % Algebraic states
    u2      = SX.sym('u2');
    u       = vertcat(u1, u2);
    
    %% xdot
    x1_dot    = SX.sym('x1_dot');     % Differential States
    x2_dot    = SX.sym('x2_dot');
    
    xdot = [x1_dot; x2_dot];
    
    %% Dynamics: implicit DAE formulation (index-1)
    % x = vertcat(xpos, ypos, alpha, vx, vy, valpha);
    % z = vertcat(ax, ay, aalpha, Fx, Fy);
    f_impl = vertcat(x1_dot+x1-0.1*z2-u1, ...
                     x2_dot+x2-0.1*z1-u2,  ...
                     z1-x1, ...
                     z2-x2);
    
    %% initial value
    %     x0 = [0.1; -0.1];
    %     z0 = [0.0, 0.0];
    %     u0 = 0;

    model.f_impl_expr = f_impl;
    model.x = x;
    model.xdot = xdot;
    model.u = u;
    model.z = z;
    model.name = model_name_prefix;
    
end
