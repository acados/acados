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

function [ model ] = export_crane_dae_model()
    %% this function generates an implicit ODE / index-1 DAE model,
    % which consists of a CasADi expression f_impl_expr
    % that depends on the symbolic CasADi variables x, xdot, u, z,
    % and a model name, which will be used as a prefix for generated C
    % functions later on;
    
    %% CasADi
    import casadi.*
    casadi_version = CasadiMeta.version();
    if ( strcmp(casadi_version(1:3),'3.4') || strcmp(casadi_version(1:3),'3.5')) % require casadi 3.4.x
        casadi_opts = struct('mex', false, 'casadi_int', 'int', 'casadi_real', 'double');
    else % old casadi versions
        error('Please provide CasADi version 3.4 or 3.5 to ensure compatibility with acados')
    end
    model_name_prefix = 'crane_dae';
    
    %% Parameters (taken from Riens ACADO model)
    tau1 = 0.012790605943772;   a1   = 0.047418203070092;
    tau2 = 0.024695192379264;   a2   = 0.034087337273386;
    g = 9.81;

    %% Set up States & Controls
    xC = SX.sym('xC');     %States
    vC = SX.sym('vC');
    xL = SX.sym('xL');     
    vL = SX.sym('vL');
    uC = SX.sym('uC');
    uL = SX.sym('uL');
    theta = SX.sym('theta');
    omega = SX.sym('omega');
    q = SX.sym('q'); % a quadrature state
    x = vertcat(xC, vC, xL, vL, uC, uL, theta, omega, q);
    xdot = SX.sym('xdot', size(x));

    uCR = SX.sym('uCR');  % Controls
    uLR = SX.sym('uLR');
    u = vertcat(uCR, uLR);

    z = SX.sym('z',2); % define an algebraic state;

    
    %% Dynamics: implicit DAE formulation (index-1)
    f_expl = vertcat(vC, ...
                  - 1/tau1 * (vC - a1 * uC), ...
                  vL,...
                  - 1/tau2 * (vL - a2 * uL), ...
                  uCR,...
                  uLR,...
                  omega, ...
                  - (a1 * uCR * cos(theta) + g* sin(theta) + 2*vL*omega) / xL, ...
                  uCR^2 + xL^2 - z(1) + cos(xL)); % dynamics of quadrature state x2;
              
    f_impl = [(f_expl - xdot);
            z - [((theta^2)/8 + xL + 8 * q + sin(uLR) );
            cos(omega + 0.1) + (xdot(3) - uCR*vL)^2 - z(2)]];

    model.f_impl_expr = f_impl;
    model.x = x;
    model.xdot = xdot;
    model.u = u;
    model.z = z;
    model.name = model_name_prefix;
    
end