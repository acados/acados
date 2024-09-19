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



function model = simple_dae_model()
    % This function generates an implicit ODE / index-1 DAE model,
    % that depends on the symbolic CasADi variables x, xdot, u, z.

    import casadi.*

    %% Set up states & controls
    x1 = SX.sym('x1');     % Differential States
    x2 = SX.sym('x2');
    x = vertcat(x1, x2);

    z1 = SX.sym('z1');     % Algebraic states
    z2 = SX.sym('z2');
    z = vertcat(z1, z2);

    u1 = SX.sym('u1');     % Controls
    u2 = SX.sym('u2');
    u = vertcat(u1, u2);

    %% xdot
    x1_dot = SX.sym('x1_dot');     % Differential States
    x2_dot = SX.sym('x2_dot');
    xdot = vertcat(x1_dot, x2_dot);

    %% Dynamics: implicit DAE formulation (index-1)
    f_impl_expr = ...
        vertcat(x1_dot-0.1*x1+0.1*z2-u1, ...
                x2_dot+x2+0.01*z1-u2,  ...
                z1-x1, ...
                z2-x2);

    model = AcadosModel();
    model.x = x;
    model.xdot = xdot;
    model.u = u;
    model.z = z;

    model.f_impl_expr = f_impl_expr;
    model.name = 'simple_dae';
end

