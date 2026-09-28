%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.



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

