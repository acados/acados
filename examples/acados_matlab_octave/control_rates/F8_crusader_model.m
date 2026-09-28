%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.

%

function model = F8_crusader_model()
    import casadi.*

    %% System parameters
    a1 = -0.877;
    a2 = -0.088;
    a3 = 0.47;
    a4 = -0.019;
    a5 = 3.846;
    a6 = -0.215;
    a7 = 0.28;
    a8 = 0.47;
    a9 = 0.63;
    c1 = -4.208;
    c2 = -0.396;
    c3 = -0.47;
    c4 = -3.564;
    c5 = -20.967;
    c6 = 6.265;
    c7 = 46;
    c8 = 61.4;

    %% Symbolic variables
    % states
    x1 = SX.sym('x1');
    x2 = SX.sym('x2');
    x3 = SX.sym('x3');
    x4 = SX.sym('x4');  % extra state for the original control input
    x = vertcat(x1, x2, x3, x4);

    % controls
    u = SX.sym('u');  % input rate becomes the input

    % state derivatives
    xdot = SX.sym('xdot', length(x), 1);

    %% The system dynamics
    % Explicit Dynamics
    f_expl_expr = vertcat(...
        a1*x1+x3+a2*x1*x3+a3*x1.^2+a4*x2.^2-x1.^2*x3+a5*x1.^3+a6*x4+a7*x1.^2*x4+a8*x1*x4.^2+a9*x4.^3, ...
        x3, ...
        c1*x1+c2*x3+c3*x1.^2+c4*x1.^3+c5*x4+c6*x1.^2*x4+c7*x1*x4.^2+c8*x4.^3, ...
        u);  % dynamics of the added state

    % Implicit Dynamics
    f_impl_expr = f_expl_expr - xdot;

    % setup an acados model
    model = AcadosModel();
    model.x = x;
    model.u = u;
    model.xdot = xdot;
    model.f_expl_expr = f_expl_expr;
    model.f_impl_expr = f_impl_expr;
    model.name = 'F8_crusader_model';

end