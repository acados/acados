%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.


function model = get_double_integrator_model()
    import casadi.*
    % set up states & controls
    p = SX.sym('p');
    v = SX.sym('v');
    x = vertcat(p, v);
    nx = length(x);

    u = SX.sym('u');

    % set up dynamics
    f_expl = vertcat(v, u);
    xdot = SX.sym('xdot', nx);
    f_impl = f_expl - xdot;

    % set up model
    model = AcadosModel();
    model.name = 'double_integrator';
    model.x = x;
    model.xdot = xdot;
    model.u = u;
    model.f_impl_expr = f_impl;
    model.f_expl_expr = f_expl;
end

