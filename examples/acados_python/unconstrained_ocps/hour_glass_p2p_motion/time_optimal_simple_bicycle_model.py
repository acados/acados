#
# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# Licensed under the 2-Clause BSD License.

from acados_template import AcadosModel
from casadi import SX, vertcat, sin, cos, tan

def export_time_optimal_simple_bicycle() -> AcadosModel:

    model_name = 'time_optimal_simple_bicycle'

    # constants
    l = 1.0

    # set up states & controls
    T = SX.sym('T')
    x = SX.sym('x')
    y = SX.sym('y')
    theta = SX.sym('theta')

    states = vertcat(T, x, y, theta)

    delta = SX.sym('delta')
    v = SX.sym('v')

    controls = vertcat(delta, v)

    # xdot
    T_dot = SX.sym('T_dot')
    x_dot = SX.sym('x_dot')
    y_dot = SX.sym('y_dot')
    theta_dot = SX.sym('theta_dot')

    states_dot = vertcat(T_dot, x_dot, y_dot, theta_dot)

    # dynamics
    f_expl = vertcat(0,
                     T*v*cos(theta),
                     T*v*sin(theta),
                     T*v/l*tan(delta)
                    )

    f_impl = states_dot - f_expl

    model = AcadosModel()

    model.f_impl_expr = f_impl
    model.f_expl_expr = f_expl
    model.x = states
    model.xdot = states_dot
    model.u = controls
    model.name = model_name

    return model

