#
# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# Licensed under the 2-Clause BSD License.

from acados_template import AcadosModel
from casadi import SX, vertcat, cos, sin

def export_mhe_ode_model() -> AcadosModel:

    model_name = 'mhe_pendulum_ode'

    # constants
    M = 1.
    m = 0.1
    g = 9.81
    l = 0.8

    # set up states
    x1      = SX.sym('x1')
    v1      = SX.sym('v1')
    theta   = SX.sym('theta')
    dtheta  = SX.sym('dtheta')

    x = vertcat(x1, theta, v1, dtheta)

    # state noise
    w_x1      = SX.sym('w_x1')
    w_v1      = SX.sym('w_v1')
    w_theta   = SX.sym('w_theta')
    w_dtheta  = SX.sym('w_dtheta')

    w = vertcat(w_x1, w_theta, w_v1, w_dtheta)

    # xdot
    x1_dot      = SX.sym('x1_dot')
    theta_dot   = SX.sym('theta_dot')
    v1_dot      = SX.sym('v1_dot')
    dtheta_dot  = SX.sym('dtheta_dot')

    xdot = vertcat(x1_dot, theta_dot, v1_dot, dtheta_dot)

    # algebraic variables
    z = []

    # parameters <= controls
    F = SX.sym('F')
    p = F

    # dynamics
    denominator = M + m - m*cos(theta)*cos(theta)
    f_expl = vertcat(v1,
                     dtheta,
                     (-m*l*sin(theta)*dtheta*dtheta + m*g*cos(theta)*sin(theta)+F)/denominator,
                     (-m*l*cos(theta)*sin(theta)*dtheta*dtheta + F*cos(theta)+(M+m)*g*sin(theta))/(l*denominator))

    # add additive state noise
    f_expl = f_expl + w
    f_impl = xdot - f_expl

    model = AcadosModel()

    model.f_impl_expr = f_impl
    model.f_expl_expr = f_expl
    model.x = x
    model.xdot = xdot
    model.u = w
    model.z = z
    model.p = p
    model.name = model_name

    return model

