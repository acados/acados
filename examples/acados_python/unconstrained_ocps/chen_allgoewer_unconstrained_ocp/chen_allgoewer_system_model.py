#
# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# Licensed under the 2-Clause BSD License.

from acados_template import AcadosModel
from casadi import SX, MX, vertcat

def export_chen_allgoewer_model(use_SX=True) -> AcadosModel:

    model_name = 'chen_allgoewer_ode'

    # constants
    mu = 0.7

    # set up states & controls
    if use_SX:
        x1 = SX.sym('x1')
        x2 = SX.sym('x2')

        u = SX.sym('u')
        # xdot
        x1_dot = SX.sym('x1_dot')
        x2_dot = SX.sym('x2_dot')
    else:
        x1 = MX.sym('x1')
        x2 = MX.sym('x2')

        u = MX.sym('u')
        # xdot
        x1_dot = MX.sym('x1_dot')
        x2_dot = MX.sym('x2_dot')

    x = vertcat(x1, x2)
    xdot = vertcat(x1_dot, x2_dot)

    # dynamics
    f_expl = vertcat(x2 + u*(mu + (1-mu)*x2),
                     x1 + u*(mu-4*(1-mu)*x2))

    f_impl = xdot - f_expl

    model = AcadosModel()

    model.f_impl_expr = f_impl
    model.f_expl_expr = f_expl
    model.x = x
    model.xdot = xdot
    model.u = u
    model.name = model_name

    return model

