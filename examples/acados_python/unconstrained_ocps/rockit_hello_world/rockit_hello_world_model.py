#
# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# Licensed under the 2-Clause BSD License.

from acados_template import AcadosModel
from casadi import SX, vertcat

def export_rockit_hello_world_model() -> AcadosModel:

    model_name = 'rockit_hello_world_model'

    # set up states & controls
    x1 = SX.sym('x1')
    x2 = SX.sym('x2')

    x = vertcat(x1, x2)

    u = SX.sym('u')

    e = 1 - x2**2

    # xdot
    x1_dot = SX.sym('x1_dot')
    x2_dot = SX.sym('x2_dot')

    xdot = vertcat(x1_dot, x2_dot)

    # dynamics
    f_expl = vertcat(e * x1 - x2 + u,
                     x1)

    f_impl = xdot - f_expl

    model = AcadosModel()

    model.f_impl_expr = f_impl
    model.f_expl_expr = f_expl
    model.x = x
    model.xdot = xdot
    model.u = u
    model.name = model_name

    return model

