#
# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# Licensed under the 2-Clause BSD License.

from acados_template import AcadosModel
from casadi import MX, external
import sys,os
from ctypes import *


def export_external_ode_model():

    model_name = 'external_ode'

    # Declare model variables
    x = MX.sym('x', 2)
    u = MX.sym('u', 1)
    xDot = MX.sym('xDot', 2)
    cdll.LoadLibrary('./test_external_lib/build/libexternal_ode_casadi.so')
    f_ext = external('libexternal_ode_casadi', 'libexternal_ode_casadi.so', {'enable_fd': True})
    f_expl = f_ext(x, u)

    f_impl = xDot - f_expl

    model = AcadosModel()

    model.f_impl_expr = f_impl
    model.f_expl_expr = f_expl
    model.x = x
    model.xdot = xDot
    model.u = u
    model.p =[]
    model.name = model_name

    return model 

