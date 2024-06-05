#
# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# The 2-Clause BSD License
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice,
# this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation
# and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.;
#

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

