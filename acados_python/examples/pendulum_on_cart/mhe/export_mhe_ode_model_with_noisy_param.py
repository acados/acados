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
from casadi import SX, vertcat, cos, sin

def export_mhe_ode_model_with_noisy_param():
    '''
    Export ode model augmented with an additional state corresponding to the
    parameter l, which is identified online
    '''

    model_name = 'mhe_pendulum_ode_with_noisy_param'

    # constants
    M = 1. # mass of the cart [kg]
    m = 0.1 # mass of the ball [kg]
    g = 9.81  # gravity constant [m/s^2]
    # l = 0.8 # length of the rod [m]  -> now estimated


    # set up states
    x1      = SX.sym('x1')
    v1      = SX.sym('v1')
    theta   = SX.sym('theta')
    dtheta  = SX.sym('dtheta')
    # add parameter l as state
    l       = SX.sym('l')

    x = vertcat(x1, theta, v1, dtheta, l)

    # state noise
    w_x1      = SX.sym('w_x1')
    w_v1      = SX.sym('w_v1')
    w_theta   = SX.sym('w_theta')
    w_dtheta  = SX.sym('w_dtheta')
    w_l       = SX.sym('w_l')

    w = vertcat(w_x1, w_theta, w_v1, w_dtheta, w_l)

    # xdot
    x1_dot      = SX.sym('x1_dot')
    theta_dot   = SX.sym('theta_dot')
    v1_dot      = SX.sym('v1_dot')
    dtheta_dot  = SX.sym('dtheta_dot')
    l_dot       = SX.sym('l_dot')

    xdot = vertcat(x1_dot, theta_dot, v1_dot, dtheta_dot, l_dot)

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
                     (-m*l*cos(theta)*sin(theta)*dtheta*dtheta + F*cos(theta)+(M+m)*g*sin(theta))/(l*denominator),
                     0)

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

