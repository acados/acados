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

from enum import Enum
import numpy as np

from acados_template import AcadosModel
from casadi import SX, vertcat, cos, sin, sqrt, sumsqr

from mpc_parameters import MPCParam


class RobotState(Enum):
    POSX = 0
    POSY = 1
    THETA = 2
    VEL = 3
    OMEGA = 4


def dxdt(x0, u):
    dx = vertcat(x0[3]*cos(x0[2]),
                 x0[3]*sin(x0[2]),
                 x0[4],
                 u[0],
                 u[1])
    return dx


def export_diff_drive_model(as_cfg:MPCParam):

    model_name = 'diff_drive_model'

    nx = as_cfg.nx
    nu = as_cfg.nu
    num_obs = as_cfg.num_obs

    x = SX.sym('x', nx, 1)
    u = SX.sym('u', nu, 1)
    xdot = SX.sym('xdot', nx, 1)

    # dynamics
    f_expl = vertcat(x[3]*cos(x[2]),
                     x[3]*sin(x[2]),
                     x[4],
                     u[0],
                     u[1])

    f_impl = xdot - f_expl

    model = AcadosModel()

    model.f_expl_expr = f_expl
    model.f_impl_expr = f_impl
    model.x = x
    model.xdot = xdot
    model.u = u
    model.name = model_name

    con_expr = SX.sym('con_expr', num_obs, 1)
    obs_param = SX.sym('obs', 2*num_obs, 1)
    # [x0, y0, x1, y1, ...]
    p_obs = [obs_param[2*k:2*(k+1)] for k in range(0, num_obs)]

    for i in range(0, num_obs):
        dist = x[0:2] - p_obs[i]
        con_expr[i] = sqrt(sumsqr(dist))
    model.p = obs_param
    model.con_h_expr = con_expr
    model.con_h_expr_e = con_expr

    return model