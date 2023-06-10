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

from acados_template import AcadosModel
from casadi import SX, vertcat, sin, cos, norm_2

import numpy as np

def export_disturbed_chain_mass_model(n_mass, m, D, L):

    model_name = 'chain_mass_ds_' + str(n_mass)
    x0 = np.array([0, 0, 0]) # fix mass (at wall)

    M = n_mass - 2 # number of intermediate massesu

    nx = (2*M + 1)*3  # differential states
    nu = 3            # control inputs

    xpos = SX.sym('xpos', (M+1)*3, 1) # position of fix mass eliminated
    xvel = SX.sym('xvel', M*3, 1)
    u = SX.sym('u', nu, 1)
    xdot = SX.sym('xdot', nx, 1)
    w = SX.sym('w', M*3, 1)

    f = SX.zeros(3*M, 1) # force on intermediate masses

    for i in range(M):
        f[3*i+2] = - 9.81

    for i in range(M+1):
        if i == 0:
            dist = xpos[i*3:(i+1)*3] - x0
        else:
            dist = xpos[i*3:(i+1)*3] - xpos[(i-1)*3:i*3]

        scale = D/m*(1-L/ norm_2(dist))
        F = scale*dist
    
        # mass on the right
        if i < M:
            f[i*3:(i+1)*3] -= F
    
        # mass on the left
        if i > 0:
            f[(i-1)*3:i*3] += F

    x = vertcat(xpos, xvel)

    # dynamics
    f_expl = vertcat(xvel, u, f+w)
    f_impl = xdot - f_expl

    model = AcadosModel()

    model.f_impl_expr = f_impl
    model.f_expl_expr = f_expl
    model.x = x
    model.xdot = xdot
    model.u = u
    model.p = w
    model.name = model_name

    return model

