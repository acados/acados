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


import numpy as np
from scipy.linalg import block_diag
from acados_template import AcadosOcp, AcadosOcpSolver
from casadi import vertcat


def export_mhe_solver(model, N, h, Q, Q0, R):

    ocp_mhe = AcadosOcp()

    ocp_mhe.model = model

    x = ocp_mhe.model.x
    u = ocp_mhe.model.u

    nx = x.size()[0]
    nu = u.size()[0]
    nparam = model.p.size()[0]

    ny_0 = 3*nx     # h(x), w and arrival cost
    ny = 2*nx     # h(x), w
    ny_e = 0

    ocp_mhe.dims.N = N

    ## set cost
    ocp_mhe.cost.cost_type_0 = 'NONLINEAR_LS' # 'LINEAR_LS'

    if ocp_mhe.cost.cost_type_0 == 'LINEAR_LS':
        ocp_mhe.cost.W_0 = block_diag(R, Q, Q0)
        ocp_mhe.cost.Vx_0 = np.zeros((ny_0, nx))
        ocp_mhe.cost.Vx_0[:nx, :] = np.eye(nx)
        ocp_mhe.cost.Vx_0[2*nx:3*nx, :] = np.eye(nx)

        ocp_mhe.cost.Vu_0 = np.zeros((ny_0, nu))
        ocp_mhe.cost.Vu_0[1*nx:2*nx, :] = np.eye(nx)

        ocp_mhe.cost.yref_0 = np.zeros((ny_0,))

    elif ocp_mhe.cost.cost_type_0 == "NONLINEAR_LS":
        ocp_mhe.cost.W_0 = block_diag(R, Q, Q0)
        ocp_mhe.model.cost_y_expr_0 = vertcat(x, u, x)
        ocp_mhe.cost.yref_0 = np.zeros((ny_0,))
    else:
        Exception('Unknown cost type')

    # intermediate
    ocp_mhe.cost.cost_type = 'NONLINEAR_LS'

    ocp_mhe.cost.W = block_diag(R, Q)
    ocp_mhe.model.cost_y_expr = vertcat(x, u)
    ocp_mhe.parameter_values = np.zeros((nparam, ))
    ocp_mhe.cost.yref = np.zeros((ny,))

    # terminal
    ocp_mhe.cost.cost_type_e = 'LINEAR_LS'
    ocp_mhe.cost.yref_e = np.zeros((0, ))

    ocp_mhe.solver_options.qp_solver = 'FULL_CONDENSING_QPOASES'
    ocp_mhe.solver_options.hessian_approx = 'GAUSS_NEWTON'
    ocp_mhe.solver_options.integrator_type = 'ERK'

    # set prediction horizon
    ocp_mhe.solver_options.tf = N*h

    ocp_mhe.solver_options.nlp_solver_type = 'SQP'
    # ocp_mhe.solver_options.nlp_solver_type = 'SQP_RTI'
    ocp_mhe.solver_options.nlp_solver_max_iter = 200

    acados_solver_mhe = AcadosOcpSolver(ocp_mhe, json_file = 'acados_ocp_mhe.json')

    return acados_solver_mhe
