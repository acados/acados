#
# Copyright 2019 Gianluca Frison, Dimitris Kouzoupis, Robin Verschueren,
# Andrea Zanelli, Niels van Duijkeren, Jonathan Frey, Tommaso Sartor,
# Branimir Novoselnik, Rien Quirynen, Rezart Qelibari, Dang Doan,
# Jonas Koenemann, Yutao Chen, Tobias Schöls, Jonas Schlagenhauf, Moritz Diehl
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
from acados_template import *


def export_mhe_solver(model, N, h, Q, Q0, R):

    # create render arguments
    ocp_mhe = AcadosOcp()

    ocp_mhe.model = model

    nx = model.x.size()[0]
    nu = model.u.size()[0]
    nparam = model.p.size()[0]

    ny = 3*nx     # h(x), w and arrival cost
    ny_e = 0

    # set ocp_nlp_dimensions
    ocp_mhe.dims.nx  = nx
    ocp_mhe.dims.ny  = ny
    ocp_mhe.dims.np  = nparam
    ocp_mhe.dims.ny_e = ny_e
    ocp_mhe.dims.nbx = 0
    ocp_mhe.dims.nbu = 0  # nu
    ocp_mhe.dims.nu  = model.u.size()[0]
    ocp_mhe.dims.N   = N

    ocp_mhe.dims.nbx_0 = 0

    # set weighting matrices

    ocp_mhe.cost.cost_type = 'NONLINEAR_LS'
    ocp_mhe.cost.cost_type_e = 'LINEAR_LS'

    ocp_mhe.cost.W = block_diag(R, Q, np.zeros((nx, nx)))

    x = ocp_mhe.model.x
    u = ocp_mhe.model.u

    ocp_mhe.model.cost_y_expr = vertcat(x, u, x)

    ocp_mhe.parameter_values = np.zeros((nparam, ))

    # set y_ref for all stages
    ocp_mhe.cost.yref  = np.zeros((ny,))
    ocp_mhe.cost.yref_e = np.zeros((ny_e, ))

    # set QP solver
    # ocp_mhe.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM'
    ocp_mhe.solver_options.qp_solver = 'FULL_CONDENSING_QPOASES'
    ocp_mhe.solver_options.hessian_approx = 'GAUSS_NEWTON'
    ocp_mhe.solver_options.integrator_type = 'ERK'

    # set prediction horizon
    ocp_mhe.solver_options.tf = N*h

    ocp_mhe.solver_options.nlp_solver_type = 'SQP'
    # ocp_mhe.solver_options.nlp_solver_type = 'SQP_RTI'
    ocp_mhe.solver_options.nlp_solver_max_iter = 200

    acados_solver_mhe = AcadosOcpSolver(ocp_mhe, json_file = 'acados_ocp.json')

    # set arrival cost weighting matrix
    acados_solver_mhe.cost_set(0, "W", block_diag(R, Q, Q0))

    return acados_solver_mhe
