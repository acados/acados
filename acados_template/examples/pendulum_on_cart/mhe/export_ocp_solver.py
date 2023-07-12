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
import scipy
from acados_template import AcadosOcpSolver, AcadosOcp
from casadi import vertcat

def export_ocp_solver(model, N, h, Q, R, Fmax=80, use_cython=False):

    ocp = AcadosOcp()

    # set model
    ocp.model = model

    Tf = N*h
    nx = model.x.size()[0]
    nu = model.u.size()[0]
    ny = nx + nu
    ny_e = nx

    ocp.dims.N = N

    # set cost
    ocp.cost.cost_type = 'NONLINEAR_LS'
    ocp.cost.cost_type_e = 'NONLINEAR_LS'

    ocp.cost.W = scipy.linalg.block_diag(R, Q)

    ocp.cost.W_e = Q

    x = ocp.model.x
    u = ocp.model.u

    ocp.model.cost_y_expr = vertcat(u, x)

    ocp.model.cost_y_expr_e = x

    ocp.cost.yref  = np.zeros((ny, ))
    ocp.cost.yref_e = np.zeros((ny_e, ))

    # setting bounds
    ocp.constraints.lbu = np.array([-Fmax])
    ocp.constraints.ubu = np.array([+Fmax])
    ocp.constraints.x0 = np.array([0.0, np.pi, 0.0, 0.0])
    ocp.constraints.idxbu = np.array([0])

    # set QP solver
    # ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM'
    ocp.solver_options.qp_solver = 'FULL_CONDENSING_QPOASES'
    ocp.solver_options.hessian_approx = 'GAUSS_NEWTON'
    ocp.solver_options.integrator_type = 'ERK'

    # set prediction horizon
    ocp.solver_options.tf = Tf
    ocp.solver_options.nlp_solver_type = 'SQP'
    # ocp.solver_options.nlp_solver_type = 'SQP_RTI'
    ocp.solver_options.nlp_solver_max_iter = 200

    if use_cython:
        AcadosOcpSolver.generate(ocp, json_file='acados_ocp.json')
        AcadosOcpSolver.build(ocp.code_export_directory, with_cython=True)
        acados_solver = AcadosOcpSolver.create_cython_solver('acados_ocp.json')
    else:
        acados_solver = AcadosOcpSolver(ocp, json_file='acados_ocp.json')

    return acados_solver