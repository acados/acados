#
# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# Licensed under the 2-Clause BSD License.

import numpy as np
import scipy
from acados_template import AcadosOcpSolver, AcadosOcp
from casadi import vertcat

def export_ocp_solver(model, N, h, Q, R, Fmax=80, use_cython=False):

    ocp = AcadosOcp()

    # set model
    ocp.model = model

    Tf = N*h
    nx = model.x.rows()
    nu = model.u.rows()
    ny = nx + nu
    ny_e = nx

    ocp.solver_options.N_horizon = N

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
        acados_solver = AcadosOcpSolver.create_cython_solver(ocp)
    else:
        acados_solver = AcadosOcpSolver(ocp)

    return acados_solver