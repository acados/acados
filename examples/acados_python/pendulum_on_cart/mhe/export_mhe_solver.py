#
# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# Licensed under the 2-Clause BSD License.


import numpy as np
from scipy.linalg import block_diag
from acados_template import AcadosModel, AcadosOcp, AcadosOcpSolver
from casadi import vertcat


def export_mhe_solver(model: AcadosModel, N: int, h, Q, Q0, R) -> AcadosOcpSolver:

    ocp_mhe = AcadosOcp()

    ocp_mhe.model = model

    x = ocp_mhe.model.x
    u = ocp_mhe.model.u

    nx = x.rows()
    nu = u.rows()
    nparam = model.p.rows()

    ny_0 = 3*nx   # h(x), w and arrival cost
    ny = 2*nx     # h(x), w

    ocp_mhe.solver_options.N_horizon = N

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
    ocp_mhe.solver_options.cost_scaling = np.ones((N+1,)) # we do not want to scale with the time step

    # set prediction horizon
    ocp_mhe.solver_options.tf = N*h

    ocp_mhe.solver_options.nlp_solver_type = 'SQP'
    # ocp_mhe.solver_options.nlp_solver_type = 'SQP_RTI'
    ocp_mhe.solver_options.nlp_solver_max_iter = 200

    acados_solver_mhe = AcadosOcpSolver(ocp_mhe, json_file = 'acados_ocp_mhe.json')

    return acados_solver_mhe
