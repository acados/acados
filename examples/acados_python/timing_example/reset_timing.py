#
# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# Licensed under the 2-Clause BSD License.

import sys
sys.path.insert(0, '../pendulum_on_cart/common')

from acados_template import AcadosOcp, AcadosOcpSolver
from pendulum_model import export_pendulum_ode_model
import numpy as np
import scipy.linalg
from utils import plot_pendulum
from casadi import vertcat
import time

def main(cost_type='NONLINEAR_LS', hessian_approximation='EXACT', ext_cost_use_num_hess=0,
         integrator_type='ERK'):
    print(f"using: cost_type {cost_type}, integrator_type {integrator_type}")
    # create ocp object to formulate the OCP
    ocp = AcadosOcp()

    # set model
    model = export_pendulum_ode_model()
    ocp.model = model

    Tf = 1.0
    nx = model.x.rows()
    nu = model.u.rows()
    ny = nx + nu
    ny_e = nx
    N = 20

    ocp.solver_options.N_horizon = N

    # set cost
    Q = 2*np.diag([1e3, 1e3, 1e-2, 1e-2])
    R = 2*np.diag([1e-2])

    x = ocp.model.x
    u = ocp.model.u

    cost_W = scipy.linalg.block_diag(Q, R)

    if cost_type == 'LS':
        ocp.cost.cost_type = 'LINEAR_LS'
        ocp.cost.cost_type_e = 'LINEAR_LS'

        ocp.cost.Vx = np.zeros((ny, nx))
        ocp.cost.Vx[:nx,:nx] = np.eye(nx)

        Vu = np.zeros((ny, nu))
        Vu[4,0] = 1.0
        ocp.cost.Vu = Vu

        ocp.cost.Vx_e = np.eye(nx)

    elif cost_type == 'NONLINEAR_LS':
        ocp.cost.cost_type = 'NONLINEAR_LS'
        ocp.cost.cost_type_e = 'NONLINEAR_LS'

        ocp.model.cost_y_expr = vertcat(x, u)
        ocp.model.cost_y_expr_e = x

    elif cost_type == 'EXTERNAL':
        ocp.cost.cost_type = 'EXTERNAL'
        ocp.cost.cost_type_e = 'EXTERNAL'

        ocp.model.cost_expr_ext_cost = vertcat(x, u).T @ cost_W @ vertcat(x, u)
        ocp.model.cost_expr_ext_cost_e = x.T @ Q @ x

    else:
        raise Exception('Unknown cost_type. Possible values are \'LS\' and \'NONLINEAR_LS\'.')

    if cost_type in ['LS', 'NONLINEAR_LS']:
        ocp.cost.yref = np.zeros((ny, ))
        ocp.cost.yref_e = np.zeros((ny_e, ))
        ocp.cost.W_e = Q
        ocp.cost.W = cost_W

    # set constraints
    Fmax = 80
    ocp.constraints.lbu = np.array([-Fmax])
    ocp.constraints.ubu = np.array([+Fmax])
    x0 = np.array([0.0, np.pi, 0.0, 0.0])
    ocp.constraints.x0 = x0
    ocp.constraints.idxbu = np.array([0])

    ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM' # FULL_CONDENSING_QPOASES
    ocp.solver_options.hessian_approx = hessian_approximation
    ocp.solver_options.regularize_method = 'CONVEXIFY'
    ocp.solver_options.integrator_type = integrator_type
    ocp.solver_options.qp_solver_cond_N = 5

    # set prediction horizon
    ocp.solver_options.tf = Tf
    ocp.solver_options.nlp_solver_type = 'SQP' # SQP_RTI
    ocp.solver_options.ext_cost_num_hess = ext_cost_use_num_hess

    ocp.code_gen_options.code_export_directory = 'c_generated_code'

    # create solver
    ocp_solver = AcadosOcpSolver.create_cython_solver(ocp)

    # time create
    Ncreate = 1
    create_time = []
    for i in range(Ncreate):
        t0 = time.time()
        # ocp_solver = AcadosOcpSolver(ocp, json_file = 'acados_ocp.json', build=False, generate=False)
        # NOTE: this needs to match the code_export_directory set above
        from c_generated_code.acados_ocp_solver_pyx import AcadosOcpSolverCython
        ocp_solver = AcadosOcpSolverCython(ocp.model.name, ocp.solver_options.nlp_solver_type, ocp.solver_options.N_horizon)
        create_time.append( time.time() - t0)
    print(f"create_time: min {np.min(create_time)*1e3:.4f} ms mean {np.mean(create_time)*1e3:.4f} ms max {np.max(create_time)*1e3:.4f} ms over {Ncreate} executions")
    # create_time: min 0.2189 ms mean 0.2586 ms max 0.6881 ms over 10000 executions

    # RESET
    Nreset = 3
    reset_time = []
    for i in range(Nreset):
        # set NaNs as input to test reset() -> NOT RECOMMENDED!!!
        # ocp_solver.options_set('print_level', 2)
        for i in range(N):
            ocp_solver.set(i, 'x', np.nan * np.ones((nx,)))
            ocp_solver.set(i, 'u', np.nan * np.ones((nu,)))
        status = ocp_solver.solve()
        ocp_solver.print_statistics() # encapsulates: stat = ocp_solver.get_stats("statistics")
        if status == 0:
            raise Exception(f'acados returned status {status}, although NaNs were given.')
        else:
            print(f'acados returned status {status}, which is expected, since NaNs were given.')

        t0 = time.time()
        ocp_solver.reset()
        reset_time.append( time.time() - t0)
    print(f"reset_time: min {np.min(reset_time)*1e3:.4f} ms mean {np.mean(reset_time)*1e3:.4f} ms max {np.max(reset_time)*1e3:.4f} ms over {Nreset} executions")

    # OLD: without HPIPM mem reset:
    # reset_time: min 0.0019 ms mean 0.0022 ms max 0.0250 ms over 10000 executions

    # NEW: with HPIPM mem reset:
    # reset_time: min 0.0036 ms mean 0.0047 ms max 0.0906 ms over 10000 executions

    for i in range(N):
        ocp_solver.set(i, 'x', x0)

    if cost_type == 'EXTERNAL':
        # NOTE: hessian is wrt [u,x]
        if ext_cost_use_num_hess:
            for i in range(N):
                ocp_solver.cost_set(i, "ext_cost_num_hess", np.diag([0.04, 4000, 4000, 0.04, 0.04, ]))
            ocp_solver.cost_set(N, "ext_cost_num_hess", np.diag([4000, 4000, 0.04, 0.04, ]))

    simX = np.zeros((N+1, nx))
    simU = np.zeros((N, nu))

    status = ocp_solver.solve()

    ocp_solver.print_statistics()
    if status != 0:
        raise Exception(f'acados returned status {status} for cost_type {cost_type}\n'
                        f'integrator_type = {integrator_type}.')

    # get solution
    iterate = ocp_solver.get_iterate()
    simX[:] = iterate.x
    simU[:] = iterate.u


if __name__ == '__main__':
    # for integrator_type in ['ERK', 'IRK']:
    # ['LIFTED_IRK']
    for integrator_type in ['ERK']:
        for cost_type in ['EXTERNAL']:
            # for cost_type in ['EXTERNAL', 'LS', 'NONLINEAR_LS']:
            hessian_approximation = 'GAUSS_NEWTON' # 'GAUSS_NEWTON, EXACT
            ext_cost_use_num_hess = 1
            main(cost_type=cost_type, hessian_approximation=hessian_approximation,
                ext_cost_use_num_hess=ext_cost_use_num_hess, integrator_type=integrator_type)
