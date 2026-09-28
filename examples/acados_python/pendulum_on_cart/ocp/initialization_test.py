#
# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# Licensed under the 2-Clause BSD License.

import sys
sys.path.insert(0, '../common')

from acados_template import AcadosOcp, AcadosOcpSolver
from pendulum_model import export_pendulum_ode_model
import numpy as np
from utils import plot_pendulum

import casadi as ca


def main(qp_solver: str = 'PARTIAL_CONDENSING_HPIPM'):
    # create ocp object to formulate the OCP
    ocp = AcadosOcp()

    # set model
    model = export_pendulum_ode_model()
    model.name = f'test_{qp_solver.replace("CONDENSING", "cond").lower()}'
    ocp.model = model

    Tf = 1.0
    nx = model.x.rows()
    nu = model.u.rows()
    N = 20

    # set prediction horizon
    ocp.solver_options.N_horizon = N
    ocp.solver_options.tf = Tf

    # cost matrices
    Q_mat = 2*np.diag([1e3, 1e3, 1e-2, 1e-2])
    R_mat = 2*np.diag([1e-2])

    # path cost
    ocp.cost.cost_type = 'NONLINEAR_LS'
    ocp.model.cost_y_expr = ca.vertcat(model.x, model.u)
    ocp.cost.yref = np.zeros((nx+nu,))
    ocp.cost.W = ca.diagcat(Q_mat, R_mat).full()

    # terminal cost
    ocp.cost.cost_type_e = 'NONLINEAR_LS'
    ocp.cost.yref_e = np.zeros((nx,))
    ocp.model.cost_y_expr_e = model.x
    ocp.cost.W_e = Q_mat

    # set constraints
    Fmax = 20
    ocp.constraints.lbu = np.array([-Fmax])
    ocp.constraints.ubu = np.array([+Fmax])
    ocp.constraints.idxbu = np.array([0])

    ocp.constraints.x0 = np.array([0.0, np.pi, 0.0, 0.0])

    # set options
    ocp.solver_options.qp_solver = qp_solver
    ocp.solver_options.hessian_approx = 'GAUSS_NEWTON'
    ocp.solver_options.integrator_type = 'IRK'
    ocp.solver_options.nlp_solver_type = 'SQP'
    ocp.solver_options.globalization = 'MERIT_BACKTRACKING'
    ocp.solver_options.qp_tol = 1e-9
    ocp.solver_options.qp_solver_mu0 = 1e2

    ocp_solver = AcadosOcpSolver(ocp, verbose=True)

    status = ocp_solver.solve()
    ocp_solver.print_statistics()

    if status != 0:
        raise Exception(f'acados returned status {status}.')

    # reference solve
    sol = ocp_solver.get_iterate()
    qp_iter = ocp_solver.get_stats("qp_iter")
    nlp_iter = ocp_solver.get_stats("nlp_iter")
    assert nlp_iter == 10, f"cold start should require 10 iterations, got {nlp_iter}"

    ocp_solver.set_iterate(sol)
    ocp_solver.solve()
    ocp_solver.print_statistics()
    nlp_iter = ocp_solver.get_stats("nlp_iter")
    assert nlp_iter == 0, f"hot start should require 0 iterations, got {nlp_iter}"

    disturbed_sol = sol
    disturbed_sol.x[0] += 0.1 * disturbed_sol.x[0]
    ocp_solver.set_iterate(disturbed_sol)
    status = ocp_solver.solve()
    ocp_solver.print_statistics()
    nlp_iter = ocp_solver.get_stats("nlp_iter")
    assert nlp_iter == 7, f"warm start should require 7 iterations, got {nlp_iter}"

    # QP warm start
    ocp_solver.set_iterate(disturbed_sol)
    ocp_solver.options_set("qp_warm_start", 2)
    ocp_solver.options_set("warm_start_first_qp_from_nlp", False)
    ocp_solver.options_set("warm_start_first_qp", True)
    # ocp_solver.options_set("qp_mu0", 1e-3)
    status = ocp_solver.solve()
    ocp_solver.print_statistics()

    # cold
    for warm_start in [0, 1, 2, 3]:
        for t0_init in [0, 1, 2]:
            print(f"testing with {warm_start=}, {t0_init=}")
            ocp_solver.set_iterate(disturbed_sol)
            ocp_solver.options_set("qp_warm_start", warm_start)
            ocp_solver.options_set("qp_t0_init", t0_init)
            status = ocp_solver.solve()
            ocp_solver.print_statistics()
            qp_iter = ocp_solver.get_stats('qp_iter')

            if warm_start in [0, 1]: # same in acados.
                if t0_init == 0:
                    assert sum(qp_iter) == 62, f"warm start should require 62 QP iterations, got {sum(qp_iter)} = sum({qp_iter})"
                elif t0_init == 1:
                    assert sum(qp_iter) == 70, f"warm start should require 70 QP iterations, got {sum(qp_iter)} = sum({qp_iter})"
                elif t0_init == 2:
                    assert sum(qp_iter) == 28, f"warm start should require 28 QP iterations, got {sum(qp_iter)} = sum({qp_iter})"
            # here t0_init doesnt matter
            elif warm_start == 2:
                assert sum(qp_iter) == 17, f"warm start should require 17 QP iterations, got {sum(qp_iter)} = sum({qp_iter})"
            elif warm_start == 3:
                assert sum(qp_iter) == 10, f"warm start should require 10 QP iterations, got {sum(qp_iter)} = sum({qp_iter})"
    del ocp_solver


if __name__ == '__main__':
    main('FULL_CONDENSING_HPIPM')
    main('PARTIAL_CONDENSING_HPIPM')
