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
from utils import plot_pendulum
import numpy as np
import casadi as ca
import scipy.linalg

from dataclasses import dataclass

@dataclass
class GlobalizationOptions:
    globalization: str
    alpha_min: float
    use_SOC: bool

def main(globalization_options: GlobalizationOptions):
    ocp = AcadosOcp()

    # set model
    model = export_pendulum_ode_model()
    ocp.model = model

    Tf = 1.0
    nx = model.x.rows()
    nu = model.u.rows()

    N = 20

    # set dimensions
    ocp.solver_options.N_horizon = N

    # set cost
    Q_mat = 2*np.diag([1e3, 1e3, 1e-2, 1e-2])
    R_mat = 2*np.diag([1e-2])

    ocp.cost.cost_type = 'EXTERNAL'
    ocp.cost.cost_type_e = 'EXTERNAL'

    cost_W = scipy.linalg.block_diag(Q_mat, R_mat)

    x = model.x
    u = model.u
    ocp.model.cost_expr_ext_cost = .5*ca.vertcat(x, u).T @ cost_W @ ca.vertcat(x, u)
    ocp.model.cost_expr_ext_cost_e = .5*x.T @ Q_mat @ x

    # add log barrier term
    p_max = 1.0
    ocp.model.cost_expr_ext_cost += -ca.log(ca.fmax(p_max - x[0], 0))

    # set constraints
    Fmax = 80
    ocp.constraints.lbu = np.array([-Fmax])
    ocp.constraints.ubu = np.array([+Fmax])
    ocp.constraints.idxbu = np.array([0])

    ocp.constraints.x0 = np.array([0.0, np.pi, 0.0, 0.0])

    # set options
    ocp.solver_options.globalization = globalization_options.globalization
    ocp.solver_options.globalization_alpha_min = globalization_options.alpha_min
    ocp.solver_options.globalization_use_SOC = globalization_options.use_SOC
    ocp.solver_options.print_level = 0

    ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM'
    ocp.solver_options.hessian_approx = 'GAUSS_NEWTON'
    ocp.solver_options.integrator_type = 'ERK'
    ocp.solver_options.nlp_solver_type = 'SQP'
    ocp.solver_options.nlp_solver_max_iter = 1000
    ocp.solver_options.globalization_full_step_dual = 0
    ocp.solver_options.qp_solver_iter_max = 2000
    ocp.solver_options.tol = 1e-4
    ocp.solver_options.qp_tol = 1e-6


    ocp.solver_options.tf = Tf

    # create solver
    ocp_solver = AcadosOcpSolver(ocp)

    simX = np.zeros((N+1, nx))
    simU = np.zeros((N, nu))

    status = ocp_solver.solve()

    ocp_solver.print_statistics()

    # get solution
    iterate = ocp_solver.get_iterate()
    simX[:] = iterate.x
    simU[:] = iterate.u

    print("cost function value", ocp_solver.get_cost())

    if globalization_options.globalization == "MERIT_BACKTRACKING" and globalization_options.use_SOC:
        expected_status = 2
    elif globalization_options.globalization == "MERIT_BACKTRACKING" and not globalization_options.use_SOC:
        expected_status = 0
    else:
        print("Warning: no expected status defined for globalization options")
        expected_status = 0

    if status != expected_status:
        raise Exception(f"Test failed: expected status {expected_status}, got {status} for setting {globalization_options}")

    if (status != 1) and (np.isnan(simX).any() or np.isnan(simU).any()):
        raise Exception(f"Test failed: Solver returned NaN, but not NaN status for setting {globalization_options}")

    if any(simX[:, 0] > p_max):
        raise Exception(f"Test failed: Position constraint violated for setting {globalization_options}")

    print(f"Test passed: expected status {expected_status}, got {status} for setting {globalization_options}")

    plot_pendulum(np.linspace(0, Tf, N+1), Fmax, simU, simX, latexify=False)
    ocp_solver = None

if __name__ == "__main__":
    settings = [
        GlobalizationOptions(globalization="MERIT_BACKTRACKING", alpha_min=1e-4, use_SOC=False),
        GlobalizationOptions(globalization="MERIT_BACKTRACKING", alpha_min=1e-4, use_SOC=True),
        # Funnel returns error: Step size gets too small. Should enter penalty phase.
        # GlobalizationOptions(globalization="FUNNEL_L1PEN_LINESEARCH", alpha_min=1e-4, use_SOC=True),
    ]
    for setting in settings:
        main(setting)
