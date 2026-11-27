#
# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# Licensed under the 2-Clause BSD License.

from acados_template import AcadosOcp, AcadosOcpSolver, AcadosModel
import numpy as np
import casadi as ca
import pytest

# tests different globalization settings on simple test problem
#
# min x^2
# s.t. x \in [-10, 10]
#
# with wrong hessian: 1.0

# Settings
TOL = 1e-6

ARMIJO_SETTINGS = [
    pytest.param(globalization, sufficient_descent, use_soc,
                 id=f"{globalization}-descent_{sufficient_descent}-soc_{use_soc}")
    for globalization in ("MERIT_BACKTRACKING", "FIXED_STEP")
    for sufficient_descent in (0, 1)
    for use_soc in (0, 1)
    if globalization != "FIXED_STEP" or not (sufficient_descent or use_soc)
]


@pytest.mark.parametrize(
    "globalization,line_search_use_sufficient_descent,globalization_use_SOC",
    ARMIJO_SETTINGS,
)
def test_armijo(globalization, line_search_use_sufficient_descent, globalization_use_SOC):

    # create ocp object to formulate the OCP
    ocp = AcadosOcp()

    # set model
    model = AcadosModel()
    x = ca.SX.sym('x')

    # dynamics: identity
    model.disc_dyn_expr = x
    model.x = x
    model.u = ca.SX.sym('u', 0, 0) # [] / None doesnt work
    model.p = []
    model.name = f'armijo_problem'
    ocp.model = model

    # discretization
    Tf = 1
    N = 0
    ocp.solver_options.N_horizon = N
    ocp.solver_options.tf = Tf

    # cost
    ocp.cost.cost_type_e = 'EXTERNAL'
    ocp.model.cost_expr_ext_cost_e = x @ x
    ocp.model.cost_expr_ext_cost_custom_hess_e = ca.DM(1.0) # 2.0 is the actual hessian

    # constarints
    ocp.constraints.idxbx_e = np.array([0])
    ocp.constraints.lbx_e = np.array([-10.0])
    ocp.constraints.ubx_e = np.array([10.0])

    # options
    ocp.solver_options.qp_solver = 'FULL_CONDENSING_QPOASES' # 'PARTIAL_CONDENSING_HPIPM' # FULL_CONDENSING_QPOASES
    ocp.solver_options.hessian_approx = 'EXACT'
    ocp.solver_options.integrator_type = 'DISCRETE'
    ocp.solver_options.print_level = 0
    ocp.solver_options.tol = TOL
    ocp.solver_options.nlp_solver_type = 'SQP' # SQP_RTI, SQP

    ocp.solver_options.globalization = globalization
    ocp.solver_options.globalization_alpha_reduction = 0.9
    ocp.solver_options.globalization_line_search_use_sufficient_descent = line_search_use_sufficient_descent
    ocp.solver_options.globalization_use_SOC = globalization_use_SOC
    ocp.solver_options.globalization_eps_sufficient_descent = 5e-1

    SQP_max_iter = 200
    ocp.solver_options.qp_solver_iter_max = 400
    ocp.solver_options.nlp_solver_max_iter = SQP_max_iter

    # create solver
    ocp_solver = AcadosOcpSolver(ocp)

    # initialize solver
    xinit = np.array([1.0])
    [ocp_solver.set(i, "x", xinit) for i in range(N+1)]

    # get stats
    status = ocp_solver.solve()
    ocp_solver.print_statistics()
    iterations = ocp_solver.get_stats('sqp_iter')
    alphas = ocp_solver.get_stats('alpha')[1:]
    qp_iters = ocp_solver.get_stats('qp_iter')
    print(f"acados ocp solver returned status {status}")

    # get solution
    solution = ocp_solver.get(0, "x")
    print(f"found solution {solution}")

    # print summary
    print(f"cost function value = {ocp_solver.get_cost()} after {iterations} SQP iterations")
    print(f"alphas: {alphas[:iterations]}")
    print(f"total number of QP iterations: {sum(qp_iters[:iterations])}")

    # compare to analytical solution
    exact_solution = np.array([0.0])
    sol_err = max(np.abs(solution - exact_solution ))
    print(f"error wrt analytical solution {sol_err}")

    # checks
    if ocp.model.cost_expr_ext_cost_custom_hess_e == 1.0:
        if globalization == 'MERIT_BACKTRACKING':
            assert sol_err <= TOL * 1e1
            assert status == 0

            if line_search_use_sufficient_descent == 1:
                assert iterations <= 22
            else:
                assert iterations >= 64

        elif globalization == 'FIXED_STEP':
            assert status == 2

    print("\n\n----------------------\n")
