# This test is an extension of the 'minimal_example_ocp_reuse_code.py' example.
#
# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# Licensed under the 2-Clause BSD License.
from acados_template import AcadosOcp, AcadosOcpSolver, ACADOS_INFTY
from pendulum_model import export_pendulum_ode_model
import numpy as np
import scipy.linalg
from utils import plot_pendulum
import pytest

PLOT = False

def run_one_sided_experiment(
        constraint_variant='one_sided',
        qp_solver='PARTIAL_CONDENSING_HPIPM',
        qp_solver_warm_start=0):

    # create ocp object to formulate the OCP
    ocp = AcadosOcp()

    # set model
    model = export_pendulum_ode_model()
    ocp.model = model

    nx = model.x.rows()
    nu = model.u.rows()
    ny = nx + nu
    ny_e = nx

    N_horizon = 20  # number of shooting nodes
    Tf = 1.0

    # set dimensions
    ocp.solver_options.N_horizon = N_horizon

    # set cost
    Q = 2 * np.diag([1e3, 1e3, 1e-2, 1e-2])
    R = 2 * np.diag([1e-2])

    ocp.cost.W_e = Q
    ocp.cost.W = scipy.linalg.block_diag(Q, R)

    ocp.cost.cost_type = 'LINEAR_LS'
    ocp.cost.cost_type_e = 'LINEAR_LS'

    ocp.cost.Vx = np.zeros((ny, nx))
    ocp.cost.Vx[:nx, :nx] = np.eye(nx)

    Vu = np.zeros((ny, nu))
    Vu[4, 0] = 1.0
    ocp.cost.Vu = Vu

    ocp.cost.Vx_e = np.eye(nx)

    ocp.cost.yref = np.zeros((ny,))
    ocp.cost.yref_e = np.zeros((ny_e,))

    # set constraints
    Fmax = 80
    ocp.constraints.lbu = np.array([-Fmax])
    ocp.constraints.ubu = np.array([+Fmax])
    ocp.constraints.idxbu = np.array([0])

    ocp.constraints.x0 = np.array([0.0, np.pi, 0.0, 0.0])

    if constraint_variant == 'one_sided':
        ocp.constraints.lbx = np.array([-ACADOS_INFTY])
        ocp.constraints.ubx = np.array([+5.0])
        ocp.constraints.idxbx = np.array([0])
        expected_status = 0
    elif constraint_variant == 'one_sided_wrong_infty':
        ocp.constraints.lbx = np.array([-0.5*ACADOS_INFTY])
        ocp.constraints.ubx = np.array([+5.0])
        ocp.constraints.idxbx = np.array([0])
        if qp_solver in ['FULL_CONDENSING_DAQP', 'PARTIAL_CONDENSING_OSQP', 'FULL_CONDENSING_QPOASES']:
            expected_status = 0
        elif qp_solver in ['FULL_CONDENSING_HPIPM', 'PARTIAL_CONDENSING_HPIPM']:
            # complementarity residual does not converge to tolerance if infty is not defined properly
            expected_status = 2
        elif qp_solver == 'PARTIAL_CONDENSING_CLARABEL':
            # QP solver fails if infty is not defined properly
            expected_status = 4

    # set options
    ocp.solver_options.qp_solver = qp_solver
    ocp.solver_options.qp_solver_iter_max = 4000
    ocp.solver_options.qp_solver_warm_start = qp_solver_warm_start
    # also warm start the first QP, i.e. the one right after the masks are changed
    ocp.solver_options.nlp_solver_warm_start_first_qp = qp_solver_warm_start > 0
    ocp.solver_options.hessian_approx = 'GAUSS_NEWTON'
    ocp.solver_options.integrator_type = 'ERK'
    ocp.solver_options.nlp_solver_type = 'SQP'
    ocp.solver_options.tol = 1e-7

    # set prediction horizon
    ocp.solver_options.tf = Tf

    # create solver
    ocp_solver = AcadosOcpSolver(ocp)

    # solve the problem defined here (original from code export), analog to 'minimal_example_ocp.py'
    simX0 = np.zeros((N_horizon + 1, nx))
    simU0 = np.zeros((N_horizon, nu))

    print(80*'-')
    print(f'solve original code with N = {N_horizon} and Tf = {Tf} s:')
    status = ocp_solver.solve()
    ocp_solver.print_statistics()

    assert status == expected_status

    # get solution
    iterate = ocp_solver.get_iterate()
    simX0[:] = iterate.x
    simU0[:] = iterate.u
    solutions = [ocp_solver.get_flat_iterate()]

    lambdas = [ocp_solver.get(i, "lam") for i in range(1, N_horizon)]
    for lam in lambdas:
        assert np.all(lam >= 0)

    # if unbounded constraint is defined properly, lambda should be zero
    i_infty = 1
    if constraint_variant == 'one_sided':
        for lam in lambdas:
            assert lam[i_infty] == 0
    elif (constraint_variant == 'one_sided_wrong_infty' and
        qp_solver in ['FULL_CONDENSING_HPIPM', 'PARTIAL_CONDENSING_HPIPM']):
        for lam in lambdas:
            assert lam[i_infty] != 0

    if constraint_variant == 'one_sided':
        stage = 1
        # check setting lambdas
        ocp_solver.set(stage, "lam", np.ones(lambdas[0].shape))
        lam = ocp_solver.get(stage, "lam")
        assert lam[i_infty] == 0

        # check updating bound
        i_infty_new = 2 # lam is ordered as lbu, lbx, ubu, ubx
        # the bound which is removed has to be active in the solution
        assert lambdas[0][i_infty_new] > 1e-3
        ocp_solver.constraints_set(stage, "lbx", -10.)
        ocp_solver.set(stage, "lam", np.ones(lambdas[0].shape))
        ocp_solver.constraints_set(stage, "ubu", ACADOS_INFTY)
        lam = ocp_solver.get(stage, "lam")

        assert lam[i_infty_new] == 0
        assert lam[i_infty] == 1.

        # check solving with updated bounds
        status = ocp_solver.solve()
        assert status == 0
        lam = ocp_solver.get(stage, "lam")
        assert lam[i_infty_new] == 0
        assert lam[i_infty] >= 0
        solutions.append(ocp_solver.get_flat_iterate())

        # check solving with original bounds again
        ocp_solver.constraints_set(stage, "lbx", -ACADOS_INFTY)
        ocp_solver.constraints_set(stage, "ubu", Fmax)
        status = ocp_solver.solve()
        assert status == 0
        lam = ocp_solver.get(stage, "lam")
        assert lam[i_infty] == 0
        assert lam[i_infty_new] >= 0
        solutions.append(ocp_solver.get_flat_iterate())

    if PLOT:
        plot_pendulum(np.linspace(0, Tf, N_horizon + 1), Fmax, simU0, simX0, latexify=False, plt_show=True, X_true_label=f'N={N_horizon}, Tf={Tf}')

    ocp_solver = None
    return solutions


@pytest.fixture(scope="module")
def reference_solutions():
    return run_one_sided_experiment()


QP_SOLVERS = [
    'FULL_CONDENSING_HPIPM',
    'PARTIAL_CONDENSING_HPIPM',
    'FULL_CONDENSING_DAQP',
    'PARTIAL_CONDENSING_OSQP',
    'PARTIAL_CONDENSING_CLARABEL',
    'FULL_CONDENSING_QPOASES',
]


@pytest.mark.parametrize("qp_solver", QP_SOLVERS)
@pytest.mark.parametrize("qp_solver_warm_start", [0, 1])
@pytest.mark.parametrize(
    "constraint_variant", ["one_sided", "one_sided_wrong_infty"])
def test_one_sided_constraints(
        constraint_variant, qp_solver, qp_solver_warm_start, reference_solutions):
    solutions = run_one_sided_experiment(
        constraint_variant=constraint_variant,
        qp_solver=qp_solver,
        qp_solver_warm_start=qp_solver_warm_start,
    )
    if constraint_variant == 'one_sided':
        for solution, reference in zip(solutions, reference_solutions):
            diff_primal, _ = (solution - reference).primal_dual_inf_norm()
            assert diff_primal < 1e-6
