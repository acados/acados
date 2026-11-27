#

from test_paths import add_example_paths

add_example_paths()
# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# Licensed under the 2-Clause BSD License.

from acados_template import AcadosOcp, AcadosOcpSolver, AcadosSimSolver, get_quadratic_penalty_expression
from pendulum_model import export_pendulum_ode_model
from utils import plot_pendulum
import numpy as np
import scipy.linalg
import pytest

SOLVER_TOL = 1e-7
N = 20
TOL_CLOSED_LOOP_TEST = 5e-5
PLOT_RESULT = False

QP_SOLVERS = ('PARTIAL_CONDENSING_HPIPM', \
                'FULL_CONDENSING_QPOASES', 'FULL_CONDENSING_HPIPM', \
                'FULL_CONDENSING_DAQP', 'PARTIAL_CONDENSING_OSQP', 'PARTIAL_CONDENSING_DAQP'
                )
# no soft constraint support:
# 'PARTIAL_CONDENSING_QPDUNES', 'PARTIAL_CONDENSING_OSQP', \


def run_closed_loop_experiment(soft_constr_type='bx', verbose=False, qp_solver='PARTIAL_CONDENSING_HPIPM', quadratic_penalty_only: bool = False):
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

    # set dimensions
    # NOTE: all dimensions but N ar detected
    ocp.solver_options.N_horizon = N

    # set cost module
    ocp.cost.cost_type = 'LINEAR_LS'
    ocp.cost.cost_type_e = 'LINEAR_LS'

    Q_mat = 2*np.diag([1e3, 1e3, 1e-2, 1e-2])
    R_mat = 2*np.diag([1e-2])

    ocp.cost.W = scipy.linalg.block_diag(Q_mat, R_mat)

    ocp.cost.Vx = np.zeros((ny, nx))
    ocp.cost.Vx[:nx,:nx] = np.eye(nx)

    Vu = np.zeros((ny, nu))
    Vu[4,0] = 1.0
    ocp.cost.Vu = Vu

    ocp.cost.Vx_e = np.eye(nx)
    ocp.cost.W_e = Q_mat

    ocp.cost.yref = np.zeros((ny, ))
    ocp.cost.yref_e = np.zeros((ny_e, ))

    # set constraints
    Fmax = 80
    vmax = 5

    x0 = np.array([0.0, np.pi, 0.0, 0.0])
    ocp.constraints.x0 = x0

    # bound on u
    ocp.constraints.lbu = np.array([-Fmax])
    ocp.constraints.ubu = np.array([+Fmax])
    ocp.constraints.idxbu = np.array([0])

    ocp.cost.zl = (0 if quadratic_penalty_only else 50)*np.ones((1,))
    ocp.cost.Zl = 10*np.ones((1,))
    ocp.cost.zu = (0 if quadratic_penalty_only else 50)*np.ones((1,))
    ocp.cost.Zu = 10*np.ones((1,))

    # soft constraint on x (either via bx or h)
    if soft_constr_type == 'bx':
        # soft bound on x
        ocp.constraints.lbx = np.array([-vmax])
        ocp.constraints.ubx = np.array([+vmax])
        ocp.constraints.idxbx = np.array([2]) # v is x[2]
        # indices of slacked constraints within bx
        ocp.constraints.idxsbx = np.array([0])

    elif soft_constr_type == 'h':
        # soft bound on x, using constraint h
        v1 = ocp.model.x[2]
        ocp.model.con_h_expr = v1

        ocp.constraints.lh = np.array([-vmax])
        ocp.constraints.uh = np.array([+vmax])
        # indices of slacked constraints within h
        ocp.constraints.idxsh = np.array([0])
    elif soft_constr_type == 'h_penalty':
        v1 = ocp.model.x[2]
        penalty_term = get_quadratic_penalty_expression(v1, -vmax, vmax, ocp.cost.Zl, ocp.cost.Zu)
        ocp.translate_cost_to_external_cost()
        ocp.model.cost_expr_ext_cost += penalty_term

        ocp.cost.zl = []
        ocp.cost.zu = []
        ocp.cost.Zu = []
        ocp.cost.Zl = []
    else:
        raise Exception(f"soft_constr_type must be 'bx', or 'h', got {soft_constr_type}.")

    # set options
    ocp.solver_options.qp_solver = qp_solver
    ocp.solver_options.hessian_approx = 'GAUSS_NEWTON'
    ocp.solver_options.integrator_type = 'ERK'
    ocp.solver_options.tf = Tf
    ocp.solver_options.nlp_solver_type = 'SQP'
    ocp.solver_options.tol = SOLVER_TOL
    ocp.solver_options.nlp_solver_ext_qp_res = 1

    ocp.solver_options.qp_solver_warm_start = 0
    ocp.solver_options.qp_solver_iter_max = 10000

    acados_ocp_solver = AcadosOcpSolver(ocp)
    acados_integrator = AcadosSimSolver(ocp)

    # closed loop
    Nsim = 20
    simX = np.zeros((Nsim+1, nx))
    simU = np.zeros((Nsim, nu))
    xcurrent = x0
    qp_iter = np.zeros(Nsim)
    sqp_iter = np.zeros(Nsim)

    for i in range(Nsim):
        simX[i,:] = xcurrent

        # solve ocp
        acados_ocp_solver.set(0, "lbx", xcurrent)
        acados_ocp_solver.set(0, "ubx", xcurrent)

        status = acados_ocp_solver.solve()
        if verbose:
            acados_ocp_solver.print_statistics()

        assert status == 0

        qp_iter[i] = np.sum(acados_ocp_solver.get_stats('qp_iter'))
        sqp_iter[i] = acados_ocp_solver.get_stats('sqp_iter')

        simU[i,:] = acados_ocp_solver.get(0, "u")

        # simulate system
        acados_integrator.set("x", xcurrent)
        acados_integrator.set("u", simU[i,:])

        status = acados_integrator.solve()
        assert status == 0

        # update state
        xcurrent = acados_integrator.get("x")

    simX[Nsim,:] = xcurrent

    # get slack values at stage 1
    sl = acados_ocp_solver.get(1, "sl")
    su = acados_ocp_solver.get(1, "su")
    print("sl", sl, "su", su)

    if PLOT_RESULT:
        dt = Tf / N
        plot_pendulum(np.linspace(0, Nsim*dt, Nsim+1), Fmax, simU, simX, latexify=False)

    print(f"\nsoft constraint example: ran formulation {soft_constr_type} with {qp_solver} successfully.\n")
    print(f"took {sqp_iter} SQP iterations and {qp_iter} QP iterations.")

    del acados_ocp_solver
    return simX, simU, sqp_iter


@pytest.fixture(scope="module")
def soft_constraint_reference():
    return run_closed_loop_experiment(
        soft_constr_type="bx",
        qp_solver="PARTIAL_CONDENSING_HPIPM",
    )


@pytest.mark.parametrize("soft_constr_type", ["bx", "h"])
@pytest.mark.parametrize("qp_solver", QP_SOLVERS)
def test_soft_constraint_formulations(
        soft_constr_type, qp_solver, soft_constraint_reference):
    result = run_closed_loop_experiment(
        soft_constr_type=soft_constr_type,
        qp_solver=qp_solver,
    )
    assert max(np.linalg.norm(value - reference)
               for value, reference in zip(result[:2], soft_constraint_reference[:2])) \
        <= TOL_CLOSED_LOOP_TEST
    np.testing.assert_array_equal(result[2], soft_constraint_reference[2])


@pytest.fixture(scope="module")
def quadratic_penalty_reference():
    return run_closed_loop_experiment(
        soft_constr_type="h",
        quadratic_penalty_only=True,
    )


@pytest.mark.parametrize("soft_constr_type", ["h", "h_penalty"])
def test_quadratic_penalty_formulations(
        soft_constr_type, quadratic_penalty_reference):
    result = run_closed_loop_experiment(
        soft_constr_type=soft_constr_type,
        quadratic_penalty_only=True,
    )
    assert max(np.linalg.norm(value - reference)
               for value, reference in zip(result[:2], quadratic_penalty_reference[:2])) \
        <= TOL_CLOSED_LOOP_TEST
    np.testing.assert_array_equal(result[2], quadratic_penalty_reference[2])