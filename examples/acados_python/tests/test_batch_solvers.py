#
# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# Licensed under the 2-Clause BSD License.

"""
Tests that verify AcadosOcpBatchSolver and AcadosSimBatchSolver produce the same
results as calling the corresponding non-batch solvers individually.

The individual comparison is done via AcadosOcpBatchSolver.ocp_solvers[i] and
AcadosSimBatchSolver.sim_solvers[i] to avoid creating additional solver instances.
"""

from test_paths import add_example_paths

add_example_paths()



import numpy as np
import scipy.linalg
import casadi as ca
import pytest

from acados_template import (
    AcadosOcp, AcadosOcpBatchSolver,
    AcadosSim, AcadosSimBatchSolver,
)
from pendulum_model import export_pendulum_ode_model
from setup_parametric_ocp import export_parametric_ocp

TOL = 1e-7
N_BATCH = 5


# ---------------------------------------------------------------------------
# OCP helpers
# ---------------------------------------------------------------------------

def setup_ocp(tol: float = TOL) -> AcadosOcp:
    ocp = AcadosOcp()
    ocp.model = export_pendulum_ode_model()

    Tf = 1.0
    N = 20
    ocp.solver_options.N_horizon = N

    Q_mat = 2 * np.diag([1e3, 1e3, 1e-2, 1e-2])
    R_mat = 2 * np.diag([1e-2])
    cost_W = scipy.linalg.block_diag(Q_mat, R_mat)

    ocp.cost.cost_type = 'NONLINEAR_LS'
    ocp.cost.cost_type_e = 'NONLINEAR_LS'
    ocp.cost.W_e = Q_mat
    ocp.cost.W = cost_W
    ocp.model.cost_y_expr = ca.vertcat(ocp.model.x, ocp.model.u)
    ocp.model.cost_y_expr_e = ocp.model.x
    ocp.cost.yref = np.zeros(ocp.model.cost_y_expr.shape).flatten()
    ocp.cost.yref_e = np.zeros(ocp.model.cost_y_expr_e.shape).flatten()

    Fmax = 80
    ocp.constraints.lbu = np.array([-Fmax])
    ocp.constraints.ubu = np.array([+Fmax])
    ocp.constraints.idxbu = np.array([0])
    ocp.constraints.x0 = np.array([0.0, np.pi, 0.0, 0.0])

    ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM'
    ocp.solver_options.hessian_approx = 'GAUSS_NEWTON'
    ocp.solver_options.integrator_type = 'IRK'
    ocp.solver_options.nlp_solver_type = 'SQP'
    ocp.solver_options.nlp_solver_tol_stat = tol
    ocp.solver_options.nlp_solver_tol_eq = tol
    ocp.solver_options.nlp_solver_tol_ineq = tol
    ocp.solver_options.nlp_solver_tol_comp = tol
    ocp.solver_options.tf = Tf
    ocp.solver_options.with_batch_functionality = True

    return ocp


# ---------------------------------------------------------------------------
# Sim helpers
# ---------------------------------------------------------------------------

def setup_sim() -> AcadosSim:
    sim = AcadosSim()
    sim.model = export_pendulum_ode_model()
    sim.solver_options.T = 0.2
    sim.solver_options.integrator_type = 'IRK'
    sim.solver_options.num_stages = 3
    sim.solver_options.num_steps = 3
    sim.solver_options.newton_iter = 10
    sim.solver_options.with_batch_functionality = True
    return sim


# ---------------------------------------------------------------------------
# Test functions
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def ocp_batch_case():
    x0_batch = np.tile(np.array([0.0, np.pi, 0.0, 0.0]), (N_BATCH, 1))
    rng = np.random.default_rng(0)
    x0_batch += 0.1 * rng.standard_normal(x0_batch.shape)

    ocp = setup_ocp()
    batch_solver = AcadosOcpBatchSolver(ocp, N_batch_init=N_BATCH, num_threads_in_batch_solve=1, verbose=False)

    # --- num_threads property ---
    assert batch_solver.num_threads_in_batch_solve == 1
    batch_solver.num_threads_in_batch_solve = 2
    assert batch_solver.num_threads_in_batch_solve == 2
    batch_solver.num_threads_in_batch_solve = 1

    batch_solver.constraints_set(0, 'lbx', x0_batch)
    batch_solver.constraints_set(0, 'ubx', x0_batch)
    batch_solver.solve()
    return ocp, batch_solver


@pytest.mark.parametrize("i", range(N_BATCH))
def test_ocp_batch_get_and_get_flat(ocp_batch_case, i):
    _, batch_solver = ocp_batch_case
    batch_x1 = batch_solver.get(1, 'x')
    batch_u0 = batch_solver.get(0, 'u')
    err_x = np.max(np.abs(batch_solver.ocp_solvers[i].get(1, 'x') - batch_x1[i]))
    err_u = np.max(np.abs(batch_solver.ocp_solvers[i].get(0, 'u') - batch_u0[i]))
    assert err_x < TOL * 100, f"get(1,'x') mismatch at index {i}: {err_x:.2e}"
    assert err_u < TOL * 100, f"get(0,'u') mismatch at index {i}: {err_u:.2e}"

    batch_xflat = batch_solver.get_flat('x')
    batch_uflat = batch_solver.get_flat('u')
    err_x = np.max(np.abs(batch_solver.ocp_solvers[i].get_flat('x') - batch_xflat[i]))
    err_u = np.max(np.abs(batch_solver.ocp_solvers[i].get_flat('u') - batch_uflat[i]))
    assert err_x < TOL * 100, f"get_flat('x') mismatch at index {i}: {err_x:.2e}"
    assert err_u < TOL * 100, f"get_flat('u') mismatch at index {i}: {err_u:.2e}"


def test_ocp_batch_rejects_oversized_request(ocp_batch_case):
    _, batch_solver = ocp_batch_case
    with pytest.raises(ValueError):
        batch_solver.get_flat('x', N_BATCH + 1)


@pytest.fixture(scope="module")
def ocp_initial_sensitivity_case():
    x0_batch = np.tile(np.array([0.0, np.pi, 0.0, 0.0]), (N_BATCH, 1))
    rng = np.random.default_rng(2)
    x0_batch += 0.1 * rng.standard_normal(x0_batch.shape)

    ocp = setup_ocp()
    batch_solver = AcadosOcpBatchSolver(ocp, N_batch_init=N_BATCH, verbose=False)
    batch_solver.constraints_set(0, 'lbx', x0_batch)
    batch_solver.constraints_set(0, 'ubx', x0_batch)
    batch_solver.solve()
    batch_solver.setup_qp_matrices_and_factorize()

    batch_d = batch_solver.eval_solution_sensitivity(
        0,
        sanity_checks=False,
        with_respect_to='initial_state',
        return_sens_x=True,
        return_sens_u=True,
    )
    return batch_solver, batch_d


@pytest.mark.parametrize("i", range(N_BATCH))
def test_ocp_batch_eval_solution_sensitivity_wrt_initial_state(
    ocp_initial_sensitivity_case, i
):
    batch_solver, batch_d = ocp_initial_sensitivity_case
    solver_i = batch_solver.ocp_solvers[i]
    solver_i.setup_qp_matrices_and_factorize()
    d_i = solver_i.eval_solution_sensitivity(
        0,
        sanity_checks=False,
        with_respect_to='initial_state',
        return_sens_x=True,
        return_sens_u=True,
    )
    err_u = np.max(np.abs(d_i['sens_u'] - batch_d['sens_u'][i]))
    err_x = np.max(np.abs(d_i['sens_x'] - batch_d['sens_x'][i]))
    assert err_u < TOL * 100, f"sens_u mismatch at index {i}: {err_u:.2e}"
    assert err_x < TOL * 100, f"sens_x mismatch at index {i}: {err_x:.2e}"


@pytest.fixture(scope="module")
def ocp_parametric_sensitivity_case():
    ocp_nx, ocp_nu = 4, 2
    learnable_params = ["A", "Q", "b"]

    x0_batch = 0.1 * (-1) ** np.arange(ocp_nx)
    x0_batch = np.tile(x0_batch, (N_BATCH, 1))
    rng = np.random.default_rng(5)
    x0_batch += 0.05 * rng.standard_normal(x0_batch.shape)

    seed_x_val = np.ones((ocp_nx, 1))
    seed_u_val = np.ones((ocp_nu, 1))

    # Build OCP with both forward and adjoint sensitivity support
    ocp = export_parametric_ocp(ocp_nx, ocp_nu, learnable_params=learnable_params)
    ocp.code_gen_options.with_solution_sens_wrt_params_forw = True
    ocp.code_gen_options.with_solution_sens_wrt_params_adj = True
    ocp.solver_options.qp_solver_ric_alg = 0
    ocp.solver_options.with_batch_functionality = True

    batch_solver = AcadosOcpBatchSolver(ocp, N_batch_init=N_BATCH, verbose=False)
    batch_solver.constraints_set(0, 'lbx', x0_batch)
    batch_solver.constraints_set(0, 'ubx', x0_batch)
    for n in range(N_BATCH):
        batch_solver.ocp_solvers[n].reset()
    batch_solver.solve()
    batch_solver.setup_qp_matrices_and_factorize()

    n_p_global = ocp.p_global_values.shape[0]

    # --- forward sensitivity wrt p_global ---
    batch_fwd = batch_solver.eval_solution_sensitivity(
        0,
        with_respect_to='p_global',
        return_sens_x=True,
        return_sens_u=True,
    )

    # --- adjoint sensitivity wrt p_global ---
    batch_seed_x = np.tile(seed_x_val[np.newaxis, :, :], (N_BATCH, 1, 1))
    batch_seed_u = np.tile(seed_u_val[np.newaxis, :, :], (N_BATCH, 1, 1))
    batch_adj = batch_solver.eval_adjoint_solution_sensitivity(
        seed_x=[(1, batch_seed_x)],
        seed_u=[(0, batch_seed_u)],
    )
    return batch_solver, batch_fwd, batch_adj, seed_x_val, seed_u_val


@pytest.mark.parametrize("i", range(N_BATCH))
def test_ocp_batch_parametric_sensitivity(ocp_parametric_sensitivity_case, i):
    batch_solver, batch_fwd, batch_adj, seed_x_val, seed_u_val = (
        ocp_parametric_sensitivity_case
    )
    solver_i = batch_solver.ocp_solvers[i]
    solver_i.setup_qp_matrices_and_factorize()

    d_i = solver_i.eval_solution_sensitivity(
        0,
        with_respect_to='p_global',
        return_sens_x=True,
        return_sens_u=True,
    )
    err_u = np.max(np.abs(d_i['sens_u'] - batch_fwd['sens_u'][i]))
    err_x = np.max(np.abs(d_i['sens_x'] - batch_fwd['sens_x'][i]))
    assert err_u < TOL * 100, f"fwd sens_u mismatch at index {i}: {err_u:.2e}"
    assert err_x < TOL * 100, f"fwd sens_x mismatch at index {i}: {err_x:.2e}"

    solver_i.setup_qp_matrices_and_factorize()
    adj_i = solver_i.eval_adjoint_solution_sensitivity(
        seed_x=[(1, seed_x_val)],
        seed_u=[(0, seed_u_val)],
    )
    err_adj = np.max(np.abs(adj_i - batch_adj[i]))
    assert err_adj < TOL * 100, f"adj sens mismatch at index {i}: {err_adj:.2e}"


@pytest.fixture(scope="module")
def sim_batch_case():
    rng = np.random.default_rng(4)
    x_batch = rng.standard_normal((N_BATCH, 4))
    u_batch = rng.standard_normal((N_BATCH, 1))

    sim = setup_sim()
    batch_integrator = AcadosSimBatchSolver(sim, N_batch_init=N_BATCH, num_threads_in_batch_solve=1, verbose=False)

    # --- num_threads property ---
    assert batch_integrator.num_threads_in_batch_solve == 1
    batch_integrator.num_threads_in_batch_solve = 2
    assert batch_integrator.num_threads_in_batch_solve == 2
    batch_integrator.num_threads_in_batch_solve = 1

    # Set inputs and run batch solve
    for n in range(N_BATCH):
        batch_integrator.sim_solvers[n].set('x', x_batch[n])
        batch_integrator.sim_solvers[n].set('u', u_batch[n])
    batch_integrator.solve()

    # Collect batch outputs
    x_next_batch = np.array([batch_integrator.sim_solvers[n].get('x') for n in range(N_BATCH)])
    return batch_integrator, x_batch, u_batch, x_next_batch


@pytest.mark.parametrize("i", range(N_BATCH))
def test_sim_batch(sim_batch_case, i):
    batch_integrator, x_batch, u_batch, x_next_batch = sim_batch_case
    x_next_ref = batch_integrator.sim_solvers[i].simulate(x=x_batch[i], u=u_batch[i])
    err = np.max(np.abs(x_next_ref - x_next_batch[i]))
    assert err < 1e-10, f"SimBatchSolver: error at index {i} {err:.2e} exceeds tolerance"


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__]))

