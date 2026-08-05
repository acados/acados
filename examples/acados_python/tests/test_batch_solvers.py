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

"""
Tests that verify AcadosOcpBatchSolver and AcadosSimBatchSolver produce the same
results as calling the corresponding non-batch solvers individually.

The individual comparison is done via AcadosOcpBatchSolver.ocp_solvers[i] and
AcadosSimBatchSolver.sim_solvers[i] to avoid creating additional solver instances.
"""

import sys
sys.path.insert(0, '../pendulum_on_cart/common')
sys.path.insert(0, '../solution_sensitivities_convex_example')

import numpy as np
import scipy.linalg
import casadi as ca

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

def test_ocp_batch_get_and_get_flat():
    """
    Create one AcadosOcpBatchSolver, solve N_BATCH problems, then verify that:
      - batch.get()      matches ocp_solvers[i].get()
      - batch.get_flat() matches ocp_solvers[i].get_flat()
      - num_threads_in_batch_solve property getter/setter works
      - requesting n_batch > N_batch_max raises ValueError
    """
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

    nx = ocp.dims.nx
    nu = ocp.dims.nu

    # --- get() ---
    batch_x1 = batch_solver.get(1, 'x')   # shape (N_BATCH, nx)
    batch_u0 = batch_solver.get(0, 'u')   # shape (N_BATCH, nu)
    for i in range(N_BATCH):
        err_x = np.max(np.abs(batch_solver.ocp_solvers[i].get(1, 'x') - batch_x1[i]))
        err_u = np.max(np.abs(batch_solver.ocp_solvers[i].get(0, 'u') - batch_u0[i]))
        assert err_x < TOL * 100, f"get(1,'x') mismatch at index {i}: {err_x:.2e}"
        assert err_u < TOL * 100, f"get(0,'u') mismatch at index {i}: {err_u:.2e}"

    # --- get_flat() ---
    batch_xflat = batch_solver.get_flat('x')   # shape (N_BATCH, nx*(N+1))
    batch_uflat = batch_solver.get_flat('u')   # shape (N_BATCH, nu*N)
    for i in range(N_BATCH):
        err_x = np.max(np.abs(batch_solver.ocp_solvers[i].get_flat('x') - batch_xflat[i]))
        err_u = np.max(np.abs(batch_solver.ocp_solvers[i].get_flat('u') - batch_uflat[i]))
        assert err_x < TOL * 100, f"get_flat('x') mismatch at index {i}: {err_x:.2e}"
        assert err_u < TOL * 100, f"get_flat('u') mismatch at index {i}: {err_u:.2e}"

    # --- n_batch > N_batch_max should raise ---
    try:
        batch_solver.get_flat('x', N_BATCH + 1)
        raise AssertionError("Expected ValueError was not raised")
    except ValueError:
        pass

    print("test_ocp_batch_get_and_get_flat: PASSED")


def test_ocp_batch_eval_solution_sensitivity_wrt_initial_state():
    """
    Create one AcadosOcpBatchSolver, solve N_BATCH problems, then verify that
    batch.eval_solution_sensitivity(with_respect_to='initial_state') agrees with
    ocp_solvers[i].eval_solution_sensitivity() called individually.
    """
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

    nx = ocp.dims.nx
    nu = ocp.dims.nu

    for i in range(N_BATCH):
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

    print("test_ocp_batch_eval_solution_sensitivity_wrt_initial_state: PASSED")


def test_ocp_batch_parametric_sensitivity():
    """
    Create one AcadosOcpBatchSolver with the parametric convex OCP, then verify:
      - eval_solution_sensitivity(with_respect_to='p_global') agrees with ocp_solvers[i]
      - eval_adjoint_solution_sensitivity() agrees with ocp_solvers[i]
    Both forward and adjoint sensitivities are compiled into the same solver
    (with_solution_sens_wrt_params_forw=True, with_solution_sens_wrt_params_adj=True).
    """
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

    # compare against individual ocp_solvers[i]
    for i in range(N_BATCH):
        solver_i = batch_solver.ocp_solvers[i]
        solver_i.setup_qp_matrices_and_factorize()

        # forward
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

        # adjoint
        solver_i.setup_qp_matrices_and_factorize()
        adj_i = solver_i.eval_adjoint_solution_sensitivity(
            seed_x=[(1, seed_x_val)],
            seed_u=[(0, seed_u_val)],
        )
        err_adj = np.max(np.abs(adj_i - batch_adj[i]))
        assert err_adj < TOL * 100, f"adj sens mismatch at index {i}: {err_adj:.2e}"

    print("test_ocp_batch_parametric_sensitivity: PASSED")


def test_sim_batch():
    """
    Create one AcadosSimBatchSolver, run N_BATCH integrations, then verify:
      - sim_solvers[i].get('x') after batch solve matches a fresh individual simulate() call
      - num_threads_in_batch_solve property getter/setter works
    The individual reference is obtained by calling sim_solvers[i].simulate() on each
    solver from the same batch object, so no additional solver is created.
    """
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

    # Individual reference: run each sim_solver sequentially via simulate()
    x_next_ref = np.zeros_like(x_next_batch)
    for n in range(N_BATCH):
        x_next_ref[n] = batch_integrator.sim_solvers[n].simulate(x=x_batch[n], u=u_batch[n])

    err = np.max(np.abs(x_next_ref - x_next_batch))
    assert err < 1e-10, f"SimBatchSolver: max error {err:.2e} exceeds tolerance"

    print("test_sim_batch: PASSED")


if __name__ == '__main__':
    test_ocp_batch_get_and_get_flat()
    test_ocp_batch_eval_solution_sensitivity_wrt_initial_state()
    test_ocp_batch_parametric_sensitivity()
    test_sim_batch()
    print("\nAll batch solver tests passed.")

