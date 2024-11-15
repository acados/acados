# -*- coding: future_fstrings -*-
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

import sys
sys.path.insert(0, '../pendulum_on_cart/solution_sensitivities')
import numpy as np
from sensitivity_utils import plot_cost_gradient_results
from acados_template import AcadosOcpBatchSolver, AcadosOcpSolver
from setup_parametric_ocp import PARAM_VALUE_DICT, export_parametric_ocp
import time


def setup_example(learnable_params, num_threads_in_batch_solve: int = 1):
    ocp = export_parametric_ocp(PARAM_VALUE_DICT, learnable_params=learnable_params, num_threads_in_batch_solve = num_threads_in_batch_solve)
    ocp.solver_options.with_solution_sens_wrt_params = True
    ocp.solver_options.nlp_solver_type = "SQP"
    ocp.solver_options.hessian_approx = "EXACT"
    return ocp


def main(N_batch: int = 64):

    learnable_params = ["A", "Q", "b"]
    x0 = np.array([0.1, -0.2])

    delta_p = 0.005
    p_nominal = np.concatenate([PARAM_VALUE_DICT[k].flatten() for k in learnable_params]).flatten()

    p_test = np.arange(0.0, 1.0, delta_p)
    np_test = p_test.shape[0]

    dp = np.ones(p_nominal.shape).reshape((-1,))
    # Q
    dp[4] = 10.0
    dp[5] = 0.0
    dp[6] = 0.0
    dp[7] = 10.0
    # b
    dp[8] = -0.2
    dp[9] = -0.2

    ocp = export_parametric_ocp(PARAM_VALUE_DICT, learnable_params = learnable_params)
    ocp.solver_options.with_value_sens_wrt_params = True
    ocp.solver_options.nlp_solver_type = "SQP"
    acados_ocp_solver = AcadosOcpBatchSolver(ocp, N_batch)

    optimal_value_grad = np.zeros((np_test,))
    optimal_value = np.zeros((np_test,))

    pi = np.zeros(np_test)
    for i, p in enumerate(p_test):
        p_val = p_nominal + p * dp
        acados_ocp_solver.set_p_global_and_precompute_dependencies(p_val)
        pi[i] = acados_ocp_solver.solve_for_x0(x0, fail_on_nonzero_status=False)[0]

        status = acados_ocp_solver.get_status()
        if status != 0:
            print(f"Solver failed with status {status} for p_val = {p_val}.")

        optimal_value[i] = acados_ocp_solver.get_cost()
        optimal_value_grad[i] = acados_ocp_solver.eval_and_get_optimal_value_gradient("p_global") @ dp

    # evaluate cost gradient
    optimal_value_grad_via_fd = np.gradient(optimal_value, delta_p)
    cost_reconstructed_np_grad = np.cumsum(optimal_value_grad_via_fd) * delta_p + optimal_value[0]
    cost_reconstructed_acados = np.cumsum(optimal_value_grad) * delta_p + optimal_value[0]

    plot_cost_gradient_results(p_test, optimal_value, optimal_value_grad,
                               optimal_value_grad_via_fd, cost_reconstructed_np_grad,
                               cost_reconstructed_acados, y_scale_log=True,
                               title=f"varying parameters {', '.join(learnable_params)} in direction {dp}",
                               xlabel=r"$\alpha$ in $p+\alpha \Delta p$")

    # checks
    test_tol = 1e-3
    median_diff = np.median(np.abs(optimal_value_grad - optimal_value_grad_via_fd))
    print(f"Median difference between value function gradient obtained by acados and via FD is {median_diff:.2e} should be < {test_tol:.2e}.")
    assert median_diff <= test_tol



def main_sequential(x0, N_sim):

    learnable_params = ["A", "Q", "b"]
    ocp = setup_example(learnable_params)
    solver = AcadosOcpSolver(ocp, verbose=False)

    nx = ocp.dims.nx
    nu = ocp.dims.nu

    simU = np.zeros((N_sim, nu))
    simX = np.zeros((N_sim+1, nx))
    simX[0,:] = x0

    adjoints = []
    x_trajs = []

    t0 = time.time()
    for i in range(N_sim):
        simU[i,:] = solver.solve_for_x0(x0_bar=simX[i, :])
        simX[i+1,:] = solver.get(1, "x")
        sens_adj = solver.eval_adjoint_solution_sensitivity([(1, np.ones((ocp.dims.nx, 1)))], [(1, np.ones((ocp.dims.nu, 1)))])
        adjoints.append(sens_adj)
        x_trajs.append(solver.get_flat('x'))

    t_elapsed = 1e3 * (time.time() - t0)
    print("main_sequential, solve, adjoints and get:", f"{t_elapsed:.3f}ms")

    return simX, simU, adjoints, x_trajs


def main_batch(Xinit, simU, x_trajs, adjoints_ref, tol, num_threads_in_batch_solve=1):

    N_batch = Xinit.shape[0] - 1

    learnable_params = ["A", "Q", "b"]
    ocp = setup_example(learnable_params, num_threads_in_batch_solve)
    batch_solver = AcadosOcpBatchSolver(ocp, N_batch, verbose=False)

    for n in range(N_batch):
        batch_solver.ocp_solvers[n].constraints_set(0, "lbx", Xinit[n])
        batch_solver.ocp_solvers[n].constraints_set(0, "ubx", Xinit[n])

    # set initial guess
    Xinit_batch = np.array([np.tile(Xinit[i], (ocp.solver_options.N_horizon+1,)) for i in range(N_batch)])
    t0 = time.time()
    batch_solver.set_flat('x', Xinit_batch)
    t_elapsed = 1e3 * (time.time() - t0)

    print(f"main_batch: with {num_threads_in_batch_solve} threads, set_flat: {t_elapsed:.3f}ms")

    # solve
    t0 = time.time()
    batch_solver.solve()
    t_elapsed = 1e3 * (time.time() - t0)

    print(f"main_batch: with {num_threads_in_batch_solve} threads, solve: {t_elapsed:.3f}ms")

    for n in range(N_batch):
        u = batch_solver.ocp_solvers[n].get(0, "u")

        if not np.linalg.norm(u-simU[n]) < tol:
            raise Exception(f"solution should match sequential call up to {tol} got error {np.linalg.norm(u-simU[n])} for {n}th batch solve")

        x_traj = batch_solver.ocp_solvers[n].get_flat('x')
        if not np.linalg.norm(x_traj-x_trajs[n]) < tol:
            raise Exception(f"solution should match sequential call up to {tol} got error {np.linalg.norm(u-simU[n])} for {n}th batch solve")

    # eval adjoint
    t0 = time.time()
    sens_adj = batch_solver.eval_adjoint_solution_sensitivity([(1, np.ones((N_batch, ocp.dims.nx, 1)))], [(1, np.ones((N_batch, ocp.dims.nu, 1)))])
    t_elapsed = 1e3 * (time.time() - t0)

    print(f"main_batch: with {num_threads_in_batch_solve} threads, adjoint solution sens: {t_elapsed:.3f}ms")

    for n in range(N_batch):

        print("---")
        # print(sens_adj[n] - adjoints_ref[n])
        print(np.max(np.abs(sens_adj[n] - adjoints_ref[n])))
        # if not np.max(np.abs(sens_adj[n] - adjoints_ref[n])) < tol*10:
        #     raise Exception(f"solution should match sequential call up to {tol*10} got error {np.linalg.norm(sens_adj[n] - adjoints_ref[n])} for {n}th batch solve")


if __name__ == "__main__":

    tol = 1e-5
    # N_batch = 256
    N_batch = 16
    x0 = np.array([0.1, -0.2])

    simX, simU, adjoints, x_trajs = main_sequential(x0=x0, N_sim=N_batch)

    print("BATCH NUM THREADS = 1")
    main_batch(Xinit=simX, simU=simU, x_trajs=x_trajs, adjoints_ref=adjoints, tol=tol, num_threads_in_batch_solve=1)

    print("BATCH NUM THREADS = 4")
    main_batch(Xinit=simX, simU=simU, x_trajs=x_trajs, adjoints_ref=adjoints, tol=tol, num_threads_in_batch_solve=4)
