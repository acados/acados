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

import numpy as np
from acados_template import AcadosOcpBatchSolver, AcadosOcpSolver
from setup_parametric_ocp import PARAM_VALUE_DICT, export_parametric_ocp
import time


def main_sequential(x0, N_sim):

    learnable_params = ["A", "Q", "b"]
    ocp = export_parametric_ocp(PARAM_VALUE_DICT, learnable_params=learnable_params)
    ocp.solver_options.with_solution_sens_wrt_params = True

    solver = AcadosOcpSolver(ocp, verbose=False)

    nx = ocp.dims.nx
    nu = ocp.dims.nu

    simU = np.zeros((N_sim, nu))
    simX = np.zeros((N_sim+1, nx))
    simX[0,:] = x0

    adjoints = []
    param_vals = []

    t0 = time.time()
    for i in range(N_sim):
        solver.reset()

        # vary global parameters
        param = ocp.p_global_values.copy()
        param[i % param.shape[0]] += 0.05
        solver.set_p_global_and_precompute_dependencies(param)

        simU[i,:] = solver.solve_for_x0(x0_bar=simX[i, :])
        simX[i+1,:] = solver.get(1, "x")

        solver.setup_qp_matrices_and_factorize()
        sens_adj = solver.eval_adjoint_solution_sensitivity([(1, np.ones((ocp.dims.nx, 1)))], [(1, np.ones((ocp.dims.nu, 1)))])
        adjoints.append(sens_adj)
        param_vals.append(param)

    t_elapsed = 1e3 * (time.time() - t0)
    print("main_sequential, reset, set p_global, solve, adjoints and get:", f"{t_elapsed:.3f} ms\n")

    return simX, simU, param_vals, adjoints


def main_batch(Xinit, simU, param_vals, adjoints_ref, tol, num_threads_in_batch_solve=1):

    N_batch = Xinit.shape[0] - 1

    learnable_params = ["A", "Q", "b"]
    ocp = export_parametric_ocp(PARAM_VALUE_DICT, learnable_params=learnable_params)
    ocp.solver_options.with_solution_sens_wrt_params = True

    batch_solver = AcadosOcpBatchSolver(ocp, N_batch, num_threads_in_batch_solve=num_threads_in_batch_solve, verbose=False)

    # reset, set bounds and p_global
    t0 = time.time()
    for n in range(N_batch):
        batch_solver.ocp_solvers[n].constraints_set(0, "lbx", Xinit[n])
        batch_solver.ocp_solvers[n].constraints_set(0, "ubx", Xinit[n])
        batch_solver.ocp_solvers[n].reset()
        batch_solver.ocp_solvers[n].set_p_global_and_precompute_dependencies(param_vals[n])
    t_elapsed = 1e3 * (time.time() - t0)

    print(f"main_batch: with {num_threads_in_batch_solve} threads, reset, set x_0 and p_global: {t_elapsed:.3f} ms")

    # solve
    t0 = time.time()
    batch_solver.solve(N_batch)
    t_elapsed = 1e3 * (time.time() - t0)

    print(f"main_batch: with {num_threads_in_batch_solve} threads, solve: {t_elapsed:.3f} ms")

    for n in range(N_batch):
        u = batch_solver.ocp_solvers[n].get(0, "u")
        diff = np.linalg.norm(u-simU[n])
        if not diff < tol:
            raise Exception(f"solution should match sequential call up to {tol} got error {diff} for {n}th batch solve")

    # actually not needed for convex problem but we want to test it
    batch_solver.setup_qp_matrices_and_factorize(N_batch)

    # eval adjoint
    t0 = time.time()
    sens_adj = batch_solver.eval_adjoint_solution_sensitivity([(1, np.ones((N_batch, ocp.dims.nx, 1)))], [(1, np.ones((N_batch, ocp.dims.nu, 1)))])
    t_elapsed = 1e3 * (time.time() - t0)

    print(f"main_batch: with {num_threads_in_batch_solve} threads, adjoint solution sens: {t_elapsed:.3f} ms\n")

    for n in range(N_batch):
        diff = np.linalg.norm(sens_adj[n] - adjoints_ref[n])
        if not diff < tol:
            raise Exception(f"solution should match sequential call up to {tol} got error {diff} for {n}th batch solve")


if __name__ == "__main__":

    tol = 1e-12
    N_batch = 128
    x0 = np.array([0.1, -0.2])

    print("main sequential")
    simX, simU, param_vals, adjoints = main_sequential(x0=x0, N_sim=N_batch)

    print("main batch")
    main_batch(Xinit=simX, simU=simU, param_vals=param_vals, adjoints_ref=adjoints, tol=tol, num_threads_in_batch_solve=1)
    main_batch(Xinit=simX, simU=simU, param_vals=param_vals, adjoints_ref=adjoints, tol=tol, num_threads_in_batch_solve=4)
