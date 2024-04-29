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
sys.path.insert(0, '../common')

from acados_template import AcadosOcp,  AcadosOcpSolver, AcadosOcpBatchSolver
from pendulum_model import export_pendulum_ode_model
import numpy as np
import time
import scipy
import casadi as ca


def setup_ocp(num_threads_in_batch_solve=1, tol=1e-7):

    ocp = AcadosOcp()

    ocp.model = export_pendulum_ode_model()

    Tf = 1.0
    N = 20

    # set dimensions
    ocp.dims.N = N

    Q_mat = 2*np.diag([1e3, 1e3, 1e-2, 1e-2])
    R_mat = 2*np.diag([1e-2])
    cost_W = scipy.linalg.block_diag(Q_mat, R_mat)

    ocp.cost.cost_type = 'NONLINEAR_LS'
    ocp.cost.cost_type_e = 'NONLINEAR_LS'

    ocp.cost.W_e = Q_mat
    ocp.cost.W = cost_W
    ocp.model.cost_y_expr = ca.vertcat(ocp.model.x, ocp.model.u)
    ocp.model.cost_y_expr_e = ocp.model.x
    ocp.cost.yref = np.zeros(ocp.model.cost_y_expr.shape).flatten()
    ocp.cost.yref_e = np.zeros(ocp.model.cost_y_expr_e.shape).flatten()

    # set constraints
    Fmax = 80
    ocp.constraints.lbu = np.array([-Fmax])
    ocp.constraints.ubu = np.array([+Fmax])
    ocp.constraints.idxbu = np.array([0])

    ocp.constraints.x0 = np.array([0.0, np.pi, 0.0, 0.0])

    ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM' # FULL_CONDENSING_QPOASES
    ocp.solver_options.hessian_approx = 'GAUSS_NEWTON' # 'GAUSS_NEWTON', 'EXACT'
    ocp.solver_options.integrator_type = 'IRK'
    ocp.solver_options.nlp_solver_type = 'SQP' # SQP_RTI, SQP
    ocp.solver_options.nlp_solver_tol_stat = tol
    ocp.solver_options.nlp_solver_tol_eq = tol
    ocp.solver_options.nlp_solver_tol_ineq = tol
    ocp.solver_options.nlp_solver_tol_comp = tol
    ocp.solver_options.num_threads_in_batch_solve = num_threads_in_batch_solve

    ocp.solver_options.tf = Tf

    return ocp


def main_sequential(x0, N_sim, tol):

    ocp = setup_ocp(tol=tol)
    solver = AcadosOcpSolver(ocp, verbose=False)

    nx = ocp.dims.nx
    nu = ocp.dims.nu

    simU = np.zeros((N_sim, nu))
    simX = np.zeros((N_sim+1, nx))
    simX[0,:] = x0

    for i in range(N_sim):
        simU[i,:] = solver.solve_for_x0(x0_bar=simX[i, :])
        simX[i+1,:] = solver.get(1, "x")

    return simX, simU


def main_batch(Xinit, simU, tol, num_threads_in_batch_solve=1):

    N_batch = Xinit.shape[0] - 1
    ocp = setup_ocp(num_threads_in_batch_solve, tol)
    batch_solver = AcadosOcpBatchSolver(ocp, N_batch, verbose=False)

    for n in range(N_batch):
        batch_solver.ocp_solvers[n].constraints_set(0, "lbx", Xinit[n])
        batch_solver.ocp_solvers[n].constraints_set(0, "ubx", Xinit[n])

        # set initial guess
        for i in range(ocp.dims.N):
            batch_solver.ocp_solvers[n].set(i, "x", Xinit[n])

    t0 = time.time()
    batch_solver.solve()
    t_elapsed = time.time() - t0
    t_elapsed *= 1000

    print("parallel:  " if num_threads_in_batch_solve else "sequential:", f"{t_elapsed:.3f}ms")

    for n in range(N_batch):
        u = batch_solver.ocp_solvers[n].get(0, "u")
        if not np.linalg.norm(u-simU[n]) < tol*10:
            print(np.linalg.norm(u-simU[n]))
            breakpoint()


if __name__ == "__main__":

    tol = 1e-7
    N_batch = 256
    x0 = np.array([0.0, np.pi, 0.0, 0.0])
    u0 = np.array([0.0])

    simX, simU = main_sequential(x0=x0, N_sim=N_batch, tol=tol)

    main_batch(Xinit=simX, simU=simU, tol=tol, num_threads_in_batch_solve=1)
    main_batch(Xinit=simX, simU=simU, tol=tol, num_threads_in_batch_solve=4)

