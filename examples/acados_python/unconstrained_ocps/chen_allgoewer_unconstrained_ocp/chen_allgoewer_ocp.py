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

from acados_template import AcadosOcp, AcadosOcpSolver
from chen_allgoewer_system_model import export_chen_allgoewer_model
import numpy as np
from utils import plot_trajectory

def main(plot_solution = False):

    # The flag denotes, if the problem should be transformed into a feasibility
    # problem, or if the unconstrained OCP should be solved.
    SOLVE_FEASIBILITY_PROBLEM = True

    # create ocp object to formulate the OCP
    ocp = AcadosOcp()

    # set model
    model = export_chen_allgoewer_model(use_SX=False)
    ocp.model = model

    Tf = 5.0
    nx = model.x.rows()
    nu = model.u.rows()
    N = 20
    M = 10 # Needed for integrator

    # set dimensions
    ocp.dims.N = N

    # set cost
    Q_mat = np.array([[0.5, 0], [0, 0.5]])
    R_mat = np.array([[0.8]])
    P_mat = np.array([[10.0, 0], [0, 10.0]])

    u_max = 1
    tau = 100
    # the 'EXTERNAL' cost type can be used to define general cost terms
    # NOTE: This leads to additional (exact) hessian contributions when using GAUSS_NEWTON hessian.

    if not SOLVE_FEASIBILITY_PROBLEM:
        ocp.cost.cost_type = 'EXTERNAL'
        ocp.cost.cost_type_e = 'EXTERNAL'
        tau_beta = tau * np.fmax(0, model.u-u_max)**2 + tau * np.fmin(0, model.u + u_max)**2
        ocp.model.cost_expr_ext_cost = 0.5 * model.x.T @ Q_mat @ model.x + 0.5 * model.u.T @ R_mat @ model.u + tau_beta
        ocp.model.cost_expr_ext_cost_e = 0.5 * model.x.T @ P_mat @ model.x

    # set constraints
    ocp.constraints.x0 = np.array([0.42, 0.45])

    if SOLVE_FEASIBILITY_PROBLEM:
        # Path constraints on control
        u_max = 1.5
        ocp.constraints.lbu = np.array([-u_max])
        ocp.constraints.ubu = np.array([+u_max])
        ocp.constraints.idxbu = np.array([0])

        # Terminal constraints
        ocp.constraints.lbx_e = np.array([0.0, 0.03])
        ocp.constraints.ubx_e  = np.array([0.0, 0.03])
        ocp.constraints.idxbx_e = np.arange(nx)


    # set options
    ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM'
    ocp.solver_options.qp_solver_cond_N = N
    ocp.solver_options.hessian_approx = 'GAUSS_NEWTON'
    ocp.solver_options.integrator_type = 'ERK'
    ocp.solver_options.sim_method_num_steps = M
    ocp.solver_options.print_level = 1
    ocp.solver_options.nlp_solver_type = 'DDP'
    ocp.solver_options.nlp_solver_max_iter = 10
    ocp.solver_options.globalization = 'MERIT_BACKTRACKING'
    ocp.solver_options.with_adaptive_levenberg_marquardt = True

    # set prediction horizon
    ocp.solver_options.tf = Tf

    if SOLVE_FEASIBILITY_PROBLEM:
        ocp.translate_to_feasibility_problem(parametric_x0=True)

    ocp_solver = AcadosOcpSolver(ocp, json_file = 'chen_allgoewer_acados.json')

    for i in range(N):
        ocp_solver.cost_set(i, "scaling", 1.0)

    sol_X = np.zeros((N+1, nx))
    sol_U = np.zeros((N, nu))

    # Load and set the initial guess
    with open('chen_allgoewer_initial_guess.npy', 'rb') as f:
        X_init = np.load(f)
        U_init = np.load(f)

    initial_conditions = [np.array([0.42, 0.45]), np.array([0.42, 0.5])]
    for initial_condition in initial_conditions:

        # Initial guess
        for i in range(N):
            ocp_solver.set(i, "x", X_init[:,i])
            ocp_solver.set(i, "u", U_init[:,i])
        ocp_solver.set(N, "x", X_init[:,N])

        ocp_solver.set(0, 'p', initial_condition)

        status = ocp_solver.solve()

        if status != 0:
            raise Exception(f'acados returned status {status}.')

        iter = ocp_solver.get_stats('nlp_iter')
        assert iter in [4,5], "DDP Solver should converge within 4 or 5 iterations!"

        # get solution
        for i in range(N):
            sol_X[i,:] = ocp_solver.get(i, "x")
            sol_U[i,:] = ocp_solver.get(i, "u")
        sol_X[N,:] = ocp_solver.get(N, "x")

        assert np.allclose(sol_X[0,:].squeeze(), initial_condition), "Initial condition does not coincide with parameter!"

        if plot_solution:
            plot_trajectory([X_init, sol_X.T], ["Initial guess", "Solution"])

if __name__ == '__main__':
    main(plot_solution = False)
