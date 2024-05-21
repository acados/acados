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
from rockit_hello_world_model import export_rockit_hello_world_model
import numpy as np

def main():

    # The flag denotes, if the problem should be transformed into a feasibility
    # problem, or if the unconstrained OCP should be solved.
    SOLVE_FEASIBILITY_PROBLEM = True
    DYNAMICALLY_FEASIBLE_INITIAL_GUESS = True

    # create ocp object to formulate the OCP
    ocp = AcadosOcp()

    # set model
    model = export_rockit_hello_world_model()
    ocp.model = model

    Tf = 1.0
    nx = model.x.rows()
    nu = model.u.rows()
    N = 2

    # set dimensions
    ocp.dims.N = N

    if not SOLVE_FEASIBILITY_PROBLEM:
        # set cost
        # the 'EXTERNAL' cost type can be used to define general cost terms
        # NOTE: This leads to additional (exact) hessian contributions when using GAUSS_NEWTON hessian.
        ocp.cost.cost_type = 'EXTERNAL'
        ocp.cost.cost_type_e = 'EXTERNAL'
        ocp.model.cost_expr_ext_cost = model.x.T @ model.x + model.u.T @ model.u
        ocp.model.cost_expr_ext_cost_e = model.x[0]* model.x[0]

    # set constraints
    ocp.constraints.x0 = np.array([0.0, 1.0])

    if SOLVE_FEASIBILITY_PROBLEM:
        # Path constraints on control
        u_max = 1.0
        ocp.constraints.lbu = np.array([-u_max])
        ocp.constraints.ubu = np.array([+u_max])
        ocp.constraints.idxbu = np.array([0])

        # Path constraint on x1
        x_min = -0.25
        ocp.constraints.lbx = np.array([x_min])
        ocp.constraints.ubx = np.array([1e12])
        ocp.constraints.idxbx = np.array([0])

        ocp.constraints.lbx_e = np.array([x_min])
        ocp.constraints.ubx_e  = np.array([1e12])
        ocp.constraints.idxbx_e = np.array([0])

    # set options
    ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM'
    ocp.solver_options.qp_solver_cond_N = N
    ocp.solver_options.hessian_approx = 'GAUSS_NEWTON'
    ocp.solver_options.integrator_type = 'ERK'
    ocp.solver_options.sim_method_num_steps = 1
    ocp.solver_options.print_level = 1
    ocp.solver_options.nlp_solver_type = 'DDP'
    ocp.solver_options.nlp_solver_max_iter = 100
    ocp.solver_options.globalization = 'MERIT_BACKTRACKING'
    ocp.solver_options.with_adaptive_levenberg_marquardt = True

    # set prediction horizon
    ocp.solver_options.tf = Tf

    if SOLVE_FEASIBILITY_PROBLEM:
        ocp.translate_to_feasibility_problem()

    ocp_solver = AcadosOcpSolver(ocp, json_file = 'rockit_hello_world.json')

    # acados multiplies all stage costs with the time step by default
    for i in range(N):
        ocp_solver.cost_set(i, "scaling", 1.0)
    sol_X = np.zeros((N+1, nx))
    sol_U = np.zeros((N, nu))

    if DYNAMICALLY_FEASIBLE_INITIAL_GUESS:
        # Load the initial guess
        with open('rockit_hello_world_initial_guess.npy', 'rb') as f:
            X_init = np.load(f)
            U_init = np.load(f)
    else:
        X_init = np.zeros((nx, N+1))
        x1s = (1/N)*np.ones((1,N+1))
        X_init[0,:] = x1s

        U_init = np.zeros((nu, N))
        us = np.linspace(0, 1.0/N, N)
        U_init[0,:] = us

    # Set initial guess
    for i in range(N):
        ocp_solver.set(i, "x", X_init[:,i])
        ocp_solver.set(i, "u", U_init[:,i])
    ocp_solver.set(N, "x", X_init[:,N])

    # Solve the problem
    status = ocp_solver.solve()

    ocp_solver.print_statistics()

    iter = ocp_solver.get_stats('nlp_iter')
    assert iter <= 6, "DDP Solver should converge within 6 iterations!"

    if status != 0:
        raise Exception(f'acados returned status {status}.')

    # get solution
    for i in range(N):
        sol_X[i,:] = ocp_solver.get(i, "x")
        sol_U[i,:] = ocp_solver.get(i, "u")
    sol_X[N,:] = ocp_solver.get(N, "x")

if __name__ == '__main__':
    main()