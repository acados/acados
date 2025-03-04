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

from acados_template import AcadosOcp, AcadosOcpSolver, AcadosModel, ACADOS_INFTY
import numpy as np
from casadi import *
from matplotlib import pyplot as plt
from itertools import product

# Problem with infeasible linearization
#
# min -x
#
# s.t.  x <= 1
#       x^2 >= 4
#
# The optimal solution is x* = -2 with lam1 = 0 and lam2 = 0.25,
#
# but started from x > 0, SQP solver converges to
# infeasible point x^ = 1 with lam1 = 2 and lam2 = 1 for the infeasibility problem.


def main():
    # run test cases
    params = {'globalization': ['FUNNEL_L1PEN_LINESEARCH'],
              'nlp_solver_type': ['SQP', 'SQP_WITH_FEASIBLE_QP'],
              'init_iterate': [np.array([-0.001]), np.array([0.5])]}

    keys, values = zip(*params.items())
    for combination in product(*values):
        setting = dict(zip(keys, combination))
        test_convergence_of_solver(setting)


def create_solver(setting):

    globalization = setting['globalization']
    nlp_solver_type = setting['nlp_solver_type']

    # create ocp object to formulate the OCP
    ocp = AcadosOcp()

    # set model
    model = AcadosModel()
    x = SX.sym('x')

    # dynamics: identity
    model.disc_dyn_expr = x
    model.x = x
    model.u = SX.sym('u', 0, 0) # [] / None doesnt work
    model.p = []
    model.name = f'inconsistent_qp_linearization'
    ocp.model = model

    # discretization
    Tf = 1
    N = 1
    ocp.solver_options.N_horizon = N
    ocp.solver_options.tf = Tf

    # cost
    ocp.cost.cost_type_0 = 'EXTERNAL'
    ocp.model.cost_expr_ext_cost_0 = -model.x[0]

    # constraints
    ocp.model.con_h_expr_0 = x**2
    ocp.constraints.lh_0 = np.array([4.0])
    ocp.constraints.uh_0 = np.array([ACADOS_INFTY])

    # add bounds on x
    nx = 1
    ocp.constraints.idxbx_0 = np.array(range(nx))
    ocp.constraints.lbx_0 = -ACADOS_INFTY * np.ones((nx))
    ocp.constraints.ubx_0 = 1 * np.ones((nx))

    # set options
    ocp.solver_options.tol = 1e-6
    ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM'
    ocp.solver_options.qp_solver_cond_N = N
    ocp.solver_options.qp_solver_iter_max = 1000
    ocp.solver_options.qp_tol = 1e-9
    ocp.solver_options.qp_solver_mu0 = 1e4
    ocp.solver_options.hessian_approx = 'EXACT'
    ocp.solver_options.regularize_method = 'MIRROR'
    ocp.solver_options.integrator_type = 'DISCRETE'
    ocp.solver_options.print_level = 1
    ocp.solver_options.nlp_solver_type = nlp_solver_type
    ocp.solver_options.globalization = globalization
    ocp.solver_options.globalization_full_step_dual = True
    ocp.solver_options.globalization_alpha_min = 1e-15
    ocp.solver_options.nlp_solver_max_iter = 20
    ocp.solver_options.search_direction_mode = "BYRD_OMOJOKUN"
    ocp.solver_options.use_exact_hessian_in_feas_qp = False
    ocp_solver = AcadosOcpSolver(ocp, json_file=f'{model.name}.json')

    return ocp, ocp_solver

def test_convergence_of_solver(setting):

    N = 1
    ocp, ocp_solver = create_solver(setting)
    xinit = setting['init_iterate']

    # initialize solver
    [ocp_solver.set(i, "x", xinit) for i in range(N+1)]

    # solve
    status = ocp_solver.solve()

    # get solution
    solution = ocp_solver.get(0, "x")

    # compare to analytical solution
    exact_solution = np.array([-2.0])

    infeasible_solution = np.array([1.0])

    if ocp.solver_options.nlp_solver_type == 'SQP':
        if np.allclose(xinit, np.array([-0.001])):
            assert status == 0, "Standard SQP should be able to solve the problem!"
            assert np.allclose(solution, exact_solution), "Optimal solution should be -2!"
        elif np.allclose(xinit, np.array([0.5])):
            assert status == 4, "QP subproblem should get infeasible for standard SQP!"
    if ocp.solver_options.nlp_solver_type == 'SQP_WITH_FEASIBLE_QP':
        if np.allclose(xinit, np.array([-0.001])):
            assert status == 0, "SQP with feasible QP should be able to solve the problem!"
            assert np.allclose(solution, exact_solution), "Optimal soluvscodetion should be -2!"
        elif np.allclose(xinit, np.array([0.5])):
            assert status == 3, "SQP with feasible QP should converge to infeasible stationary point with min step!"
            assert np.allclose(solution, infeasible_solution), "Optimal solution should be 1!"

if __name__ == '__main__':
    main()
