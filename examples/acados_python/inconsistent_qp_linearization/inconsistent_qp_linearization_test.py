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
    params = {'nlp_solver_type': ['SQP_WITH_FEASIBLE_QP'],
              'search_direction_mode':['NOMINAL_QP'],
              'max_iter':[20],
              'init_iterate': [np.array([-0.001])]}
            #   'init_iterate': [np.array([-0.001]), np.array([0.5])]}

    keys, values = zip(*params.items())
    for combination in product(*values):
        setting = dict(zip(keys, combination))
        test_convergence_of_solver(setting)

def test_nominal_qp():
    params = {'nlp_solver_type': 'SQP_WITH_FEASIBLE_QP',
              'search_direction_mode':'NOMINAL_QP',
              'max_iter':1,
              'init_iterate': np.array([-0.001])}

    N = 1
    ocp, ocp_solver = create_solver(params)
    xinit = params['init_iterate']

    # initialize solver
    [ocp_solver.set(i, "x", xinit) for i in range(N+1)]

    # solve
    _ = ocp_solver.solve()

    iter0 = ocp_solver.get_iterate(0)
    iter1 = ocp_solver.get_iterate(1)

    # solution is d = -1999.9995
    d = iter1.x_traj[0] - iter0.x_traj[0]
    assert np.allclose(d, -1999.9995), f"Solution should be -1999.9995, got {d}"

def test_byrd_omojokun_qps():
    params = {'nlp_solver_type': 'SQP_WITH_FEASIBLE_QP',
              'search_direction_mode':'BYRD_OMOJOKUN',
              'max_iter':1,
              'init_iterate': np.array([-0.001])}

    N = 1
    ocp, ocp_solver = create_solver(params)
    xinit = params['init_iterate']

    # initialize solver
    [ocp_solver.set(i, "x", xinit) for i in range(N+1)]

    # solve
    status = ocp_solver.solve()
    last_qp = ocp_solver.get_last_qp()
    last_relaxed_qp = ocp_solver.get_last_relaxed_qp()

    iter0 = ocp_solver.get_iterate(0)
    iter1 = ocp_solver.get_iterate(1)

    # feasibility QP solution should be (d = -10, s = 3.98)
    # here should be a test of the bounds

    # nominal QP solution should be d= -10
    d = iter1.x_traj[0] - iter0.x_traj[0]
    assert np.allclose(d, -10), f"Solution should be -10, got {d}"


def create_solver_opts(N=1,
                       Tf=1,
                       nlp_solver_type ='SQP_WITH_FEASIBLE_QP',
                       max_iter: int = 20,
                       search_direction_mode='NOMINAL_QP'):

    solver_options = AcadosOcp().solver_options

    # set options
    solver_options.N_horizon = N
    solver_options.tol = 1e-6
    solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM'
    solver_options.qp_solver_cond_N = N
    solver_options.qp_solver_iter_max = 1000
    solver_options.qp_tol = 1e-9
    solver_options.qp_solver_mu0 = 1e4
    solver_options.hessian_approx = 'EXACT'
    solver_options.regularize_method = 'MIRROR'
    solver_options.integrator_type = 'DISCRETE'
    solver_options.print_level = 1
    solver_options.nlp_solver_type = nlp_solver_type
    solver_options.globalization = 'FUNNEL_L1PEN_LINESEARCH'
    solver_options.globalization_full_step_dual = True
    solver_options.globalization_alpha_min = 1e-15
    solver_options.nlp_solver_max_iter = max_iter
    solver_options.search_direction_mode = search_direction_mode
    solver_options.use_constraint_hessian_in_feas_qp = False
    solver_options.store_iterates = True

    # set prediction horizon
    solver_options.tf = Tf

    return solver_options

def create_solver(setting):
    print(setting)

    nlp_solver_type = setting['nlp_solver_type']
    search_direction_mode = setting['search_direction_mode']
    max_iter = setting['max_iter']

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

    ocp.solver_options = create_solver_opts(N, Tf, nlp_solver_type, max_iter, search_direction_mode)
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
    ocp_solver.print_statistics()

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
            assert np.allclose(solution, exact_solution), "Optimal solution should be -2!"
        elif np.allclose(xinit, np.array([0.5])):
            assert status == 3, "SQP with feasible QP should converge to infeasible stationary point with min step!"
            assert np.allclose(solution, infeasible_solution), "Optimal solution should be 1!"

if __name__ == '__main__':
    test_nominal_qp()
    test_byrd_omojokun_qps()
    main()
