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
# Test problem hs015 from Hock & Schittkowsky test
#
# This is Charlie Vanarets favorite problem since filterSQP it is small, but
# filterSQP requires the feasibility restoration to solve the problem
#
# Problem is in two dimension
#
# min 100*(x2 - x1^2)^2 + (1 - x1)^2
#
# s.t.  x1*x2 >= 1
#       x1 + x2^2 >= 0
#
# Optimal solution is x* = (0.5, 2), and optimal objective is f^*=306.5
# x0 = (-2, 1)


def main():
    # run test cases
    params = {'globalization': ['FUNNEL_L1PEN_LINESEARCH'],
              'nlp_solver_type': ['SQP', 'SQP_WITH_FEASIBLE_QP']}

    keys, values = zip(*params.items())
    for combination in product(*values):
        setting = dict(zip(keys, combination))
        solve_infeasible_linearization(setting)

def solve_infeasible_linearization(setting):

    globalization = setting['globalization']
    nlp_solver_type = setting['nlp_solver_type']

    # create ocp object to formulate the OCP
    ocp = AcadosOcp()

    # set model
    model = AcadosModel()
    x = SX.sym('x', 2)

    # dynamics: identity
    model.disc_dyn_expr = x
    model.x = x
    model.u = SX.sym('u', 0, 0) # [] / None doesnt work
    model.p = []
    model.name = f'hs_015'
    ocp.model = model

    # discretization
    Tf = 1
    N = 1
    ocp.solver_options.N_horizon = N
    ocp.solver_options.tf = Tf

    # cost
    ocp.cost.cost_type_e = 'EXTERNAL'
    ocp.model.cost_expr_ext_cost_e = 100*(x[1] - x[0]**2)**2 + (1 - x[0])**2

    # constraints
    ocp.model.con_h_expr_0 = vertcat(x[0]*x[1], x[0] + x[1]**2)
    ocp.constraints.lh_0 = np.array([1.0, 0.0])
    ocp.constraints.uh_0 = np.array([ACADOS_INFTY, ACADOS_INFTY])

    # add bounds on x
    ocp.constraints.idxbx_0 = np.array(range(1))
    ocp.constraints.ubx_0 = np.array([0.5])
    ocp.constraints.lbx_0 = np.array([-ACADOS_INFTY])

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
    ocp.solver_options.nlp_solver_max_iter = 100
    ocp.solver_options.search_direction_mode = "BYRD_OMOJOKUN"
    ocp.solver_options.use_constraint_hessian_in_feas_qp = False

    ocp.code_export_directory = f'c_generated_code_{model.name}'
    ocp_solver = AcadosOcpSolver(ocp, json_file=f'{model.name}.json')

    # initialize solver
    xinit = np.array([-2, 1])
    [ocp_solver.set(i, "x", xinit) for i in range(N+1)]

    # solve
    status = ocp_solver.solve()

    # get solution
    solution = ocp_solver.get(0, "x")
    cost = ocp_solver.get_cost()

    # compare to analytical solution
    exact_solution = np.array([0.5, 2])
    optimal_objective = 306.5

    if ocp.solver_options.nlp_solver_type == 'SQP':
        assert status == 4, "As expected the standard SQP method should not be able to solve hs015!"
    if ocp.solver_options.nlp_solver_type == 'SQP_WITH_FEASIBLE_QP':
        assert status == 0, "SQP_WITH_FEASIBLE_QP method should converge!"
        assert np.allclose(solution, exact_solution), f"Found optimal solution should be (0.5,2), got {solution}!"
        assert np.allclose(cost, optimal_objective), f"Found cost should be 306.5, got {cost}!"

if __name__ == '__main__':
    main()
