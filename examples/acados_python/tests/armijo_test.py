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

from acados_template import AcadosOcp, AcadosOcpSolver, AcadosModel
import numpy as np
from casadi import *
from itertools import product

# tests different globalization settings on simple test problem
#
# min x^2
# s.t. x \in [-10, 10]
#
# with wrong hessian: 1.0

# Settings
TOL = 1e-6

def main():
    # run test cases
    params = {'globalization': ['MERIT_BACKTRACKING', 'FIXED_STEP'],
              'line_search_use_sufficient_descent' : [0, 1],
              'globalization_use_SOC' : [0, 1] }

    keys, values = zip(*params.items())
    for combination in product(*values):
        setting = dict(zip(keys, combination))
        if setting['globalization'] == 'FIXED_STEP' and \
          (setting['globalization_use_SOC'] or setting['line_search_use_sufficient_descent']):
            # skip some equivalent settings
            pass
        else:
            solve_armijo_problem_with_setting(setting)


def solve_armijo_problem_with_setting(setting):
    globalization = setting['globalization']
    line_search_use_sufficient_descent = setting['line_search_use_sufficient_descent']
    globalization_use_SOC = setting['globalization_use_SOC']

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
    model.name = f'armijo_problem'
    ocp.model = model

    # discretization
    Tf = 1
    N = 1
    ocp.dims.N = N
    ocp.solver_options.tf = Tf

    # cost
    ocp.cost.cost_type_e = 'EXTERNAL'
    ocp.model.cost_expr_ext_cost_e = x @ x
    ocp.model.cost_expr_ext_cost_custom_hess_e = 1.0 # 2.0 is the actual hessian

    # constarints
    ocp.constraints.idxbx = np.array([0])
    ocp.constraints.lbx = np.array([-10.0])
    ocp.constraints.ubx = np.array([10.0])

    # options
    ocp.solver_options.qp_solver = 'FULL_CONDENSING_QPOASES' # 'PARTIAL_CONDENSING_HPIPM' # FULL_CONDENSING_QPOASES
    ocp.solver_options.hessian_approx = 'EXACT'
    ocp.solver_options.integrator_type = 'DISCRETE'
    ocp.solver_options.print_level = 0
    ocp.solver_options.tol = TOL
    ocp.solver_options.nlp_solver_type = 'SQP' # SQP_RTI, SQP
    ocp.solver_options.globalization = globalization
    ocp.solver_options.alpha_reduction = 0.9
    ocp.solver_options.line_search_use_sufficient_descent = line_search_use_sufficient_descent
    ocp.solver_options.globalization_use_SOC = globalization_use_SOC
    ocp.solver_options.eps_sufficient_descent = 5e-1
    SQP_max_iter = 200
    ocp.solver_options.qp_solver_iter_max = 400
    ocp.solver_options.nlp_solver_max_iter = SQP_max_iter

    # create solver
    ocp_solver = AcadosOcpSolver(ocp, json_file=f'{model.name}.json')

    # initialize solver
    xinit = np.array([1.0])
    [ocp_solver.set(i, "x", xinit) for i in range(N+1)]

    # get stats
    status = ocp_solver.solve()
    ocp_solver.print_statistics()
    iter = ocp_solver.get_stats('sqp_iter')
    alphas = ocp_solver.get_stats('alpha')[1:]
    qp_iters = ocp_solver.get_stats('qp_iter')
    print(f"acados ocp solver returned status {status}")

    # get solution
    solution = ocp_solver.get(0, "x")
    print(f"found solution {solution}")

    # print summary
    print(f"solved Armijo test problem with settings {setting}")
    print(f"cost function value = {ocp_solver.get_cost()} after {iter} SQP iterations")
    print(f"alphas: {alphas[:iter]}")
    print(f"total number of QP iterations: {sum(qp_iters[:iter])}")

    # compare to analytical solution
    exact_solution = np.array([0.0])
    sol_err = max(np.abs(solution - exact_solution ))
    print(f"error wrt analytical solution {sol_err}")

    # checks
    if ocp.model.cost_expr_ext_cost_custom_hess_e == 1.0:
        if globalization == 'MERIT_BACKTRACKING':
            if sol_err > TOL*1e1:
                raise Exception(f"error of numerical solution wrt exact solution = {sol_err} > tol = {TOL*1e1}")
            else:
                print(f"matched analytical solution with tolerance {TOL}")
            if status != 0:
                raise Exception(f"acados solver returned status {status} != 0.")

            if line_search_use_sufficient_descent == 1:
                if iter > 22:
                    raise Exception(f"acados ocp solver took {iter} iterations." + \
                        "Expected <= 22 iterations for globalized SQP method with aggressive eps_sufficient_descent condition on Armijo test problem.")
            else:
                if iter < 64:
                    raise Exception(f"acados ocp solver took {iter} iterations." + \
                        "Expected > 64 iterations for globalized SQP method without sufficient descent condition on Armijo test problem.")

        elif globalization == 'FIXED_STEP':
            if status != 2:
                raise Exception(f"acados solver returned status {status} != 2. Expected maximum iterations for full-step SQP on Armijo test problem.")
            else:
                print(f"Sucess: Expected maximum iterations for full-step SQP on Armijo test problem.")

    print(f"\n\n----------------------\n")

if __name__ == '__main__':
    main()
