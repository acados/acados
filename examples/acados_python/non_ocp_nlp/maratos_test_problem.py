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
from matplotlib import pyplot as plt
from itertools import product
# Simplest NLP with Maratos effect
#
# min x_1
#
# s.t. x_1^2 + x_2^2 = 1

# Settings
PLOT = False
TOL = 1e-6

def main():
    # run test cases
    params = {'globalization': ['MERIT_BACKTRACKING', 'FIXED_STEP', 'FUNNEL_L1PEN_LINESEARCH'],
              'line_search_use_sufficient_descent' : [0, 1],
              'globalization_use_SOC' : [0, 1] }

    keys, values = zip(*params.items())
    for combination in product(*values):
        setting = dict(zip(keys, combination))
        if setting['globalization'] == 'FIXED_STEP' and \
          (setting['globalization_use_SOC'] or setting['line_search_use_sufficient_descent']):
            # skip some equivalent settings
            pass
        elif setting['globalization'] == 'FUNNEL_L1PEN_LINESEARCH' and \
          (setting['globalization_use_SOC'] or setting['line_search_use_sufficient_descent']):
            # skip some equivalent settings
            pass
        else:
            solve_maratos_problem_with_setting(setting)
            # exit(1)
            # pass
    # setting = {"globalization": "MERIT_BACKTRACKING",
            #   "line_search_use_sufficient_descent": 0,
            #   "globalization_use_SOC": 1}
    # solve_maratos_problem_with_setting(setting)


def solve_maratos_problem_with_setting(setting):

    globalization = setting['globalization']
    line_search_use_sufficient_descent = setting['line_search_use_sufficient_descent']
    globalization_use_SOC = setting['globalization_use_SOC']

    print(f"running maratos test problem with settings {setting}")

    # create ocp object to formulate the OCP
    ocp = AcadosOcp()

    # set model
    x1 = SX.sym('x1')
    x2 = SX.sym('x2')
    x = vertcat(x1, x2)

    # dynamics: identity
    ocp.model.x = x
    ocp.model.name = f'maratos_problem'

    # discretization
    N = 1
    ocp.solver_options.N_horizon = N

    if N == 0:
        # cost
        ocp.cost.cost_type_e = 'EXTERNAL'
        ocp.model.cost_expr_ext_cost_e = x1

        # constraints
        ocp.model.con_h_expr_e = x1 ** 2 + x2 ** 2
        ocp.constraints.lh_e = np.array([1.0])
        ocp.constraints.uh_e = np.array([1.0])
    elif N == 1:
        # dynamics: identity
        ocp.model.disc_dyn_expr = x
        ocp.model.u = SX.sym('u', 0, 0) # [] / None doesnt work

        # discretization
        ocp.solver_options.tf = 1.0
        ocp.solver_options.integrator_type = 'DISCRETE'

        # cost
        ocp.cost.cost_type_e = 'EXTERNAL'
        ocp.model.cost_expr_ext_cost_e = x1

        # constarints
        ocp.model.con_h_expr_0 = x1 ** 2 + x2 ** 2
        ocp.constraints.lh_0 = np.array([1.0])
        ocp.constraints.uh_0 = np.array([1.0])
    else:
        raise NotImplementedError('N > 1 not implemented')
    # # soften
    # ocp.constraints.idxsh_e = np.array([0])
    # ocp.cost.zl_e = 1e5 * np.array([1])
    # ocp.cost.zu_e = 1e5 * np.array([1])
    # ocp.cost.Zl_e = 1e5 * np.array([1])
    # ocp.cost.Zu_e = 1e5 * np.array([1])

    # add bounds on x
    # nx = 2
    # ocp.constraints.idxbx_e = np.array(range(nx))
    # ocp.constraints.lbx_e = -2 * np.ones((nx))
    # ocp.constraints.ubx_e = 2 * np.ones((nx))

    # set options
    ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM' # TODO: check difference wrt FULL_CONDENSING
    ocp.solver_options.hessian_approx = 'EXACT'
    # ocp.solver_options.print_level = 2
    ocp.solver_options.tol = TOL
    ocp.solver_options.nlp_solver_type = 'SQP' # SQP_RTI, SQP
    ocp.solver_options.levenberg_marquardt = 1e-1 # / (N+1)
    SQP_max_iter = 300
    ocp.solver_options.qp_solver_iter_max = 400
    ocp.solver_options.qp_tol = 5e-7
    ocp.solver_options.regularize_method = 'MIRROR'
    # ocp.solver_options.exact_hess_constr = 0
    ocp.solver_options.globalization = globalization
    ocp.solver_options.globalization_alpha_min = 1e-2
    ocp.solver_options.globalization_line_search_use_sufficient_descent = line_search_use_sufficient_descent
    ocp.solver_options.globalization_use_SOC = globalization_use_SOC
    ocp.solver_options.globalization_eps_sufficient_descent = 1e-1
    ocp.solver_options.store_iterates = True

    ocp.solver_options.nlp_solver_max_iter = SQP_max_iter
    ocp_solver = AcadosOcpSolver(ocp, verbose=False)

    # initialize solver
    rad_init = 0.1 #0.1 #np.pi / 4
    xinit = np.array([np.cos(rad_init), np.sin(rad_init)])
    # xinit = np.array([0.82120912, 0.58406911])
    [ocp_solver.set(i, "x", xinit) for i in range(N+1)]

    # solve
    ocp_solver.solve()
    ocp_solver.print_statistics()
    iter = ocp_solver.get_stats('sqp_iter')
    alphas = ocp_solver.get_stats('alpha')[1:]
    qp_iters = ocp_solver.get_stats('qp_iter')
    residuals = ocp_solver.get_stats('statistics')[1:5,1:iter]
    iterates = ocp_solver.get_iterates()
    x_iterates = iterates.as_array('x')
    if N > 0:
        xdiff = x_iterates[:, 0, :] - x_iterates[:, 1, :]
        xdiff = np.linalg.norm(xdiff, axis=1)
        print(f"xdiff = {xdiff}")

    # get solution
    solution = ocp_solver.get(0, "x")

    # print summary
    print(f"solved Maratos test problem with settings {setting}")
    print(f"cost function value = {ocp_solver.get_cost()} after {iter} SQP iterations")
    print(f"alphas: {alphas[:iter]}")
    print(f"total number of QP iterations: {sum(qp_iters[:iter])}")
    max_infeasibility = np.max(residuals[1:3])
    print(f"max infeasibility: {max_infeasibility}")

    # compare to analytical solution
    exact_solution = np.array([-1, 0])
    sol_err = max(np.abs(solution - exact_solution ))

    # checks
    if sol_err > TOL*1e1:
        print(f"error of numerical solution wrt exact solution = {sol_err} > tol = {TOL*1e1}")
    else:
        print(f"matched analytical solution with tolerance {TOL}")

    try:
        if globalization == 'FIXED_STEP':
            if max_infeasibility < 5.0:
                raise Exception(f"Expected max_infeasibility > 5.0 when using full step SQP on Maratos problem")
            if iter != 10:
                raise Exception(f"Expected 10 SQP iterations when using full step SQP on Maratos problem, got {iter}")
            if any(alphas[:iter] != 1.0):
                raise Exception(f"Expected all alphas = 1.0 when using full step SQP on Maratos problem")
        elif globalization == 'MERIT_BACKTRACKING':
            if max_infeasibility > 0.5:
                raise Exception(f"Expected max_infeasibility < 0.5 when using globalized SQP on Maratos problem")
            elif globalization_use_SOC == 0:
                if iter not in range(56, 61):
                    raise Exception(f"Expected 56 to 60 SQP iterations when using globalized SQP without SOC on Maratos problem, got {iter}")
            elif line_search_use_sufficient_descent == 1:
                if iter not in range(29, 37):
                    # NOTE: got 29 locally and 36 on Github actions.
                    # On Github actions the inequality constraint was numerically violated in the beginning.
                    # This leads to very different behavior, since the merit gradient is so different.
                    # Github actions:  merit_grad = -1.669330e+00, merit_grad_cost = -1.737950e-01, merit_grad_dyn = 0.000000e+00, merit_grad_ineq = -1.495535e+00
                    # Jonathan Laptop: merit_grad = -1.737950e-01, merit_grad_cost = -1.737950e-01, merit_grad_dyn = 0.000000e+00, merit_grad_ineq = 0.000000e+00
                    raise Exception(f"Expected SQP iterations in range(29, 37) when using globalized SQP with SOC on Maratos problem, got {iter}")
            else:
                if iter != 16:
                    raise Exception(f"Expected 16 SQP iterations when using globalized SQP with SOC on Maratos problem, got {iter}")
        elif globalization == 'FUNNEL_L1PEN_LINESEARCH':
            if iter > 12:
                raise Exception(f"Expected not more than 12 SQP iterations when using Funnel Method SQP, got {iter}")

    except Exception as inst:
        if N == 0:
            print(f"Exceptions in this file are tailored to formulation with N=1, difference should be investigated.")
            print(f"got Exception {inst} in test with settings {setting}")
        else:
            raise(inst)

    if PLOT:
        plt.figure()
        axs = plt.plot(solution[0], solution[1], 'x', label='solution')

        cm = plt.cm.get_cmap('RdYlBu')
        axs = plt.scatter(x_iterates[:iter+1, 0, 0], x_iterates[:iter+1, 0, 1], c=range(iter+1), s=35, cmap=cm, label='iterates')
        plt.colorbar(axs)

        ts = np.linspace(0,2*np.pi,100)
        plt.plot(1 * np.cos(ts)+0,1 * np.sin(ts)-0, 'r')
        plt.axis('square')
        plt.legend()
        plt.title(f"Maratos problem with N = {N}, x formulation, SOC {globalization_use_SOC}")
        plt.show()

    print(f"\n\n----------------------\n")

if __name__ == '__main__':
    main()
