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
# Simplest NLP with Marathos effect
#
# min x_1
#
# s.t. x_1^2 + x_2^2 = 1

# Settings
PLOT = False
FOR_LOOPING = False # call solver in for loop to get all iterates
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
            solve_marathos_problem_with_setting(setting)


def solve_marathos_problem_with_setting(setting):

    globalization = setting['globalization']
    line_search_use_sufficient_descent = setting['line_search_use_sufficient_descent']
    globalization_use_SOC = setting['globalization_use_SOC']

    # create ocp object to formulate the OCP
    ocp = AcadosOcp()

    # set model
    model = AcadosModel()
    x1 = SX.sym('x1')
    x2 = SX.sym('x2')
    x = vertcat(x1, x2)

    # dynamics: identity
    model.disc_dyn_expr = x
    model.x = x
    model.u = SX.sym('u', 0, 0) # [] / None doesnt work
    model.p = []
    model.name = f'marathos_problem'
    ocp.model = model

    # discretization
    Tf = 1
    N = 1
    ocp.dims.N = N
    ocp.solver_options.tf = Tf

    # cost
    ocp.cost.cost_type_e = 'EXTERNAL'
    ocp.model.cost_expr_ext_cost_e = x1

    # constarints
    ocp.model.con_h_expr_0 = x1 ** 2 + x2 ** 2
    ocp.constraints.lh_0 = np.array([1.0])
    ocp.constraints.uh_0 = np.array([1.0])
    # # soften
    # ocp.constraints.idxsh = np.array([0])
    # ocp.cost.zl = 1e5 * np.array([1])
    # ocp.cost.zu = 1e5 * np.array([1])
    # ocp.cost.Zl = 1e5 * np.array([1])
    # ocp.cost.Zu = 1e5 * np.array([1])

    # add bounds on x
    # nx = 2
    # ocp.constraints.idxbx_0 = np.array(range(nx))
    # ocp.constraints.lbx_0 = -2 * np.ones((nx))
    # ocp.constraints.ubx_0 = 2 * np.ones((nx))

    # set options
    ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM' # FULL_CONDENSING_QPOASES
    # PARTIAL_CONDENSING_HPIPM, FULL_CONDENSING_QPOASES, FULL_CONDENSING_HPIPM,
    # PARTIAL_CONDENSING_QPDUNES, PARTIAL_CONDENSING_OSQP
    ocp.solver_options.hessian_approx = 'EXACT'
    ocp.solver_options.integrator_type = 'DISCRETE'
    # ocp.solver_options.print_level = 1
    ocp.solver_options.tol = TOL
    ocp.solver_options.nlp_solver_type = 'SQP' # SQP_RTI, SQP
    ocp.solver_options.globalization = globalization
    ocp.solver_options.alpha_min = 1e-2
    # ocp.solver_options.__initialize_t_slacks = 0
    # ocp.solver_options.regularize_method = 'CONVEXIFY'
    ocp.solver_options.levenberg_marquardt = 1e-1
    # ocp.solver_options.print_level = 2
    SQP_max_iter = 300
    ocp.solver_options.qp_solver_iter_max = 400
    ocp.solver_options.regularize_method = 'MIRROR'
    # ocp.solver_options.exact_hess_constr = 0
    ocp.solver_options.line_search_use_sufficient_descent = line_search_use_sufficient_descent
    ocp.solver_options.globalization_use_SOC = globalization_use_SOC
    ocp.solver_options.eps_sufficient_descent = 1e-1
    ocp.solver_options.qp_tol = 5e-7

    if FOR_LOOPING: # call solver in for loop to get all iterates
        ocp.solver_options.nlp_solver_max_iter = 1
        ocp_solver = AcadosOcpSolver(ocp, json_file=f'{model.name}.json')
    else:
        ocp.solver_options.nlp_solver_max_iter = SQP_max_iter
        ocp_solver = AcadosOcpSolver(ocp, json_file=f'{model.name}.json')

    # initialize solver
    rad_init = 0.1 #0.1 #np.pi / 4
    xinit = np.array([np.cos(rad_init), np.sin(rad_init)])
    # xinit = np.array([0.82120912, 0.58406911])
    [ocp_solver.set(i, "x", xinit) for i in range(N+1)]

    # solve
    if FOR_LOOPING: # call solver in for loop to get all iterates
        iterates = np.zeros((SQP_max_iter+1, 2))
        iterates[0, :] = xinit
        alphas = np.zeros((SQP_max_iter,))
        qp_iters = np.zeros((SQP_max_iter,))
        iter = SQP_max_iter
        residuals = np.zeros((4, SQP_max_iter))

        # solve
        for i in range(SQP_max_iter):
            status = ocp_solver.solve()
            ocp_solver.print_statistics() # encapsulates: stat = ocp_solver.get_stats("statistics")
            # print(f'acados returned status {status}.')
            iterates[i+1, :] = ocp_solver.get(0, "x")
            if status in [0, 4]:
                iter = i
                break
            alphas[i] = ocp_solver.get_stats('alpha')[1]
            qp_iters[i] = ocp_solver.get_stats('qp_iter')[1]
            residuals[:, i] = ocp_solver.get_stats('residuals')

    else:
        ocp_solver.solve()
        ocp_solver.print_statistics()
        iter = ocp_solver.get_stats('sqp_iter')
        alphas = ocp_solver.get_stats('alpha')[1:]
        qp_iters = ocp_solver.get_stats('qp_iter')
        residuals = ocp_solver.get_stats('statistics')[1:5,1:iter]

    # get solution
    solution = ocp_solver.get(0, "x")

    # print summary
    print(f"solved Marathos test problem with settings {setting}")
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
        raise Exception(f"error of numerical solution wrt exact solution = {sol_err} > tol = {TOL*1e1}")
    else:
        print(f"matched analytical solution with tolerance {TOL}")

    try:
        if globalization == 'FIXED_STEP':
            # import pdb; pdb.set_trace()
            if max_infeasibility < 5.0:
                raise Exception(f"Expected max_infeasibility > 5.0 when using full step SQP on Marathos problem")
            if iter != 10:
                raise Exception(f"Expected 10 SQP iterations when using full step SQP on Marathos problem, got {iter}")
            if any(alphas[:iter] != 1.0):
                raise Exception(f"Expected all alphas = 1.0 when using full step SQP on Marathos problem")
        elif globalization == 'MERIT_BACKTRACKING':
            if max_infeasibility > 0.5:
                raise Exception(f"Expected max_infeasibility < 0.5 when using globalized SQP on Marathos problem")
            if globalization_use_SOC == 0:
                if FOR_LOOPING and iter != 57:
                    raise Exception(f"Expected 57 SQP iterations when using globalized SQP without SOC on Marathos problem, got {iter}")
            elif line_search_use_sufficient_descent == 1:
                if iter not in range(29, 37):
                    # NOTE: got 29 locally and 36 on Github actions.
                    # On Github actions the inequality constraint was numerically violated in the beginning.
                    # This leads to very different behavior, since the merit gradient is so different.
                    # Github actions:  merit_grad = -1.669330e+00, merit_grad_cost = -1.737950e-01, merit_grad_dyn = 0.000000e+00, merit_grad_ineq = -1.495535e+00
                    # Jonathan Laptop: merit_grad = -1.737950e-01, merit_grad_cost = -1.737950e-01, merit_grad_dyn = 0.000000e+00, merit_grad_ineq = 0.000000e+00
                    raise Exception(f"Expected SQP iterations in range(29, 37) when using globalized SQP with SOC on Marathos problem, got {iter}")
            else:
                if iter != 12:
                    raise Exception(f"Expected 12 SQP iterations when using globalized SQP with SOC on Marathos problem, got {iter}")
    except Exception as inst:
        if FOR_LOOPING and globalization == "MERIT_BACKTRACKING":
            print("\nAcados globalized OCP solver behaves different when for looping due to different merit function weights.",
            "Following exception is not raised\n")
            print(inst, "\n")
        else:
            raise(inst)

    if PLOT:
        plt.figure()
        axs = plt.plot(solution[0], solution[1], 'x', label='solution')

        if FOR_LOOPING: # call solver in for loop to get all iterates
            cm = plt.cm.get_cmap('RdYlBu')
            axs = plt.scatter(iterates[:iter+1,0], iterates[:iter+1,1], c=range(iter+1), s=35, cmap=cm, label='iterates')
            plt.colorbar(axs)

        ts = np.linspace(0,2*np.pi,100)
        plt.plot(1 * np.cos(ts)+0,1 * np.sin(ts)-0, 'r')
        plt.axis('square')
        plt.legend()
        plt.title(f"Marathos problem with N = {N}, x formulation, SOC {globalization_use_SOC}")
        plt.show()

    print(f"\n\n----------------------\n")

if __name__ == '__main__':
    main()
