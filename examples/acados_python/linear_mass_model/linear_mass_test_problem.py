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
import scipy.linalg
from linear_mass_model import *
from itertools import product

## SETTINGS:
OBSTACLE = True
SOFTEN_OBSTACLE = False
SOFTEN_TERMINAL = True
INITIALIZE = True
PLOT = False
OBSTACLE_POWER = 2

# an OCP to test Maratos effect an second order correction

def main():
    # run test cases

    # all setting
    params = {'globalization': ['FIXED_STEP', 'MERIT_BACKTRACKING', 'FUNNEL_L1PEN_LINESEARCH'], # MERIT_BACKTRACKING, FIXED_STEP
              'globalization_line_search_use_sufficient_descent' : [0, 1],
              'qp_solver' : ['FULL_CONDENSING_HPIPM', 'PARTIAL_CONDENSING_HPIPM', 'FULL_CONDENSING_QPOASES'],
              'globalization_use_SOC' : [0, 1] }

    keys, values = zip(*params.items())
    for combination in product(*values):
        setting = dict(zip(keys, combination))
        if (setting['globalization'] == 'FIXED_STEP' or setting['globalization'] == 'FUNNEL_L1PEN_LINESEARCH') and \
          (setting['globalization_use_SOC'] or setting['globalization_line_search_use_sufficient_descent']):
            # skip some equivalent settings
            pass
        else:
            solve_maratos_ocp(setting)

    setting = {}
    setting['globalization'] = 'MERIT_BACKTRACKING'
    setting['globalization_use_SOC'] = 0
    setting['qp_solver'] = 'PARTIAL_CONDENSING_HPIPM'
    setting['globalization_line_search_use_sufficient_descent'] = 1
    print("Use deprecated options!")
    solve_maratos_ocp(setting, use_deprecated_options=True)


def solve_maratos_ocp(setting, use_deprecated_options=False):

    globalization = setting['globalization']
    globalization_line_search_use_sufficient_descent = setting['globalization_line_search_use_sufficient_descent']
    globalization_use_SOC = setting['globalization_use_SOC']
    qp_solver = setting['qp_solver']

    # create ocp object to formulate the OCP
    ocp = AcadosOcp()

    # set model
    model = export_linear_mass_model()
    ocp.model = model

    nx = model.x.rows()
    nu = model.u.rows()
    ny = nu

    # discretization
    Tf = 2
    N = 20
    shooting_nodes = np.linspace(0, Tf, N+1)
    ocp.solver_options.N_horizon = N

    # set cost
    Q = 2*np.diag([])
    R = 2*np.diag([1e1, 1e1])

    ocp.cost.W_e = Q
    ocp.cost.W = scipy.linalg.block_diag(Q, R)

    ocp.cost.cost_type = 'LINEAR_LS'
    ocp.cost.cost_type_e = 'LINEAR_LS'

    ocp.cost.Vx = np.zeros((ny, nx))

    Vu = np.eye((nu))
    ocp.cost.Vu = Vu
    ocp.cost.yref = np.zeros((ny, ))

    # set constraints
    Fmax = 5
    ocp.constraints.lbu = -Fmax * np.ones((nu,))
    ocp.constraints.ubu = +Fmax * np.ones((nu,))
    ocp.constraints.idxbu = np.array(range(nu))
    x0 = np.array([1e-1, 1.1, 0, 0])
    ocp.constraints.x0 = x0

    # terminal constraint
    x_goal = np.array([0, -1.1, 0, 0])
    ocp.constraints.idxbx_e = np.array(range(nx))
    ocp.constraints.lbx_e = x_goal
    ocp.constraints.ubx_e = x_goal

    if SOFTEN_TERMINAL:
        ocp.constraints.idxsbx_e = np.array(range(nx))
        ocp.cost.zl_e = 1e4 * np.ones(nx)
        ocp.cost.zu_e = 1e4 * np.ones(nx)
        ocp.cost.Zl_e = 1e6 * np.ones(nx)
        ocp.cost.Zu_e = 1e6 * np.ones(nx)

    # add obstacle
    if OBSTACLE:
        obs_rad = 1.0; obs_x = 0.0; obs_y = 0.0
        circle = (obs_x, obs_y, obs_rad)
        ocp.constraints.uh = np.array([100.0]) # doenst matter
        ocp.constraints.lh = np.array([obs_rad**2])
        x_square = model.x[0] ** OBSTACLE_POWER + model.x[1] ** OBSTACLE_POWER
        ocp.model.con_h_expr = x_square
        # copy for terminal
        ocp.constraints.uh_e = ocp.constraints.uh
        ocp.constraints.lh_e = ocp.constraints.lh
        ocp.model.con_h_expr_e = ocp.model.con_h_expr
    else:
        circle = None

    # soften
    if OBSTACLE and SOFTEN_OBSTACLE:
        ocp.constraints.idxsh = np.array([0])
        ocp.constraints.idxsh_e = np.array([0])
        Zh = 1e6 * np.ones(1)
        zh = 1e4 * np.ones(1)
        # initial: no obstacle constraint, no addtional slack
        ocp.cost.zl_0 = ocp.cost.zl
        ocp.cost.zu_0 = ocp.cost.zu
        ocp.cost.Zl_0 = ocp.cost.Zl
        ocp.cost.Zu_0 = ocp.cost.Zu
        # path & terminal: slacked obstacle constraint
        ocp.cost.zl = np.concatenate((ocp.cost.zl, zh))
        ocp.cost.zu = np.concatenate((ocp.cost.zu, zh))
        ocp.cost.Zl = np.concatenate((ocp.cost.Zl, Zh))
        ocp.cost.Zu = np.concatenate((ocp.cost.Zu, Zh))
        ocp.cost.zl_e = np.concatenate((ocp.cost.zl_e, zh))
        ocp.cost.zu_e = np.concatenate((ocp.cost.zu_e, zh))
        ocp.cost.Zl_e = np.concatenate((ocp.cost.Zl_e, Zh))
        ocp.cost.Zu_e = np.concatenate((ocp.cost.Zu_e, Zh))

    # set options
    ocp.solver_options.qp_solver = qp_solver # FULL_CONDENSING_QPOASES
    # PARTIAL_CONDENSING_HPIPM, FULL_CONDENSING_QPOASES, FULL_CONDENSING_HPIPM,
    # PARTIAL_CONDENSING_QPDUNES, PARTIAL_CONDENSING_OSQP
    ocp.solver_options.hessian_approx = 'GAUSS_NEWTON'
    ocp.solver_options.integrator_type = 'ERK'
    # ocp.solver_options.print_level = 1
    ocp.solver_options.nlp_solver_type = 'SQP' # SQP_RTI, SQP
    ocp.solver_options.globalization = globalization
    if use_deprecated_options:
        ocp.solver_options.alpha_min = 0.01
        ocp.solver_options.globalization_use_SOC = setting['globalization_use_SOC']
        ocp.solver_options.full_step_dual = 1
        ocp.solver_options.eps_sufficient_descent = 1e-4
        ocp.solver_options.line_search_use_sufficient_descent = setting['globalization_line_search_use_sufficient_descent']
    else:
        ocp.solver_options.globalization_alpha_min = 0.01
    # ocp.solver_options.levenberg_marquardt = 1e-2
    ocp.solver_options.qp_solver_cond_N = 0
    ocp.solver_options.print_level = 1
    ocp.solver_options.nlp_solver_max_iter = 200
    ocp.solver_options.qp_solver_iter_max = 400
    # NOTE: this is needed for PARTIAL_CONDENSING_HPIPM to get expected behavior
    qp_tol = 5e-7
    ocp.solver_options.qp_solver_tol_stat = qp_tol
    ocp.solver_options.qp_solver_tol_eq = qp_tol
    ocp.solver_options.qp_solver_tol_ineq = qp_tol
    ocp.solver_options.qp_solver_tol_comp = qp_tol
    ocp.solver_options.qp_solver_ric_alg = 1

    # set prediction horizon
    ocp.solver_options.tf = Tf

    ocp_solver = AcadosOcpSolver(ocp, json_file=f'{model.name}_ocp.json', verbose=False)

    if globalization == "FUNNEL_L1PEN_LINESEARCH":
        # Test the options setters
        ocp_solver.options_set('globalization_funnel_init_increase_factor', 15.0)
        ocp_solver.options_set('globalization_funnel_init_upper_bound', 1.0)
        ocp_solver.options_set('globalization_funnel_sufficient_decrease_factor', 0.9)
        ocp_solver.options_set('globalization_funnel_kappa', 0.9)
        ocp_solver.options_set('globalization_funnel_fraction_switching_condition', 1e-3)
        ocp_solver.options_set('globalization_funnel_initial_penalty_parameter', 1.0)
    if not use_deprecated_options:
        ocp_solver.options_set('globalization_line_search_use_sufficient_descent', globalization_line_search_use_sufficient_descent)
        ocp_solver.options_set('globalization_use_SOC', globalization_use_SOC)
        ocp_solver.options_set('globalization_full_step_dual', 1)

    if INITIALIZE:# initialize solver
        # [ocp_solver.set(i, "x", x0 + (i/N) * (x_goal-x0)) for i in range(N+1)]
        [ocp_solver.set(i, "x", x0) for i in range(N+1)]
        # [ocp_solver.set(i, "u", 2*(np.random.rand(2) - 0.5)) for i in range(N)]

    # solve
    status = ocp_solver.solve()
    ocp_solver.print_statistics() # encapsulates: stat = ocp_solver.get_stats("statistics")
    sqp_iter = ocp_solver.get_stats('sqp_iter')
    print(f'acados returned status {status}.')

    # ocp_solver.store_iterate(f'it{ocp.solver_options.nlp_solver_max_iter}_{model.name}.json')

    # get solution
    simX = np.array([ocp_solver.get(i,"x") for i in range(N+1)])
    simU = np.array([ocp_solver.get(i,"u") for i in range(N)])
    pi_multiplier = [ocp_solver.get(i, "pi") for i in range(N)]

    # print summary
    print(f"solved Maratos test problem with settings {setting}")
    print(f"cost function value = {ocp_solver.get_cost()} after {sqp_iter} SQP iterations")
    # print(f"alphas: {alphas[:iter]}")
    # print(f"total number of QP iterations: {sum(qp_iters[:iter])}")
    # max_infeasibility = np.max(residuals[1:3])
    # print(f"max infeasibility: {max_infeasibility}")

    # checks
    if status != 0:
        raise Exception(f"acados solver returned status {status} != 0.")
    if globalization == "FIXED_STEP":
        if sqp_iter != 18:
            raise Exception(f"acados solver took {sqp_iter} iterations, expected 18.")
    elif globalization == "MERIT_BACKTRACKING" and not use_deprecated_options:
        if sqp_iter not in range(17, 23):
            raise Exception(f"acados solver took {sqp_iter} iterations, expected range(17, 23).")

    if PLOT:
        plot_linear_mass_system_X_state_space(simX, circle=circle, x_goal=x_goal)
        plot_linear_mass_system_U(shooting_nodes, simU)
    print(f"\n\n----------------------\n")

if __name__ == '__main__':
    main()
