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
import scipy.linalg
from linear_mass_model import *
from itertools import product

# an OCP to test Maratos effect an second order correction

def main():

    # params = {'use_merit_fun_only': [True, False],
    params = {'use_merit_fun_only': [False],
              'initial_obj_multiplier': [1e0, 1e-1]}

    # # SETTINGS:
    SOFTEN_CONTROLS = True
    SOFTEN_OBSTACLE = False
    SOFTEN_TERMINAL = True
    PLOT = True
    solve_maratos_ocp(SOFTEN_OBSTACLE, SOFTEN_TERMINAL, SOFTEN_CONTROLS, PLOT, {'use_merit_fun_only': False,
                                                                                'initial_obj_multiplier':1e0})

    # SOFTEN_CONTROLS = False
    # SOFTEN_OBSTACLE = False
    # SOFTEN_TERMINAL = False
    # PLOT = True

    # keys, values = zip(*params.items())
    # for combination in product(*values):
    #     setting = dict(zip(keys, combination))
    #     solve_maratos_ocp(SOFTEN_OBSTACLE, SOFTEN_TERMINAL, SOFTEN_CONTROLS, PLOT, setting)

def feasible_qp_dims_test(SOFTEN_OBSTACLE, SOFTEN_TERMINAL, SOFTEN_CONTROLS, N, ocp_solver: AcadosOcpSolver):
    """
    The dynamics has four state variables and two control variables
    """
    dims = ocp_solver.acados_ocp.dims

    for i in range(N+1):
        idxs = ocp_solver.get_from_qp_in(i, "relaxed_idxs")
        idxb = ocp_solver.get_from_qp_in(i, "relaxed_idxb")

        # Initial stage
        if i == 0:
            assert len(idxb) == dims.nbu + dims.nbx_0, f"We should have {dims.nbu+dims.nbx} indices for bounds on x and u, but got {len(idxb)}"

            if not SOFTEN_CONTROLS:
                assert len(idxs) == 0, f"i=0, NOT SOFTEN_CONTROLS: The initial condition should have 0 slacks, got {len(idxs)}!"
            else:
                assert len(idxs) == dims.nbu, f"i=0, SOFTEN_CONTROLS: The initial condition should have {dims.nbu} slacks, got {len(idxs)}!"

        if i > 0 and i < N:
            assert len(idxb) == dims.nbu, f"We should have {dims.nbu} indices for bounds on u, but got {len(idxb)}"

            if not SOFTEN_CONTROLS:
                assert len(idxs) == dims.nh, f"i=0, NOT SOFTEN_CONTROLS: The initial condition should have {dims.nh} slacks, got {len(idxs)}!"
            else:
                assert len(idxs) == dims.nh + dims.nbu, f"i=0: SOFTEN_CONTROLS: The initial condition should have {dims.nh + dims.nbu} slacks, got {len(idxs)}!"

        # if not SOFTEN_CONTROLS and not SOFTEN_OBSTACLE and SOFTEN_TERMINAL:
        if i == N:
            # We slack the obstacle constraint and the terminal constraints
            assert len(idxs) == dims.nh_e + dims.nbx_e, f"i=N+1: Everything should be slacked, but got only {len(idxs)} slacks"

def feasible_qp_index_test(SOFTEN_OBSTACLE, SOFTEN_TERMINAL, SOFTEN_CONTROLS, N, ocp_solver: AcadosOcpSolver):
    """
    The dynamics has four state variables and two control variables
    """
    dims = ocp_solver.acados_ocp.dims

    for i in range(N+1):
        idxs = ocp_solver.get_from_qp_in(i, "relaxed_idxs").squeeze()
        idxb = ocp_solver.get_from_qp_in(i, "relaxed_idxb").squeeze()

        # Initial stage
        if i == 0:
            assert np.allclose(idxb, np.arange(dims.nbx_0 + dims.nbu)) , f"We should have {dims.nbx} bounds on x and u, but got {len(idxb)}"

            if not SOFTEN_CONTROLS:
                assert np.allclose(idxs,np.arange(0)), f"i=0, NOT SOFTEN_CONTROLS: The initial condition should have 0 slacks, got {len(idxs)}!"
            else:
                assert np.allclose(idxs, np.arange(dims.nbu)), f"i=0, SOFTEN_CONTROLS: The initial stage should have slack indices {np.arange(dims.nbu)} slacks, got {idxs})!"

        if i > 0 and i < N:
            assert np.allclose(idxb, np.arange(dims.nbu)), f"We should have {dims.nbu} indices for bounds on u, but got {len(idxb)}"

            if not SOFTEN_CONTROLS:
                assert np.allclose(idxs, np.arange(dims.nbx + dims.nbu, dims.nbx + dims.nbu + dims.nh)), f"i=0, NOT SOFTEN_CONTROLS: The initial condition should have {dims.nh} slacks, got {len(idxs)}!"
            else:
                assert np.allclose(idxs, np.arange(dims.nbx + dims.nbu + dims.nh)), f"i=0: SOFTEN_CONTROLS: The initial condition should have {dims.nh + dims.nbu} slacks, got {len(idxs)}!"

        # if not SOFTEN_CONTROLS and not SOFTEN_OBSTACLE and SOFTEN_TERMINAL:
        if i == N:
            # We slack the obstacle constraint and the terminal constraints
            assert np.allclose(idxs, np.arange(dims.nh_e + dims.nbx_e)), f"i=N+1: Everything should be slacked"

def create_solver_opts(setting: dict, N=4, Tf=2):

    use_merit_fun_only = setting['use_merit_fun_only']
    initial_obj_multiplier = setting['initial_obj_multiplier']

    solver_options = AcadosOcp().solver_options

    # set options
    solver_options.N_horizon = N
    solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM'
    # qp_tol = 5e-7
    # solver_options.qp_solver_tol_stat = qp_tol
    # solver_options.qp_solver_tol_eq = qp_tol
    # solver_options.qp_solver_tol_ineq = qp_tol
    # solver_options.qp_solver_tol_comp = qp_tol
    solver_options.qp_solver_ric_alg = 1
    solver_options.qp_solver_mu0 = 1e4
    solver_options.qp_solver_warm_start = 1
    solver_options.qp_solver_iter_max = 400
    # solver_options.hessian_approx = 'EXACT'
    # solver_options.regularize_method = 'MIRROR'
    solver_options.hessian_approx = 'GAUSS_NEWTON'
    # solver_options.levenberg_marquardt = 1e-2
    solver_options.integrator_type = 'ERK'
    solver_options.nlp_solver_type = 'SQP_WITH_FEASIBLE_QP'
    solver_options.globalization = 'FUNNEL_L1PEN_LINESEARCH'
    solver_options.globalization_full_step_dual = True
    solver_options.print_level = 1
    solver_options.nlp_solver_max_iter = 10
    solver_options.search_direction_mode = "BYRD_OMOJOKUN"

    solver_options.globalization_funnel_use_merit_fun_only = use_merit_fun_only
    solver_options.initial_objective_multiplier = initial_obj_multiplier
    if False:
        solver_options.use_exact_hessian_in_feas_qp = True
    else:
        solver_options.use_exact_hessian_in_feas_qp = False

    # set prediction horizon
    solver_options.tf = Tf

    return solver_options

def solve_maratos_ocp(SOFTEN_OBSTACLE, SOFTEN_TERMINAL, SOFTEN_CONTROLS, PLOT, setting):

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
    N = 4
    shooting_nodes = np.linspace(0, Tf, N+1)

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

    # Slack the controls
    if SOFTEN_CONTROLS:
        ocp.constraints.idxsbu = np.array(range(nu))

    x0 = np.array([1e-1, 1.1, 0, 0])
    ocp.constraints.x0 = x0

    # terminal constraint
    x_goal = np.array([0, -1.1, 0, 0])
    ocp.constraints.idxbx_e = np.array(range(nx))
    ocp.constraints.lbx_e = x_goal
    ocp.constraints.ubx_e = x_goal

    if SOFTEN_TERMINAL:
        ocp.constraints.idxsbx_e = np.array(range(nx))
        ocp.cost.zl_e = 42 * 1e5 * np.ones(nx)
        ocp.cost.zu_e = 42 * 1e5 * np.ones(nx)
        ocp.cost.Zl_e = 0 * np.ones(nx)
        ocp.cost.Zu_e = 0 * np.ones(nx)

    # add cost for slacks
    if SOFTEN_CONTROLS:
        ocp.cost.zl = 1e1 * np.ones(nu)
        ocp.cost.zu = 1e1 * np.ones(nu)
        ocp.cost.Zl = 1e1 * np.ones(nu)
        ocp.cost.Zu = 1e1 * np.ones(nu)

    # add obstacle
    obs_rad = 1.0; obs_x = 0.0; obs_y = 0.0
    circle = (obs_x, obs_y, obs_rad)
    ocp.constraints.lh = -np.array([ACADOS_INFTY])
    ocp.constraints.uh = -np.array([obs_rad**2])
    x_square = model.x[0] ** 2 + model.x[1] ** 2
    ocp.model.con_h_expr = -x_square
    # copy for terminal
    ocp.constraints.uh_e = ocp.constraints.uh
    ocp.constraints.lh_e = ocp.constraints.lh
    ocp.model.con_h_expr_e = ocp.model.con_h_expr

    # # soften
    if SOFTEN_OBSTACLE:
        ocp.constraints.idxsh = np.array([0])
        ocp.constraints.idxsh_e = np.array([0])
        Zh = 1e6 * np.ones(1)
        zh = 1e4 * np.ones(1)
        ocp.cost.zl = zh
        ocp.cost.zu = zh
        ocp.cost.Zl = Zh
        ocp.cost.Zu = Zh
        ocp.cost.zl_e = np.concatenate((ocp.cost.zl_e, zh))
        ocp.cost.zu_e = np.concatenate((ocp.cost.zu_e, zh))
        ocp.cost.Zl_e = np.concatenate((ocp.cost.Zl_e, Zh))
        ocp.cost.Zu_e = np.concatenate((ocp.cost.Zu_e, Zh))

    # load options
    ocp.solver_options = create_solver_opts(setting, N, Tf)
    # create ocp solver
    ocp_solver = AcadosOcpSolver(ocp, json_file=f'{model.name}_ocp.json', verbose=False)

    # # initialize
    for i in range(N+1):
        ocp_solver.set(i, "x", (N+1-i)/(N+1) * x0 + i/(N+1) * x_goal)

    # solve
    status = ocp_solver.solve()

    sqp_iter = ocp_solver.get_stats('sqp_iter')
    print(f'acados returned status {status}.')

    if ocp.solver_options.nlp_solver_type == 'SQP_WITH_FEASIBLE_QP':
        feasible_qp_dims_test(SOFTEN_OBSTACLE, SOFTEN_TERMINAL, SOFTEN_CONTROLS, N, ocp_solver)
        feasible_qp_index_test(SOFTEN_OBSTACLE, SOFTEN_TERMINAL, SOFTEN_CONTROLS, N, ocp_solver)

    # get solution
    simX = np.array([ocp_solver.get(i,"x") for i in range(N+1)])
    simU = np.array([ocp_solver.get(i,"u") for i in range(N)])
    pi_multiplier = [ocp_solver.get(i, "pi") for i in range(N)]

    # We should put the optimal solution here ....

    # print summary
    print(f"cost function value = {ocp_solver.get_cost()} after {sqp_iter} SQP iterations")
    print(f"solved sqp_wfqp problem with settings SOFTEN_OBSTACLE = {SOFTEN_OBSTACLE},SOFTEN_TERMINAL = {SOFTEN_TERMINAL}, SOFTEN_CONTROL = {SOFTEN_CONTROLS}, use_merit_fun_only = {setting['use_merit_fun_only']}, funnel_init_penalty_parameter = {setting['initial_obj_multiplier']}")

    if PLOT:
        plot_linear_mass_system_X_state_space(simX, circle=circle, x_goal=x_goal)

    if ocp.solver_options.globalization_funnel_use_merit_fun_only:
        assert status == 0, "Merit function should always converge!"
    # if not ocp.solver_options.globalization_funnel_use_merit_fun_only and\
    #     ocp.solver_options.initial_objective_multiplier == 1e0\
    #     and not SOFTEN_CONTROLS:
    #     assert status in [3,8], "Funnel should converge to an infeasible point at the moment for unsoftened problem and initial penalty parameter 1e0, got status {status}!"
    elif not ocp.solver_options.globalization_funnel_use_merit_fun_only and\
        ocp.solver_options.initial_objective_multiplier != 1e0:
        assert status == 0, "Funnel should find solution!"

    print(f"\n\n----------------------\n")

if __name__ == '__main__':
    main()
