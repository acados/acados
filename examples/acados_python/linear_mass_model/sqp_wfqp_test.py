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

from acados_template import AcadosOcp, AcadosOcpSolver, ACADOS_INFTY
import numpy as np
import scipy.linalg
from linear_mass_model import export_linear_mass_model, plot_linear_mass_system_X_state_space

# an OCP to test the behavior of the SQP_WITH_FEASIBLE_QP functionalities


def feasible_qp_dims_test(soften_obstacle, soften_terminal, soften_controls, N, ocp_solver: AcadosOcpSolver):
    """
    The dynamics has four state variables and two control variables
    """
    dims = ocp_solver.acados_ocp.dims

    for i in range(N+1):
        idxs = ocp_solver.get_from_qp_in(i, "relaxed_idxs")
        idxb = ocp_solver.get_from_qp_in(i, "relaxed_idxb")

        # Initial stage
        if i == 0:
            assert len(idxb) == dims.nbu + dims.nbx_0, f"We should have {dims.nbu+dims.nbx} indices for bounds on x and u at stage {i}, but got {len(idxb)}"

            if not soften_controls:
                assert len(idxs) == 0, f"i=0, NOT soften_controls: The initial condition should have 0 slacks, got {len(idxs)}!"
            else:
                assert len(idxs) == dims.nbu, f"i=0, soften_controls: The initial condition should have {dims.nbu} slacks, got {len(idxs)}!"

        if i > 0 and i < N:
            assert len(idxb) == dims.nbu, f"We should have {dims.nbu} indices for bounds on u, but got {len(idxb)}"

            if not soften_controls:
                assert len(idxs) == dims.nh, f"i=0, NOT soften_controls: The initial condition should have {dims.nh} slacks, got {len(idxs)}!"
            else:
                assert len(idxs) == dims.nh + dims.nbu, f"i=0: soften_controls: The initial condition should have {dims.nh + dims.nbu} slacks, got {len(idxs)}!"

        # if not soften_controls and not soften_obstacle and soften_terminal:
        if i == N:
            # We slack the obstacle constraint and the terminal constraints
            assert len(idxs) == dims.nh_e + dims.nbx_e, f"i=N+1: Everything should be slacked, but got only {len(idxs)} slacks"

def feasible_qp_index_test(soften_obstacle, soften_terminal, soften_controls, N, ocp_solver: AcadosOcpSolver):
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

            if not soften_controls:
                assert np.allclose(idxs,np.arange(0)), f"i=0, NOT soften_controls: The initial condition should have 0 slacks, got {len(idxs)}!"
            else:
                assert np.allclose(idxs, np.arange(dims.nbu)), f"i=0, soften_controls: The initial stage should have slack indices {np.arange(dims.nbu)} slacks, got {idxs})!"

        if i > 0 and i < N:
            assert np.allclose(idxb, np.arange(dims.nbu)), f"We should have {dims.nbu} indices for bounds on u, but got {len(idxb)}"

            if not soften_controls:
                assert np.allclose(idxs, np.arange(dims.nbx + dims.nbu, dims.nbx + dims.nbu + dims.nh)), f"i=0, NOT soften_controls: The initial condition should have {dims.nh} slacks, got {len(idxs)}!"
            else:
                assert np.allclose(idxs, np.arange(dims.nbx + dims.nbu + dims.nh)), f"i=0: soften_controls: The initial condition should have {dims.nh + dims.nbu} slacks, got {len(idxs)}!"

        # TODO: rework here!
        # if not soften_controls and not soften_obstacle and soften_terminal:
        if i == N:
            # We slack the obstacle constraint and the terminal constraints
            assert np.allclose(idxs, np.arange(dims.nh_e + dims.nbx_e)), f"i=N+1: Everything should be slacked"

def create_solver_opts(N=4, Tf=2, nlp_solver_type = 'SQP_WITH_FEASIBLE_QP', allow_switching_modes=True):

    solver_options = AcadosOcp().solver_options

    # set options
    solver_options.N_horizon = N
    solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM'
    qp_tol = 5e-7
    solver_options.qp_tol = qp_tol
    solver_options.qp_solver_ric_alg = 1
    solver_options.qp_solver_mu0 = 1e4
    solver_options.qp_solver_warm_start = 1
    solver_options.qp_solver_iter_max = 400
    solver_options.hessian_approx = 'GAUSS_NEWTON'
    solver_options.integrator_type = 'ERK'
    solver_options.nlp_solver_type = nlp_solver_type
    solver_options.globalization = 'FUNNEL_L1PEN_LINESEARCH'
    solver_options.globalization_full_step_dual = True
    solver_options.print_level = 1
    solver_options.nlp_solver_max_iter = 20
    solver_options.use_constraint_hessian_in_feas_qp = False

    if not allow_switching_modes:
        solver_options.search_direction_mode = 'BYRD_OMOJOKUN'
        solver_options.allow_direction_mode_switch_to_nominal = False

    # set prediction horizon
    solver_options.tf = Tf

    return solver_options

def create_solver(solver_name: str, soften_obstacle: bool, soften_terminal: bool,
                  soften_controls: bool, nlp_solver_type: str = 'SQP_WITH_FEASIBLE_QP',
                  allow_switching_modes: bool = True):

    # create ocp object to formulate the OCP
    ocp = AcadosOcp()

    # set model
    model = export_linear_mass_model()
    ocp.model = model
    ocp.model.name += solver_name

    nx = model.x.rows()
    nu = model.u.rows()
    ny = nu

    # discretization
    Tf = 2
    N = 4

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
    if soften_controls:
        ocp.constraints.idxsbu = np.array(range(nu))

    x0 = np.array([1e-1, 1.1, 0, 0])
    ocp.constraints.x0 = x0

    # terminal constraint
    x_goal = np.array([0, -1.1, 0, 0])
    ocp.constraints.idxbx_e = np.array(range(nx))
    ocp.constraints.lbx_e = x_goal
    ocp.constraints.ubx_e = x_goal

    if soften_terminal:
        ocp.constraints.idxsbx_e = np.array(range(nx))
        ocp.cost.zl_e = 42 * 1e5 * np.ones(nx)
        ocp.cost.zu_e = 42 * 1e5 * np.ones(nx)
        ocp.cost.Zl_e = 0 * np.ones(nx)
        ocp.cost.Zu_e = 0 * np.ones(nx)

    # add cost for slacks
    if soften_controls:
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
    if soften_obstacle:
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

    # load options
    ocp.solver_options = create_solver_opts(N, Tf, nlp_solver_type, allow_switching_modes)

    # create ocp solver
    ocp_solver = AcadosOcpSolver(ocp, json_file=f'{model.name}_{solver_name}_ocp.json', verbose=False)

    # # initialize
    for i in range(N+1):
        ocp_solver.set(i, "x", (N+1-i)/(N+1) * x0 + i/(N+1) * x_goal)

    return ocp, ocp_solver

def standard_test(ocp: AcadosOcp, ocp_solver: AcadosOcpSolver, soften_obstacle: bool,
                  soften_terminal: bool, soften_controls: bool, plot: bool):
    # solve
    status = ocp_solver.solve()

    N = ocp.solver_options.N_horizon

    sqp_iter = ocp_solver.get_stats('sqp_iter')
    print(f'acados returned status {status}.')

    if ocp.solver_options.nlp_solver_type == 'SQP_WITH_FEASIBLE_QP':
        feasible_qp_dims_test(soften_obstacle, soften_terminal, soften_controls, N, ocp_solver)
        feasible_qp_index_test(soften_obstacle, soften_terminal, soften_controls, N, ocp_solver)

    # get solution
    sol_X = np.array([ocp_solver.get(i,"x") for i in range(N+1)])

    # print summary
    print(f"cost function value = {ocp_solver.get_cost()} after {sqp_iter} SQP iterations")
    print(f"solved sqp_wfqp problem with settings soften_obstacle = {soften_obstacle},soften_terminal = {soften_terminal}, SOFTEN_CONTROL = {soften_controls}")

    if plot:
        obs_rad = 1.0; obs_x = 0.0; obs_y = 0.0
        circle = (obs_x, obs_y, obs_rad)
        x_goal = x_goal = np.array([0, -1.1, 0, 0])
        plot_linear_mass_system_X_state_space(sol_X, circle=circle, x_goal=x_goal)

    print(f"\n\n----------------------\n")

def test_same_behavior_sqp_and_sqp_wfqp():
    # # SETTINGS:
    soften_controls = True
    soften_obstacle = False
    soften_terminal = True

    # SQP solver
    _, ocp_solver1 = create_solver("v1", soften_obstacle, soften_terminal, soften_controls, nlp_solver_type='SQP_WITH_FEASIBLE_QP')
    status1 = ocp_solver1.solve()

    _, ocp_solver2 = create_solver("v2", soften_obstacle, soften_terminal, soften_controls, nlp_solver_type='SQP')
    status2 = ocp_solver2.solve()

    assert status1 == status2, "both solvers should converge"

    # check residuals
    res_solver1 = ocp_solver1.get_residuals()
    res_solver2 = ocp_solver2.get_residuals()
    assert np.array_equal(res_solver1, res_solver2), "both solvers should have identical residual stats"

    # check solutions
    sol_1 = ocp_solver1.store_iterate_to_flat_obj()
    sol_2 = ocp_solver2.store_iterate_to_flat_obj()
    if sol_1.allclose(sol_2):
        print("Both solvers have the same solution.")
    else:
        raise ValueError("Solutions of solvers differ!")

    print(f"\n\n----------------------\n")

def sqp_wfqp_test_same_matrices():
    # # SETTINGS:
    soften_controls = False
    soften_obstacle = False
    soften_terminal = False

    # SQP solver
    _, ocp_solver1 = create_solver("v1", soften_obstacle, soften_terminal, soften_controls, nlp_solver_type='SQP_WITH_FEASIBLE_QP', allow_switching_modes=False)
    _ = ocp_solver1.solve()

    qp = ocp_solver1.get_last_qp()
    relaxed_qp = ocp_solver1.get_last_relaxed_qp()

    dynamics_and_bu = ['b_', 'A_', 'B_', 'lbu_', 'ubu_']
    for prefix in dynamics_and_bu:
        for i in range(ocp_solver1.N):
            assert np.equal(qp[prefix+str(i)], relaxed_qp['relaxed_'+prefix+str(i)]).all(), f" matrices do not coincide for {prefix}{i},"

    constraints =  ['C_', 'D_']
    for prefix in constraints:
        for i in range(1, ocp_solver1.N):
            assert np.equal(qp[prefix+str(i)], relaxed_qp['relaxed_'+prefix+str(i)]).all(), f" matrices do not coincide for {prefix}{i},"

    idxb =  ['idxb_']
    for prefix in idxb:
        for i in range(0, ocp_solver1.N+1):
            assert np.equal(qp[prefix+str(i)], relaxed_qp['relaxed_'+prefix+str(i)]).all(), f" matrices do not coincide for {prefix}{i},"

    print(f"\n\n----------------------\n")


def main_test():
    # SETTINGS:
    soften_controls = True
    soften_obstacle = False
    soften_terminal = True
    plot = True
    ocp, ocp_solver1 = create_solver("1", soften_obstacle, soften_terminal, soften_controls)
    standard_test(ocp, ocp_solver1, soften_obstacle, soften_terminal, soften_controls, plot)

    soften_controls = False
    soften_obstacle = False
    soften_terminal = False
    plot = True
    ocp, ocp_solver2 = create_solver("2", soften_obstacle, soften_terminal, soften_controls)
    standard_test(ocp, ocp_solver2, soften_obstacle, soften_terminal, soften_controls, plot)


if __name__ == '__main__':
    main_test()
    test_same_behavior_sqp_and_sqp_wfqp()
    sqp_wfqp_test_same_matrices()
