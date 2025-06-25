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

from acados_template import AcadosOcp, AcadosOcpSolver, ACADOS_INFTY, AcadosOcpIterate
import numpy as np
import scipy.linalg
from linear_mass_model import export_linear_mass_model


def create_solver_opts(N=4, Tf=2, nlp_solver_type = 'SQP_WITH_FEASIBLE_QP', allow_switching_modes=True, globalization= 'FUNNEL_L1PEN_LINESEARCH'):

    solver_options = AcadosOcp().solver_options

    # set options
    solver_options.N_horizon = N
    solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM'
    solver_options.qp_tol = 1e-9
    solver_options.qp_solver_ric_alg = 1
    solver_options.qp_solver_mu0 = 1e4
    solver_options.qp_solver_warm_start = 1
    solver_options.qp_solver_iter_max = 400
    solver_options.hessian_approx = 'GAUSS_NEWTON'
    solver_options.integrator_type = 'ERK'
    solver_options.nlp_solver_type = nlp_solver_type
    solver_options.globalization = globalization
    solver_options.globalization_full_step_dual = True
    solver_options.print_level = 1
    solver_options.nlp_solver_max_iter = 30
    solver_options.use_constraint_hessian_in_feas_qp = False
    solver_options.nlp_solver_ext_qp_res = 0

    if not allow_switching_modes:
        solver_options.search_direction_mode = 'BYRD_OMOJOKUN'
        solver_options.allow_direction_mode_switch_to_nominal = False

    # set prediction horizon
    solver_options.tf = Tf

    return solver_options

def create_solver(solver_name: str, soften_obstacle: bool, soften_terminal: bool,
                  soften_controls: bool, nlp_solver_type: str = 'SQP_WITH_FEASIBLE_QP',
                  globalization= 'FUNNEL_L1PEN_LINESEARCH',
                  allow_switching_modes: bool = True,
                  use_qp_scaling: bool = False):

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
    Fmax = 2
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
        ocp.cost.zl_e = 42 * 1e3 * np.ones(nx)
        ocp.cost.zu_e = 42 * 1e3 * np.ones(nx)
        ocp.cost.Zl_e = 0 * np.ones(nx)
        ocp.cost.Zu_e = 0 * np.ones(nx)

    # add cost for slacks
    if soften_controls:
        ocp.cost.zl = 1e1 * np.ones(nu)
        ocp.cost.zu = 1e1 * np.ones(nu)
        ocp.cost.Zl = 1e1 * np.ones(nu)
        ocp.cost.Zu = 1e1 * np.ones(nu)

    # add obstacle
    obs_rad = 1.0
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
    ocp.solver_options = create_solver_opts(N, Tf, nlp_solver_type, allow_switching_modes, globalization)
    if use_qp_scaling:
        ocp.solver_options.qpscaling_scale_constraints = "INF_NORM"
        ocp.solver_options.qpscaling_scale_objective = "OBJECTIVE_GERSHGORIN"

    # create ocp solver
    ocp_solver = AcadosOcpSolver(ocp, json_file=f'{model.name}_{solver_name}_ocp.json', verbose=False)

    # # initialize
    for i in range(N+1):
        ocp_solver.set(i, "x", (N+1-i)/(N+1) * x0 + i/(N+1) * x_goal)

    return ocp, ocp_solver



def check_qp_scaling(ocp_solver: AcadosOcpSolver):
    if ocp_solver.acados_ocp.solver_options.qpscaling_scale_constraints == "NO_CONSTRAINT_SCALING":
        try:
            constraint_scaling = ocp_solver.get_qp_scaling_constraints(0)
        except Exception as e:
            print(f"constraint scaling not done as expected.")
            return

    for i in range(ocp_solver.N+1):
        constraint_scaling = ocp_solver.get_qp_scaling_constraints(i)
        print(f"Constraint scaling at stage {i}: {constraint_scaling}")

    if ocp_solver.acados_ocp.solver_options.qpscaling_scale_objective != "NO_OBJECTIVE_SCALING":
        objective_scaling = ocp_solver.get_qp_scaling_objective()
        print(f"Objective scaling: {objective_scaling}")


def call_solver(ocp: AcadosOcp, ocp_solver: AcadosOcpSolver, soften_obstacle: bool,
                  soften_terminal: bool, soften_controls: bool, plot: bool) -> AcadosOcpIterate:
    # solve
    status = ocp_solver.solve()
    ocp_solver.print_statistics()

    sqp_iter = ocp_solver.get_stats('sqp_iter')
    if status != 0:
        raise RuntimeError(f"acados returned status {status} after {sqp_iter} SQP iterations.")
        # print(f'acados returned status {status}.')

    # print summary
    print(f"cost function value = {ocp_solver.get_cost()} after {sqp_iter} SQP iterations")
    print(f"solved sqp_wfqp problem with settings soften_obstacle = {soften_obstacle},soften_terminal = {soften_terminal}, SOFTEN_CONTROL = {soften_controls}")

def check_residual_solutions(stat1: np.ndarray, stat2: np.ndarray):
    n_rows1 = len(stat1[0])
    n_rows2 = len(stat2[0])

    assert n_rows1 == n_rows2, f"Both solvers should take the same number of iterations!, got {n_rows1} for solver 1, and {n_rows2} for solver 2"

    for jj in range(n_rows1):
        # res_stat
        assert np.allclose(stat1[1][jj], stat2[1][jj]), f"res_stat differs in iter {jj}"
        # res_eq
        assert np.allclose(stat1[2][jj], stat2[2][jj]), f"res_eq differs in iter {jj}"
        # res_ineq
        assert np.allclose(stat1[3][jj], stat2[3][jj]), f"res_ineq differs in iter {jj}"
        # res_comp
        assert np.allclose(stat1[4][jj], stat2[4][jj]), f"res_comp differs in iter {jj}"

def test_qp_scaling(nlp_solver_type = 'SQP', globalization = 'FUNNEL_L1PEN_LINESEARCH'):
    print(f"\n\nTesting solver={nlp_solver_type} with globalization={globalization}")
    # SETTINGS:
    soften_controls = True
    soften_obstacle = True
    soften_terminal = True

    # reference
    ocp_1, ocp_solver_1 = create_solver("1", soften_obstacle, soften_terminal, soften_controls, nlp_solver_type=nlp_solver_type, globalization=globalization, allow_switching_modes=False, use_qp_scaling=False)
    sol_1 = call_solver(ocp_1, ocp_solver_1, soften_obstacle, soften_terminal, soften_controls, plot=False)
    check_qp_scaling(ocp_solver_1)
    sol_1 = ocp_solver_1.store_iterate_to_obj()
    stats_1 = ocp_solver_1.get_stats("statistics")

    # test QP scaling
    ocp_2, ocp_solver_2 = create_solver("2", soften_obstacle, soften_terminal, soften_controls, nlp_solver_type=nlp_solver_type, allow_switching_modes=False, use_qp_scaling=True)
    sol_2 = call_solver(ocp_2, ocp_solver_2, soften_obstacle, soften_terminal, soften_controls, plot=False)
    check_qp_scaling(ocp_solver_2)
    sol_2 = ocp_solver_2.store_iterate_to_obj()
    ocp_solver_2.get_from_qp_in(1, "idxs_rev")
    stats_2 = ocp_solver_2.get_stats("statistics")

    check_residual_solutions(stats_1, stats_2)

    # check solutions
    for field in ["x_traj", "u_traj", "sl_traj", "su_traj", "lam_traj", "pi_traj"]:
        v1 = getattr(sol_1, field)
        v2 = getattr(sol_2, field)
        for i in range(len(v1)):
            if not np.allclose(v1[i], v2[i], atol=1e-6):
                print(f"Field {field} differs at index {i}: max diff = {np.max(np.abs(v1[i] - v2[i]))}")
                print(f"got difference {v1[i] - v2[i]}")
            else:
                pass
                # print(f"Field {field} is the same at index {i}.")

    # equivalent check
    if sol_1.allclose(sol_2, atol=1e-6):
        print("Both solvers have the same solution.")
    else:
        raise ValueError("Solutions of solvers differ!")

    print("\n\n---------------------------------------------------------")

if __name__ == '__main__':
    test_qp_scaling(nlp_solver_type = 'SQP', globalization = 'FUNNEL_L1PEN_LINESEARCH')
    test_qp_scaling(nlp_solver_type = 'SQP', globalization= 'MERIT_BACKTRACKING')
    test_qp_scaling(nlp_solver_type = 'SQP_WITH_FEASIBLE_QP', globalization = 'FUNNEL_L1PEN_LINESEARCH')
    test_qp_scaling(nlp_solver_type = 'SQP_WITH_FEASIBLE_QP', globalization= 'MERIT_BACKTRACKING')

