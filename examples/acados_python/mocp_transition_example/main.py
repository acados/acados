# -*- coding: future_fstrings -*-
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

import numpy as np
import casadi as ca

from acados_template import AcadosModel, AcadosOcp, AcadosMultiphaseOcp, AcadosOcpSolver, casadi_length, is_empty, latexify_plot
import matplotlib.pyplot as plt

latexify_plot()

X0 = np.array([2.0, 0.0])
PENALTY_X = 1e0
T_HORIZON = 1.0
N_HORIZON = 20

L2_COST_V = 1e-1
L2_COST_P = 1e0
L2_COST_A = 1e-3

def get_double_integrator_model() -> AcadosModel:

    # set up states & controls
    p = ca.SX.sym('p')
    v = ca.SX.sym('v')
    x = ca.vertcat(p, v)
    nx = casadi_length(x)

    u = ca.SX.sym('u')

    # set up dynamics
    f_expl = ca.vertcat(v, u)
    xdot = ca.SX.sym('xdot', nx)
    f_impl = f_expl - xdot

    # set up model
    model = AcadosModel()
    model.name = 'double_integrator'
    model.x = x
    model.xdot = xdot
    model.u = u
    model.f_impl_expr = f_impl
    model.f_expl_expr = f_expl

    return model


def get_single_integrator_model() -> AcadosModel:

    # set up states & controls
    p = ca.SX.sym('p')
    x = ca.vertcat(p)
    nx = casadi_length(x)

    u = ca.SX.sym('v')

    # set up dynamics
    f_expl = ca.vertcat(u)
    xdot = ca.SX.sym('xdot', nx)
    f_impl = f_expl - xdot

    # set up model
    model = AcadosModel()
    model.name = 'single_integrator'
    model.x = x
    model.xdot = xdot
    model.u = u
    model.f_impl_expr = f_impl
    model.f_expl_expr = f_expl

    return model


def get_transition_model() -> AcadosModel:

    # set up states & controls
    p = ca.SX.sym('p')
    v = ca.SX.sym('v')
    x = ca.vertcat(p, v)

    # set up dynamics
    new_x = p

    # set up model
    model = AcadosModel()
    model.name = 'transition_model'
    model.x = x
    model.u = ca.SX.sym('u', 0, 0)
    model.disc_dyn_expr = new_x

    return model


def formulate_double_integrator_ocp() -> AcadosOcp:
    ocp = AcadosOcp()

    ocp.model = get_double_integrator_model()

    ocp.dims.N = N_HORIZON

    ocp.cost.cost_type = 'NONLINEAR_LS'
    ocp.cost.cost_type_e = 'NONLINEAR_LS'
    ocp.cost.W = np.diag([L2_COST_P, L2_COST_V, L2_COST_A])
    ocp.cost.W_e = np.diag([1e1, 1e1])
    ocp.cost.yref = np.array([0.0, 0.0, 0.0])
    ocp.cost.yref_e = np.array([0.0, 0.0])

    ocp.model.cost_y_expr = ca.vertcat(ocp.model.x, ocp.model.u)
    ocp.model.cost_y_expr_e = ocp.model.x

    u_max = 50.0
    ocp.constraints.lbu = np.array([-u_max])
    ocp.constraints.ubu = np.array([u_max])
    ocp.constraints.idxbu = np.array([0])

    ocp.constraints.x0 = X0

    return ocp


def formulate_single_integrator_ocp() -> AcadosOcp:
    ocp = AcadosOcp()

    ocp.model = get_single_integrator_model()

    ocp.dims.N = N_HORIZON

    ocp.cost.cost_type = 'NONLINEAR_LS'
    ocp.cost.cost_type_e = 'NONLINEAR_LS'
    ocp.cost.W = np.diag([L2_COST_P, L2_COST_V])
    ocp.cost.W_e = np.diag([1e1])
    ocp.cost.yref = np.array([0.0, 0.0])
    ocp.cost.yref_e = np.array([0.0])

    ocp.model.cost_y_expr = ca.vertcat(ocp.model.x, ocp.model.u)
    ocp.model.cost_y_expr_e = ocp.model.x

    u_max = 5.0
    ocp.constraints.lbu = np.array([-u_max])
    ocp.constraints.ubu = np.array([u_max])
    ocp.constraints.idxbu = np.array([0])

    return ocp


def main_multiphase_ocp():

    N_list = [10, 1, 15]
    ocp = AcadosMultiphaseOcp(N_list=N_list)

    phase_0 = formulate_double_integrator_ocp()
    ocp.set_phase(phase_0, 0)

    phase_1 = AcadosOcp()
    phase_1.model = get_transition_model()
    phase_1.cost.cost_type = 'NONLINEAR_LS'
    phase_1.model.cost_y_expr = phase_1.model.x
    # TODO: how to choose this transition cost
    phase_1.cost.W = np.diag([L2_COST_P, 1e-1 * L2_COST_V])
    phase_1.cost.yref = np.array([0., 0.])

    ocp.set_phase(phase_1, 1)

    phase_2 = formulate_single_integrator_ocp()
    ocp.set_phase(phase_2, 2)

    ocp.solver_options.nlp_solver_type = 'SQP'

    # the transition stage uses discrete dynamics!
    ocp.mocp_opts.integrator_type = ['IRK', 'DISCRETE', 'IRK']

    # the time step for the transition stage is set to 1 to correct for the scaling of the stage cost with the time step
    ocp.solver_options.tf = T_HORIZON + 1.0
    T_HORIZON_1 = 0.4 * T_HORIZON
    T_HORIZON_2 = T_HORIZON - T_HORIZON_1
    ocp.solver_options.time_steps = np.array(N_list[0] * [T_HORIZON_1/N_list[0]] + [1.0] + N_list[2] * [T_HORIZON_2/N_list[2]] )

    acados_ocp_solver = AcadosOcpSolver(ocp, verbose=False)
    acados_ocp_solver.solve_for_x0(X0)

    n_phases = len(N_list)

    x_traj_phases = n_phases*[None]
    u_traj_phases = n_phases*[None]
    t_grid_phases = n_phases*[None]

    for i_phase in range(n_phases):
        x_traj_phases[i_phase] = [acados_ocp_solver.get(i, 'x') for i in range(ocp.start_idx[i_phase], ocp.end_idx[i_phase]+1)]
        u_traj_phases[i_phase] = [acados_ocp_solver.get(i, 'u') for i in range(ocp.start_idx[i_phase], ocp.end_idx[i_phase])]
        t_grid_phases[i_phase] = ocp.solver_options.shooting_nodes[ocp.start_idx[i_phase]: ocp.end_idx[i_phase]+1]
        print(f"Phase {i_phase}:\nt grid \n {t_grid_phases[i_phase]} \nx traj\n {x_traj_phases[i_phase]} \nu traj {u_traj_phases[i_phase]}")
        print("-----------------------------------")

    # plot solution
    t_grid_2_plot = t_grid_phases[2] - 1.0

    fig, ax = plt.subplots(3, 1, sharex=True)

    p_traj_0 = [x[0] for x in x_traj_phases[0]]
    ax[0].plot(t_grid_phases[0], p_traj_0, color='C0', label='phase 1')

    p_traj_2 = [x[0] for x in x_traj_phases[2]]
    ax[0].plot(t_grid_2_plot, p_traj_2, color='C1', label='phase 2')

    v_traj_0 = [x[1] for x in x_traj_phases[0]]
    ax[1].plot(t_grid_phases[0], v_traj_0, color='C0')

    v_traj_2 = [u_traj_phases[2][0][0]] + [x[0] for x in u_traj_phases[2]]
    ax[1].step(t_grid_2_plot, v_traj_2, '-', color='C1')

    a_traj = [u_traj_phases[0][0][0]] + [x[0] for x in u_traj_phases[0]]
    ax[2].step(t_grid_phases[0], a_traj, color='C0')

    for i, l in enumerate(['$p$', '$v$', '$a$']):
        ax[i].grid()
        ax[i].set_ylabel(l)

    ax[0].set_xlim([0, T_HORIZON])
    ax[0].legend()


def main_standard_ocp():
    ocp = formulate_double_integrator_ocp()
    ocp.solver_options.tf = T_HORIZON
    ocp.solver_options.nlp_solver_type = 'SQP'

    acados_ocp_solver = AcadosOcpSolver(ocp, verbose=False)

    status = acados_ocp_solver.solve_for_x0(X0)
    acados_ocp_solver.print_statistics()

    x_traj = [acados_ocp_solver.get(i, 'x') for i in range(ocp.dims.N+1)]
    u_traj = [acados_ocp_solver.get(i, 'u') for i in range(ocp.dims.N)]

    t_grid = ocp.solver_options.shooting_nodes

    # plot solution
    nx = len(x_traj[0])
    fig, ax = plt.subplots(nx, 1, sharex=True)

    for i, label in enumerate(["$p$", "$v$"]):
        ax[i].plot(t_grid, np.array(x_traj)[:,i], color='C0')
        ax[i].grid(True)
        ax[i].set_ylabel(label)

    ax[-1].set_xlabel('$t$ [s]')


if __name__ == "__main__":
    main_multiphase_ocp()
    main_standard_ocp()
    plt.show()