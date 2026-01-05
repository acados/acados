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

from acados_template import AcadosModel, AcadosOcp, AcadosOcpSolver, AcadosSim, AcadosSimSolver, mpc_utils
import matplotlib.pyplot as plt

X0 = np.array([2.0, 0.0])
T_HORIZON = 2.0
N_HORIZON = 50
NC_HORIZON = 5

L2_COST_V = 1e-1
L2_COST_P = 1e0
L2_COST_A = 1e-3

T_SIM = 2
TS_SIM = 0.01

def get_double_integrator_model() -> AcadosModel:
    model = AcadosModel()
    model.name = 'double_integrator'
    p = ca.SX.sym('p')
    v = ca.SX.sym('v')
    model.x = ca.vertcat(p, v)
    model.u = ca.SX.sym('u')
    model.f_expl_expr = ca.vertcat(v, model.u)
    return model

def double_integrator_ocp(qp_solver) -> AcadosOcp:
    ocp = AcadosOcp()
    ocp.model = get_double_integrator_model()
    
    ocp.cost.cost_type = 'NONLINEAR_LS'
    ocp.cost.W = np.diag([L2_COST_P, L2_COST_V, L2_COST_A])
    ocp.cost.yref = np.array([0.0, 0.0, 0.0])
    ocp.model.cost_y_expr = ca.vertcat(ocp.model.x, ocp.model.u)

    ocp.constraints.lbu = np.array([-50.])
    ocp.constraints.ubu = np.array([50.])
    ocp.constraints.idxbu = np.array([0])
    ocp.constraints.x0 = X0

    ocp.solver_options.tf = T_HORIZON
    ocp.solver_options.N_horizon = N_HORIZON
    ocp.solver_options.integrator_type = 'ERK'
    ocp.solver_options.nlp_solver_type = 'SQP'
    ocp.solver_options.hessian_approx = 'GAUSS_NEWTON'
    ocp.solver_options.qp_solver = qp_solver

    return ocp

def setup_integrator(dt) -> AcadosSimSolver:
    sim = AcadosSim()
    sim.model = get_double_integrator_model()
    sim.solver_options.T = dt
    return AcadosSimSolver(sim)

def main_ocp_with_ctrl_hor():
    ocp = double_integrator_ocp(qp_solver='FULL_CONDENSING_HPIPM')
    ocp_solver = AcadosOcpSolver(ocp)

    # Initialize Simulation
    sim_solver = setup_integrator(TS_SIM)

    Nsim = int(T_SIM / TS_SIM)
    simT = np.linspace(0, Nsim * TS_SIM, Nsim, endpoint=False)
    simU = np.zeros((Nsim, ocp.dims.nu))
    simX = np.zeros((Nsim+1, ocp.dims.nx))
    simX[0,:] = X0

    time_tot = np.zeros((Nsim))

    # Simulate closed loop (Control Horizon = Prediction Horizon)
    for i in range(Nsim):
        simU[i,:] = ocp_solver.solve_for_x0(simX[i, :])
        time_tot[i] = 1e3 * ocp_solver.get_stats('time_tot')
        simX[i+1, :] = sim_solver.simulate(simX[i, :], simU[i,:])

    # Plot first results
    fig, axs = plt.subplots(4, 1, figsize=(8, 10))
    fig.suptitle('Control Horizon Example (N = ' + str(N_HORIZON) + ')', fontsize=14, weight='bold')
    fig.set_tight_layout(True)

    color_pred_hor = 'tab:red'
    label_pred_hor = 'Nc = ' + str(N_HORIZON)
    axs[0].plot(simT, simX[:-1,0], label=label_pred_hor, color=color_pred_hor)
    axs[1].plot(simT, simX[:-1,1], label=label_pred_hor, color=color_pred_hor)
    axs[2].step(simT, simU, label=label_pred_hor, color=color_pred_hor, where='post')
    axs[3].plot(simT, time_tot, 'o', label=label_pred_hor, color=color_pred_hor, alpha=0.3)

    # Create mocp with control horizon from regular ocp
    ocp = double_integrator_ocp(qp_solver='FULL_CONDENSING_HPIPM')
    mocp = mpc_utils.create_ocp_with_control_horizon(ocp, NC_HORIZON)
    mocp_solver = AcadosOcpSolver(mocp)

    # Simulate closed loop (Control Horizon < Prediction Horizon)
    for i in range(Nsim):
        simU[i,:] = mocp_solver.solve_for_x0(simX[i, :])
        time_tot[i] = 1e3 * mocp_solver.get_stats('time_tot')
        simX[i+1, :] = sim_solver.simulate(simX[i, :], simU[i,:])

        # Plot horizon in first time step
        if i == 0:
            # Get data from mocp
            t = np.linspace(0, T_HORIZON, N_HORIZON+1)
            u = np.empty(N_HORIZON)
            x = np.empty_like(t)
            for j in [k for k in range(N_HORIZON+2) if k != NC_HORIZON-1]:
                if j < NC_HORIZON:
                    x[j] = mocp_solver.get(j, "x")[0]
                    u[j] = mocp_solver.get(j, "u")[0]
                else:
                    x[j-1] = mocp_solver.get(j, "x")[0]
                    if j <= N_HORIZON:
                        u[j-1] = mocp_solver.get(j, "x")[2]

            fig_hor = plt.figure(figsize=(6, 4))
            fig_hor.set_tight_layout(True)

            plt.plot(t, x, label='State', color='tab:blue')
            plt.plot(t, np.zeros_like(t), label='Reference', color='tab:orange', linestyle='--')
            plt.step(t[:-1], u, label='Control', color='tab:red', where='post')

            plt.title('Horizon in first step (t = 0s)')
            plt.xlabel('Time (s)')
            plt.ylabel('Value (m, m/s)')
            plt.xlim([0, T_HORIZON])
            plt.legend()
            plt.grid()

    # Plot final results
    color_ctrl_hor = 'tab:blue'
    label_ctrl_hor = 'Nc = ' + str(NC_HORIZON)
    axs[0].plot(simT, simX[:-1,0], label=label_ctrl_hor, color=color_ctrl_hor)
    axs[1].plot(simT, simX[:-1,1], label=label_ctrl_hor, color=color_ctrl_hor)
    axs[2].step(simT, simU, label=label_ctrl_hor, color=color_ctrl_hor, where='post')
    axs[3].plot(simT, time_tot, 'o', label=label_ctrl_hor, color=color_ctrl_hor, alpha=0.3)

    # Describe axes
    axs[0].set_title('Positions')
    axs[0].set_ylabel('Position (m)')

    axs[1].set_title('Velocities')
    axs[1].set_ylabel('Velocity (m/s)')

    axs[2].set_title('Accelerations')
    axs[2].set_ylabel('Acceleration (m/s²)')

    axs[3].set_title('Execution Times')
    axs[3].set_ylabel('Exec. Time (ms)')

    for ax in axs:
        ax.set_xlabel('Time (s)')
        ax.set_xlim([0, T_SIM])
        ax.grid(True)
        ax.legend()

if __name__ == "__main__":
    main_ocp_with_ctrl_hor()
    plt.show()
