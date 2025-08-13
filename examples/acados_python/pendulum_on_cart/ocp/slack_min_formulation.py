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

import sys
sys.path.insert(0, '../common')

from acados_template import AcadosOcp, AcadosOcpSolver, plot_trajectories, AcadosModel, ACADOS_INFTY
from pendulum_model import export_pendulum_ode_model
import numpy as np
import casadi as ca

def main(formulation='s_slack', plot_traj=True):
    # create ocp object to formulate the OCP
    ocp = AcadosOcp()

    # set model
    model: AcadosModel = export_pendulum_ode_model()
    ocp.model = model

    Tf = 1.0
    N = 10

    # set prediction horizon
    ocp.solver_options.N_horizon = N
    ocp.solver_options.tf = Tf

    # cost matrices
    Q_mat = 2*np.diag([1e3, 1e3, 1e-2, 1e-2])
    R_mat = 2*np.diag([1e-2])

    # path cost
    ocp.cost.cost_type = 'EXTERNAL'
    W_mat = ca.diagcat(Q_mat, R_mat).full()
    ocp.model.cost_expr_ext_cost = .5*ca.vertcat(model.x, model.u).T @ W_mat @ ca.vertcat(model.x, model.u)
    # terminal cost
    ocp.cost.cost_type_e = 'EXTERNAL'
    ocp.model.cost_expr_ext_cost_e = .5 * model.x.T @ Q_mat @ model.x

    # set constraints
    Fmax = 40
    ocp.constraints.lbu = np.array([-Fmax])
    ocp.constraints.ubu = np.array([+Fmax])
    ocp.constraints.idxbu = np.array([0])

    # initial condition
    ocp.constraints.x0 = np.array([0.0, 0.2 * np.pi, 0.0, 0.0])

    print(f"using formulation {formulation}")
    # add cost term min(x[0], x[2]) to cost
    # via slack: s <= x[0], s <= x[2]
    if formulation == 'u_slack':
        # add u
        new_u = ca.SX.sym('new_u', 1, 1)
        ocp.model.u = ca.vertcat(model.u, new_u)
        ocp.model.u_labels.append('new_u')
        # add constraints u <= x_i
        #  <=> -inf <= u - x_i <= 0
        ocp.model.con_h_expr = ca.vertcat(new_u - model.x[0], new_u - model.x[3])
        ocp.constraints.uh = np.zeros((2, 1))
        ocp.constraints.lh = - ACADOS_INFTY * np.ones((2, 1))

        ocp.model.con_h_expr_0 = ocp.model.con_h_expr
        ocp.constraints.uh_0 = ocp.constraints.uh
        ocp.constraints.lh_0 = ocp.constraints.lh

        # add cost -u
        ocp.model.cost_expr_ext_cost -= new_u
    elif formulation == 'u_slack2':
        # add u
        new_u = ca.SX.sym('new_u', 1, 1)
        ocp.model.u = ca.vertcat(model.u, new_u)
        ocp.model.u_labels.append('new_u')
        # add constraints u <= x_i
        #  <=> -inf <= u - x_i <= 0
        ocp.model.con_h_expr = ca.vertcat(new_u + model.x[0], new_u + model.x[3])
        ocp.constraints.uh = ACADOS_INFTY * np.ones((2, 1))
        ocp.constraints.lh = 0 * np.ones((2, 1))

        ocp.model.con_h_expr_0 = ocp.model.con_h_expr
        ocp.constraints.uh_0 = ocp.constraints.uh
        ocp.constraints.lh_0 = ocp.constraints.lh

        # add cost u
        ocp.model.cost_expr_ext_cost += new_u
    elif formulation == "s_slack":
        # add s
        ns = 1
        # add constraints: s <= x_i
        ocp.model.con_h_expr = ca.vertcat(model.x[0], model.x[3])
        ocp.constraints.uh = ACADOS_INFTY * np.ones((2, 1))
        ocp.constraints.lh = np.zeros((2, 1))
        ocp.constraints.idxs_rev = np.array([-1, 0, 0])
        ocp.constraints.ls = -ACADOS_INFTY * np.ones((ns, ))
        ocp.constraints.us = 0 * np.ones((ns, ))
        ocp.cost.zl = np.array([1.0])
        ocp.cost.Zl = np.array([-0.0])
        ocp.cost.zu = np.array([1.0])
        ocp.cost.Zu = np.array([0.0])

        ocp.model.con_h_expr_0 = ocp.model.con_h_expr
        ocp.constraints.uh_0 = ocp.constraints.uh
        ocp.constraints.lh_0 = ocp.constraints.lh
        ocp.cost.zl_0 = ocp.cost.zl
        ocp.cost.Zl_0 = ocp.cost.Zl
        ocp.cost.zu_0 = ocp.cost.zu
        ocp.cost.Zu_0 = ocp.cost.Zu
        nbx_0 = ocp.constraints.lbx_0.shape[0]
        nbu = ocp.constraints.lbu.shape[0]
        ocp.constraints.idxs_rev_0 = np.array((nbx_0+nbu) * [-1] + [0, 0])
        ocp.constraints.ls_0 = ocp.constraints.ls
        ocp.constraints.us_0 = ocp.constraints.us
    ocp.solver_options.qp_solver_t0_init = 0

    # set options
    ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM'
    ocp.solver_options.hessian_approx = 'GAUSS_NEWTON'
    ocp.solver_options.integrator_type = 'IRK'
    ocp.solver_options.nlp_solver_type = 'SQP'
    # ocp.solver_options.print_level = 5
    # ocp.solver_options.nlp_solver_max_iter = 2


    nx = model.x.rows()
    nu = model.u.rows()

    ocp_solver = AcadosOcpSolver(ocp, verbose=False)

    xtraj = np.zeros((N+1, nx))
    utraj = np.zeros((N, nu))

    status = ocp_solver.solve()
    ocp_solver.print_statistics()

    if status != 0:
        raise Exception(f'acados returned status {status}.')

    # get solution
    for i in range(N):
        xtraj[i,:] = ocp_solver.get(i, "x")
        utraj[i,:] = ocp_solver.get(i, "u")
    xtraj[N,:] = ocp_solver.get(N, "x")

    min_x_vals = np.minimum(xtraj[:, 0], xtraj[:, 3])
    if formulation == 'u_slack':
        slack_vals = utraj[:, 1]
        assert np.allclose(min_x_vals[:-1], slack_vals, atol=1e-6)
    elif formulation == 'u_slack2':
        slack_vals = utraj[:, 1]
        assert np.allclose(min_x_vals[:-1], -slack_vals, atol=1e-6)
    elif formulation == 's_slack':
        slack_vals = np.zeros((N, ))
        unused_slack_vals = np.zeros((N, ))
        for i in range(N):
            slack_vals[i] = ocp_solver.get(i, "sl")[0]
            unused_slack_vals[i] = ocp_solver.get(i, "su")[0]
        assert np.allclose(min_x_vals[:-1], -slack_vals, atol=1e-6)
        print(f"{unused_slack_vals=}")
        # plot slacks
        utraj = np.append(utraj, np.atleast_2d(slack_vals).transpose(), axis=1)
        model.u_labels.append('slack')

    if plot_traj:
        plot_trajectories(
            x_traj_list=[xtraj],
            u_traj_list=[utraj],
            time_traj_list=[np.linspace(0, Tf, N+1)],
            time_label=model.t_label,
            labels_list=['OCP result'],
            x_labels=model.x_labels,
            u_labels=model.u_labels,
            idxbu=ocp.constraints.idxbu,
            lbu=ocp.constraints.lbu,
            ubu=ocp.constraints.ubu,
            X_ref=None,
            U_ref=None,
            # fig_filename='pendulum_ocp.png',
            x_min=None,
            x_max=None,
        )

    return xtraj


if __name__ == '__main__':
    formulations = ['u_slack', 'u_slack2', 's_slack']
    # formulations = ['s_slack']
    xtraj_list = []
    for i, formulation in enumerate(formulations):
        xtraj = main(formulation, plot_traj=False)
        if i == 0:
            xtraj_ref = xtraj
        else:
            diff_x = np.max(np.abs(xtraj - xtraj_ref))
            print(f"diff xtraj {formulation} vs ref: {diff_x}")
            if diff_x > 1e-6:
                raise Exception(f"xtraj {formulation} differs from reference by {diff_x}, expected to be close to zero.")
