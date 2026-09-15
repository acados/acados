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

from acados_template import AcadosModel, AcadosOcp, AcadosOcpSolver, casadi_length, latexify_plot, plot_trajectories
import matplotlib.pyplot as plt

latexify_plot()

X0 = np.array([2.0, 0.0])
PENALTY_X = 1e0
T_HORIZON = 1.0
N_HORIZON = 20

L2_COST_V = 1e-1
L2_COST_P = 1e0
L2_COST_A = 1e-3

V_MAX = 1.5
U_MAX = 3.0

def get_double_integrator_model() -> AcadosModel:
    model = AcadosModel()
    model.name = 'double_integrator'
    p = ca.SX.sym('p')
    v = ca.SX.sym('v')
    model.x = ca.vertcat(p, v)
    model.u = ca.SX.sym('u')
    model.f_expl_expr = ca.vertcat(v, model.u)
    return model


def formulate_double_integrator_ocp() -> AcadosOcp:
    ocp = AcadosOcp()
    nx = 2

    ocp.model = get_double_integrator_model()

    ocp.cost.cost_type = 'NONLINEAR_LS'
    ocp.cost.W = np.diag([L2_COST_P, L2_COST_V, L2_COST_A])
    ocp.cost.yref = np.array([0.0, 0.0, 0.0])

    ocp.model.cost_y_expr = ca.vertcat(ocp.model.x, ocp.model.u)

    ocp.constraints.lbu = np.array([-U_MAX])
    ocp.constraints.ubu = np.array([U_MAX])
    ocp.constraints.idxbu = np.array([0])

    ocp.constraints.idxbx = np.array([1])
    ocp.constraints.lbx = np.array([-V_MAX])
    ocp.constraints.ubx = np.array([V_MAX])

    # slack u at initial stage
    ocp.constraints.idxs_rev_0 = np.array([0] + nx * [-1])
    ocp.cost.Zl_0 = 2 * np.array([1.0])
    ocp.cost.zl_0 = 0 * np.array([1.0])
    ocp.cost.Zu_0 = 2 * np.array([1.0])
    ocp.cost.zu_0 = 0 * np.array([1.0])

    # slack u and x at intermediate stages
    ocp.constraints.idxs_rev = np.array([0, 1])
    ocp.cost.Zl = 2 * np.array(2*[1.0])
    ocp.cost.zl = 0 * np.array(2*[1.0])
    ocp.cost.Zu = 2 * np.array(2*[1.0])
    ocp.cost.zu = 0 * np.array(2*[1.0])

    ocp.constraints.x0 = X0
    return ocp


def main(qp_solver = 'PARTIAL_CONDENSING_HPIPM'):
    solutions = []
    variants = []
    nx, nu = 2, 1

    ocp = formulate_double_integrator_ocp()
    ocp.name = 'slack_formulations'
    ocp.solver_options.N_horizon = N_HORIZON
    ocp.solver_options.tf = T_HORIZON
    ocp.solver_options.nlp_solver_max_iter = 1
    ocp.solver_options.tol = 1e-10
    ocp.solver_options.eval_residual_at_max_iter = True
    ocp.solver_options.qp_solver = qp_solver

    ocp_solver = AcadosOcpSolver(ocp)

    variants = ['xu individual slack', 'all hard', 'xu joint slack', 'u hard', 'x hard']
    for variant in variants:
        idx_unused_slack = []
        if variant == 'xu individual slack':
            ocp_solver.constraints_set(0, 'idxs_rev', ocp.constraints.idxs_rev_0)
            for stage in range(1, N_HORIZON):
                ocp_solver.constraints_set(stage, 'idxs_rev', ocp.constraints.idxs_rev)

        elif variant == 'all hard':
            ocp_solver.constraints_set(0, 'idxs_rev', np.array((nx+nu) * [-1]))
            for stage in range(1, N_HORIZON):
                ocp_solver.constraints_set(stage, 'idxs_rev', np.array(2 * [-1]))
            idx_unused_slack = [0, 1]

        elif variant == 'xu joint slack':
            ocp_solver.constraints_set(0, 'idxs_rev', ocp.constraints.idxs_rev_0)
            for stage in range(1, N_HORIZON):
                ocp_solver.constraints_set(stage, 'idxs_rev', np.array(2 * [0]))
            idx_unused_slack = [1]

        elif variant == 'u hard':
            ocp_solver.constraints_set(0, 'idxs_rev', np.array((nx+nu) * [-1]))
            for stage in range(1, N_HORIZON):
                ocp_solver.constraints_set(stage, 'idxs_rev', np.array([-1, 0]))
            idx_unused_slack = [1]

        elif variant == 'x hard':
            ocp_solver.constraints_set(0, 'idxs_rev', ocp.constraints.idxs_rev_0)
            for stage in range(1, N_HORIZON):
                ocp_solver.constraints_set(stage, 'idxs_rev', np.array([0, -1]))
            idx_unused_slack = [1]

        # solve and store results
        print(f"\nSolving with variant: {variant}\n")

        # ocp_solver.options_set('print_level', 4)
        # ocp_solver.reset()
        ocp_solver.solve()
        ocp_solver.print_statistics()

        sol = ocp_solver.get_iterate()
        solutions.append(sol)
        # print(f"{sol.su=}")
        # print(f"{sol.sl=}")
        # print(f"{sol.x=}")
        # print(f"{sol.lam=}")

        ## sanity checks
        assert ocp_solver.status == 0, f"Got status {ocp_solver.status}, should be 0."

        # u bound
        u_traj = np.array(sol.u)
        if variant in ['all hard', 'u hard']:
            np.testing.assert_array_less(u_traj, U_MAX)
            np.testing.assert_array_less(-U_MAX, u_traj)
        else:
            assert not (np.all(u_traj<U_MAX) or np.all(-U_MAX<u_traj)), "some u bound softening should be exploited, otherwise test problem is bad"

        # x bound
        v_traj = np.array([sol.x[i][1] for i in range(1, N_HORIZON)])
        if variant in ['all hard', 'x hard']:
            np.testing.assert_array_less(v_traj, V_MAX)
            np.testing.assert_array_less(-V_MAX, v_traj)
        else:
            if not (np.all(v_traj<V_MAX) or np.all(-V_MAX<v_traj)):
                raise ValueError(f"some v bound softening should be exploited, otherwise test problem is bad. Got {v_traj=}")

        # check stuff at intermediate nodes
        nb = ocp.dims.nbx + ocp.dims.nbu
        ns = ocp.dims.ns
        for stage in range(1, N_HORIZON):
            lam_l = sol.lam[stage][:nb]
            lam_u = sol.lam[stage][nb:2*nb]
            lam_sl = sol.lam[stage][2*nb:2*nb+ns]
            lam_su = sol.lam[stage][2*nb+ns:2*nb+2*ns]

            # check unused slacks are 0
            for i in idx_unused_slack:
                assert np.abs(sol.sl[stage][i]) < 1e-8, f"unused sl[{stage}][{i}] = {sol.sl[stage][i]} != 0 for variant {variant}"
                assert np.abs(sol.su[stage][i]) < 1e-8, f"unused su[{stage}][{i}] = {sol.su[stage][i]} != 0 for variant {variant}"

                assert np.abs(lam_sl[i]) < 1e-8, f"unused slack bound multiplier lam_sl[{stage}][{i}] = {lam_sl[i]} != 0 for variant {variant}"
                assert np.abs(lam_su[i]) < 1e-8, f"unused slack bound multiplier lam_su[{stage}][{i}] = {lam_su[i]} != 0 for variant {variant}"

            if variant == 'xu joint slack':
                # and stage in [1]:
                u_val = sol.u[stage][0]
                v_val = sol.x[stage][1]

                max_lower_violation = max(-U_MAX - u_val, -V_MAX - v_val, 0.0)
                max_upper_violation = max(u_val - U_MAX, v_val - V_MAX, 0.0)
                sl_difference = sol.sl[stage][0] - max_lower_violation
                su_difference = sol.su[stage][0] - max_upper_violation

                # print(f"\nChecking stage {stage}, x: {sol.x[stage]}, u: {sol.u[stage]},\n sl: {sol.sl[stage]}, {max_lower_violation=}, {sl_difference=},\n su: {sol.su[stage]}, {max_upper_violation=}, {su_difference=}")

                # NOTE: sometimes slack take larger values than constraint violation
                # But all conditions are satisfied.
                tol_s = 5e-5
                np.testing.assert_allclose(sol.sl[stage][0], max_lower_violation, atol=tol_s, rtol=tol_s, err_msg=f"sl should match max lower violation at stage {stage}")
                np.testing.assert_allclose(sol.su[stage][0], max_upper_violation, atol=tol_s, rtol=tol_s, err_msg=f"su should match max upper violation at stage {stage}\n")

                sl_stat = ocp.solver_options.cost_scaling[stage] * (ocp.cost.zl[0] + ocp.cost.Zl[0] * sol.sl[stage][0]) - sum(lam_l) - lam_sl[0]
                su_stat = ocp.solver_options.cost_scaling[stage] * (ocp.cost.zu[0] + ocp.cost.Zu[0] * sol.su[stage][0]) - sum(lam_u) - lam_su[0]

                if np.abs(sl_stat) > ocp.solver_options.tol or np.abs(su_stat) > ocp.solver_options.tol:
                    raise ValueError(f"Stationarity wrt slack variables not satisfied got: {sl_stat=}, {su_stat=}")

    # compare and eval
    nvariants = len(variants)

    plot_trajectories(
        x_traj_list=[np.array(s.x) for s in solutions],
        u_traj_list=[np.array(s.u) for s in solutions],
        time_traj_list=nvariants*[np.linspace(0, ocp.solver_options.tf, ocp.solver_options.N_horizon+1)],
        time_label=ocp.model.t_label,
        labels_list=variants,
        x_labels=ocp.model.x_labels,
        u_labels=ocp.model.u_labels,
        idxbu=ocp.constraints.idxbu,
        lbu=ocp.constraints.lbu,
        ubu=ocp.constraints.ubu,
        linestyle_list=3*['-', '--', ':', '-.'],
        X_ref=None,
        U_ref=None,
        x_min=None,
        x_max=None,
    )

if __name__ == "__main__":
    main('FULL_CONDENSING_HPIPM')
    plt.show()
    main('PARTIAL_CONDENSING_HPIPM')
