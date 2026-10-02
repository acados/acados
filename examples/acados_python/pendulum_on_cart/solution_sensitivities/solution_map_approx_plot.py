#
# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# Licensed under the 2-Clause BSD License.

import sys
import os
sys.path.insert(0, '../common')

from acados_template import AcadosOcp, AcadosOcpSolver, AcadosOcpQp, AcadosOcpQpOptions, plot_trajectories, AcadosOcpQpSolver, latexify_plot
from pendulum_model import export_pendulum_ode_model
import numpy as np
import scipy.linalg
import matplotlib.pyplot as plt

FMAX = 80
T_HORIZON = 2.0
N = 20

def create_ocp() -> AcadosOcp:
    # create ocp object to formulate the OCP
    ocp = AcadosOcp()

    # set model
    model = export_pendulum_ode_model()
    ocp.model = model

    nx = model.x.rows()
    nu = model.u.rows()
    ny = nx + nu
    ny_e = nx

    # set dimensions
    ocp.solver_options.N_horizon = N

    # set cost
    Q = 2*np.diag([1e1, 1e2, 1e-2, 5e-3])
    R = 2*np.diag([1e-2])

    ocp.cost.W_e = Q
    ocp.cost.W = scipy.linalg.block_diag(Q, R)

    ocp.cost.cost_type = 'LINEAR_LS'
    ocp.cost.cost_type_e = 'LINEAR_LS'

    ocp.cost.Vx = np.zeros((ny, nx))
    ocp.cost.Vx[:nx,:nx] = np.eye(nx)

    Vu = np.zeros((ny, nu))
    Vu[4,0] = 1.0
    ocp.cost.Vu = Vu

    ocp.cost.Vx_e = np.eye(nx)

    ocp.cost.yref = np.zeros((ny, ))
    ocp.cost.yref_e = np.zeros((ny_e, ))

    # set constraints
    ocp.constraints.lbu = np.array([-FMAX])
    ocp.constraints.ubu = np.array([+FMAX])
    ocp.constraints.idxbu = np.array([0])

    VMAX = 1.
    ocp.constraints.lbx = np.array([-VMAX])
    ocp.constraints.ubx = np.array([+VMAX])
    ocp.constraints.idxbx = np.array([2])

    ocp.constraints.idxsbx = np.array([0])
    ns = 1
    ocp.cost.zl = 5e1 * np.ones((ns,))
    ocp.cost.zu = 5e1 * np.ones((ns,))
    ocp.cost.Zl = 0 * np.ones((ns,))
    ocp.cost.Zu = 0 * np.ones((ns,))


    ocp.constraints.lbx_e = np.array([-VMAX])
    ocp.constraints.ubx_e = np.array([+VMAX])
    ocp.constraints.idxbx_e = np.array([2])

    ocp.constraints.idxsbx_e = np.array([0])
    ocp.cost.zl_e = 5e1 * np.ones((ns,))
    ocp.cost.zu_e = 5e1 * np.ones((ns,))
    ocp.cost.Zl_e = 0 * np.ones((ns,))
    ocp.cost.Zu_e = 0 * np.ones((ns,))

    ocp.constraints.x0 = np.array([0.0, np.pi, 0.0, 0.0])

    # set options
    ocp.solver_options.qp_solver = 'FULL_CONDENSING_HPIPM'
    ocp.solver_options.hessian_approx = 'GAUSS_NEWTON'
    ocp.solver_options.integrator_type = 'ERK'
    ocp.solver_options.nlp_solver_type = 'SQP'
    ocp.solver_options.qp_solver_iter_max = 1000
    ocp.solver_options.nlp_solver_max_iter = 1000
    ocp.solver_options.nlp_qp_tol_reduction_factor = 1e-3
    ocp.solver_options.qp_solver_mu0 = 1e2
    ocp.solver_options.tf = T_HORIZON

    ocp.code_gen_options.code_export_directory = 'c_generated_code'

    return ocp

def plot_solution_maps(p_list, sol_list, labels, special_points=None,
                      special_point_labels=None, fig_filename=None):
    latexify_plot()

    plt.figure()
    linestyles = 2*['-', '--', '--', ':', ':']
    for p, sol, label, linestyle in zip(p_list, sol_list, labels, linestyles):
        plt.plot(p, sol, label=label, linestyle=linestyle)
    if special_points is not None:
        for i, (p, u) in enumerate(special_points):
            label = special_point_labels[i] if special_point_labels else None
            plt.scatter(p, np.asarray(u).squeeze(), label=label, zorder=3, marker='X', color='red')
    plt.xlabel('$x_0$ projected')
    plt.ylabel('$u_0$')
    plt.xlim(min(np.min(p) for p in p_list), max(np.max(p) for p in p_list))
    plt.grid(True)
    plt.legend()

    if fig_filename is not None:
        plt.savefig(fig_filename)


def solve_for_x0_vals(solver, x0_vals, label):
    n_vals = x0_vals.shape[0]
    nu = 1
    u0_vals = np.zeros((n_vals, nu))
    for i in range(n_vals):
        x0 = x0_vals[i, :]
        u0 = solver.solve_for_x0(x0, fail_on_nonzero_status=False)
        if solver.status != 0:
            raise Exception(f"solver {label} failed for instance {i}")
        u0_vals[i, :] = u0
    return u0_vals


def solve_qp_for_delta_x0_vals(solver, x0_vals, label):
    n_vals = x0_vals.shape[0]
    nu = 1
    u0_vals = np.zeros((n_vals, nu))
    for i in range(n_vals):
        x0 = x0_vals[i, :]
        solver.set(0, 'lbx', x0)
        solver.set(0, 'ubx', x0)
        solver.solve()
        u0 = solver.get(0, 'u')
        if solver.status != 0:
            raise Exception(f"solver {label} failed for instance {i}")
        u0_vals[i, :] = u0
    return u0_vals

def main(with_sbqp_approx=False,
        with_qp_approx=False,
        with_bqp_approx=False,
        with_linearization_error=False,
        plot_traj=False,
        fig_filename=None):
    ocp = create_ocp()
    ocp_solver = AcadosOcpSolver(ocp)

    nx, nu = 4, 1
    n_vals = 1000
    x0_vals = np.zeros((n_vals, nx))
    x0_vals[:, 0] = np.linspace(-3, 8, n_vals)
    x0_vals[:, 2] = -0.6
    x0_vals[:, 1] = .2
    p_vals = x0_vals[:, 0]

    u0_list = []
    labels = []

    # exact
    label = 'exact'
    u0_vals = solve_for_x0_vals(ocp_solver, x0_vals, label)
    u0_list.append(u0_vals)
    labels.append(label)

    # SBQP
    if with_sbqp_approx:
        for tau_min in [0.05, 0.1]:
            ocp_solver.options_set('tau_min', tau_min)
            label = r'SBQP $\tau_{\mathrm{min}} = ' + f"{tau_min}" + '$'
            u0_vals = solve_for_x0_vals(ocp_solver, x0_vals, label)
            u0_list.append(u0_vals)
            labels.append(label)

    # setup QP approx
    xlin_idx = np.argmin(np.abs(x0_vals[:, 0] - 1.2))
    x_lin = x0_vals[xlin_idx, :]
    # ocp_solver.options_set("levenberg_marquardt", 10.)
    if with_linearization_error:
        ocp_solver.options_set("nlp_solver_max_iter", 3)
    u0_offset = ocp_solver.solve_for_x0(x_lin, fail_on_nonzero_status=False)

    ocp_solver.print_statistics()
    qp_dict = ocp_solver.get_last_qp()
    qp = AcadosOcpQp.from_dict(qp_dict)
    qp_opts = AcadosOcpQpOptions()
    delta_x0_vals = x0_vals - x_lin

    # QP approx
    if with_qp_approx:
        label = 'QP approx.'
        qp_solver = AcadosOcpQpSolver(qp, qp_opts)

        u0_vals = np.zeros((n_vals, nu))

        u0_vals = solve_qp_for_delta_x0_vals(qp_solver, delta_x0_vals, label)
        u0_vals += u0_offset
        u0_list.append(u0_vals)
        labels.append(label)
        u_lin = u0_vals[xlin_idx]

    # BQP approx
    if with_bqp_approx:
        for tau_min in [0.1]:
            label = r'BQP approx. $\tau_{\mathrm{min}} = ' + f"{tau_min}" + '$'
            qp_opts = AcadosOcpQpOptions()
            # qp_opts.print_level = 1
            qp_solver = AcadosOcpQpSolver(qp, qp_opts)
            qp_solver.opts_set('tau_min', tau_min)

            u0_vals = solve_qp_for_delta_x0_vals(qp_solver, delta_x0_vals, label)
            u0_vals += u0_offset
            u0_list.append(u0_vals)
            labels.append(label)

    p_vals_list = [p_vals] * len(u0_list)
    if with_bqp_approx or with_qp_approx:
        special_points = [(x_lin[0], u_lin)]
        special_point_labels = ['linearization']
    else:
        special_points, special_point_labels = None, None
    plot_solution_maps(p_vals_list, u0_list, labels,
                       special_points=special_points,
                       special_point_labels=special_point_labels,
                       fig_filename=fig_filename)
    sol = ocp_solver.get_iterate()

    if plot_traj:
        plot_trajectories(
            x_traj_list=[np.array(sol.x)],
            u_traj_list=[np.array(sol.u)],
            time_traj_list=[np.linspace(0, T_HORIZON, N+1)],
            time_label=ocp.model.t_label,
            labels_list=['OCP result'],
            x_labels=ocp.model.x_labels,
            u_labels=ocp.model.u_labels,
            idxbu=ocp.constraints.idxbu,
            lbu=ocp.constraints.lbu,
            ubu=ocp.constraints.ubu,
            X_ref=None,
            U_ref=None,
            x_min=None,
            x_max=None,
        )


if __name__ == "__main__":
    main(with_sbqp_approx=True, with_linearization_error=True, with_bqp_approx=True, with_qp_approx=True)
    # main(with_sbqp_approx=True, fig_filename='solution_map_sbqp.pdf')
    # main(with_linearization_error=False, with_qp_approx=True, with_bqp_approx=True,
    #     fig_filename='solution_map_qp.pdf')
    # main(with_linearization_error=True, with_qp_approx=True, with_bqp_approx=True,
    #     fig_filename='solution_map_qp_error.pdf')

    plt.show()
