#
# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# Licensed under the 2-Clause BSD License.

import numpy as np
import casadi as ca
from acados_template import AcadosOcp, AcadosModel, AcadosOcpSolver, latexify_plot, ACADOS_INFTY
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
latexify_plot()

PROBLEM_NAME = "non_ocp_slacks"

def create_parametric_nlp(mode='forw', one_sided=False) -> AcadosOcp:

    model = AcadosModel()
    model.x = ca.SX.sym("x", 1)
    model.p_global = ca.SX.sym("p_global", 1)
    model.cost_expr_ext_cost_e = (model.x - model.p_global)**2
    model.name = "non_ocp"
    ocp = AcadosOcp()
    ocp.model = model
    ocp.name = f"{PROBLEM_NAME}_{one_sided}_{mode}"

    ocp.constraints.lbx_e = np.array([-1.0])
    ocp.constraints.ubx_e = np.array([1.0])
    if one_sided:
        ocp.constraints.ubx_e = np.array([ACADOS_INFTY])
    ocp.constraints.idxbx_e = np.array([0])

    ocp.constraints.idxs_rev_e = np.array([0])
    ocp.model.cost_zl_expr_e = 1 + 0.2 * ca.fabs(model.p_global)**2
    ocp.model.cost_zu_expr_e = 1 + 0.2 * ca.fabs(model.p_global)**2
    ocp.model.cost_Zl_expr_e = 1 + 0.2 * ca.fabs(model.p_global)**2
    ocp.model.cost_Zu_expr_e = 1 + 0.2 * ca.fabs(model.p_global)**2

    # ocp.cost.zl_e = np.array([1.0])
    # ocp.cost.zu_e = np.array([1.0])
    # ocp.cost.Zl_e = np.array([0.0])
    # ocp.cost.Zu_e = np.array([0.0])

    ocp.cost.cost_type_e = "EXTERNAL"
    ocp.solver_options.qp_solver = "FULL_CONDENSING_HPIPM"
    ocp.solver_options.qp_solver_ric_alg = 0
    ocp.solver_options.hessian_approx = "EXACT"
    ocp.solver_options.N_horizon = 0

    ocp.p_global_values = np.zeros((1,))
    if mode == 'forw':
        ocp.code_gen_options.with_solution_sens_wrt_params_forw = True
    else:
        ocp.code_gen_options.with_solution_sens_wrt_params_adj = True
    ocp.code_gen_options.with_value_sens_wrt_params = True
    ocp.solver_options.nlp_solver_ext_qp_res = 1

    return ocp

def solve_and_compute_sens(p_test, tau, mode='adj', one_sided=False):
    np_test = p_test.shape[0]

    ocp = create_parametric_nlp(mode, one_sided=one_sided)
    ocp.solver_options.tau_min = tau
    ocp.solver_options.qp_solver_t0_init = 0
    ocp.solver_options.nlp_solver_ext_qp_res = 1
    ocp.solver_options.nlp_solver_max_iter = 2 # QP should converge in one iteration

    ocp_solver = AcadosOcpSolver(ocp, verbose=False)

    sens_x = np.zeros(np_test)
    solution = np.zeros(np_test)

    seed = [(0, np.array([[1.0]]))]

    for i, p in enumerate(p_test):
        p_val = np.array([p])

        ocp_solver.set_p_global_and_precompute_dependencies(p_val)
        status = ocp_solver.solve()
        solution[i] = ocp_solver.get(0, "x")[0]

        if status != 0:
            ocp_solver.print_statistics()
            raise Exception(f"OCP solver returned status {status} at {i}th p value {p}, {tau=}.")
            # print(f"OCP solver returned status {status} at {i}th p value {p}, {tau=}.")
            # breakpoint()

        status = ocp_solver.setup_qp_matrices_and_factorize()
        if status != 0:
            ocp_solver.print_statistics()
            raise Exception(f"OCP solver returned status {status} in setup_qp_matrices_and_factorize at {i}th p value {p}, {tau=}.")

        # Calculate the policy gradient
        if mode == 'forw':
            out_dict = ocp_solver.eval_solution_sensitivity(0, "p_global", return_sens_x=True, return_sens_u=False)
            sens_x[i] = out_dict['sens_x'].item()
        elif mode == 'adj':
            sens_adj = ocp_solver.eval_adjoint_solution_sensitivity(seed_x=seed, seed_u=None)
            sens_x[i] = sens_adj.item()

    return solution, sens_x

P_MIN = -4.0
P_MAX = 4.0
def main(mode='adj', one_sided=False):
    if mode not in ['forw', 'adj']:
        raise ValueError(f"mode must be forw, or adj, got {mode}.")
    print(f"computing sensitivities via {mode} mode.")
    delta_p = 0.002
    p_test = np.arange(P_MIN, P_MAX, delta_p)
    sens_list = []
    labels_list = []
    sol_list = []
    test_tol = 1e-3
    tau_vals = [1e-6, 1e-4, 1e-2]
    for tau in tau_vals:
        sol_tau, sens_x_tau = solve_and_compute_sens(p_test, tau, mode=mode, one_sided=one_sided)

        # Compare to numerical gradients; derivatives may differ at active-set changes.
        sens_x_fd = np.gradient(sol_tau, delta_p)
        q95_diff = np.quantile(np.abs(sens_x_tau - sens_x_fd), 0.95)
        print(f"95th-percentile difference for tau={tau} between policy gradient obtained by acados and via FD is {q95_diff:.2e} should be < {test_tol:.2e}.")
        assert q95_diff <= test_tol

        sens_list.append(sens_x_tau)
        labels_list.append(r"$\tau = 10^{" + f"{int(np.log10(tau))}" + r"}$")
        # labels_list.append(r"$\tau =" + f"{tau}" + r"$")
        sol_list.append(sol_tau)

    plot_solution_sensitivities_results(p_test, sol_list, sens_list, labels_list,
                 title=None, parameter_name=r"$\theta$", fig_filename=f"solution_sens_{PROBLEM_NAME}.pdf")


def main_value_gradient(one_sided=False):
    p_test = np.linspace(P_MIN, P_MAX, 401)
    ocp = create_parametric_nlp(mode='forw', one_sided=one_sided)
    ocp.solver_options.tau_min = 0.0
    ocp.solver_options.qp_solver_t0_init = 0
    ocp.solver_options.nlp_solver_max_iter = 2
    ocp_solver = AcadosOcpSolver(ocp, verbose=False)

    optimal_value = np.zeros(p_test.size)
    optimal_value_grad = np.zeros(p_test.size)
    for i, p in enumerate(p_test):
        ocp_solver.set_p_global_and_precompute_dependencies(np.array([p]))
        status = ocp_solver.solve()
        if status != 0:
            ocp_solver.print_statistics()
            raise RuntimeError(f"OCP solver returned status {status} at p={p}, tau=0.")

        optimal_value[i] = ocp_solver.get_cost()
        optimal_value_grad[i] = ocp_solver.eval_and_get_optimal_value_gradient("p_global").item()

    optimal_value_grad_fd = np.gradient(optimal_value, p_test)
    mean_diff = np.mean(np.abs(optimal_value_grad - optimal_value_grad_fd))
    test_tol = 1e-4
    print(f"mean difference between value sensitivities at tau=0 and finite differences is "
          f"{mean_diff:.2e} should be < {test_tol:.2e}.")

    _, axes = plt.subplots(nrows=2, ncols=1, sharex=True, figsize=(6.5, 5))
    axes[0].plot(p_test, optimal_value)
    axes[0].set_ylabel(r"optimal value $V^*(p)$")
    axes[1].plot(p_test, optimal_value_grad, label="acados")
    axes[1].plot(p_test, optimal_value_grad_fd, linestyle=":", label="finite differences")
    axes[1].set_ylabel(r"value gradient $\partial_p V^*(p)$")
    axes[1].set_xlabel(r"parameter $p$")
    for axis in axes:
        axis.grid(True)
    axes[1].legend()
    axes[1].set_xlim(P_MIN, P_MAX)
    plt.tight_layout()
    plt.show()

    assert mean_diff <= test_tol


def plot_solution_sensitivities_results(p_test, sol_list, sens_list, labels_list, title=None, parameter_name="", fig_filename=None, horizontal_plot=False):
    linestyles = ["--", "-.", "--", ":", "-.", ":"]

    nsub = 2
    if horizontal_plot:
        _, ax = plt.subplots(nrows=1, ncols=nsub, sharex=False, figsize=(12, 3.0))
    else:
        _, ax = plt.subplots(nrows=nsub, ncols=1, sharex=True, figsize=(6.5,5))

    # plot solutions
    isub = 0
    for i, sol in enumerate(sol_list):
        ax[isub].plot(p_test, sol, label=labels_list[i], linestyle=linestyles[i])
    ax[isub].set_xlim([p_test[0], p_test[-1]])
    ax[isub].set_ylabel(r"solution $x^{\star}$")
    if title is not None:
        ax[isub].set_title(title)
    ax[isub].legend()

    # plot sensitivities
    isub += 1
    for i, sens_x_tau in enumerate(sens_list):
        ax[isub].plot(p_test, sens_x_tau, color=f"C{i}", linestyle="--")
        sens_fd = np.gradient(sol_list[i], p_test)
        ax[isub].plot(p_test, sens_fd, color=f"C{i}", linestyle=":")
    ax[isub].set_xlim([P_MIN, P_MAX])
    ax[isub].set_ylabel(r"derivative $\partial_\theta x^{\star}$")
    legend_handles = [
        Line2D([], [], color="black", linestyle="--", label="acados"),
        Line2D([], [], color="black", linestyle=":", label="finite differences"),
    ]
    ax[isub].legend(handles=legend_handles)

    for i in range(nsub):
        ax[i].grid(True)
        if horizontal_plot:
            ax[i].set_xlabel(f"{parameter_name}")
    ax[-1].set_xlabel(f"{parameter_name}")

    plt.tight_layout()

    if fig_filename is not None:
        plt.savefig(fig_filename)
        print(f"stored figure as {fig_filename}")
    plt.show()


if __name__ == "__main__":
    main(mode='forw')
    main(mode='adj')
    main(mode='forw', one_sided=True)
    main(mode='adj', one_sided=True)
    main_value_gradient(one_sided=True)
