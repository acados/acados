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
from acados_template import AcadosOcp, AcadosModel, AcadosOcpSolver, latexify_plot
import matplotlib.pyplot as plt
latexify_plot()

P_SQUARED = False
if P_SQUARED:
    PROBLEM_NAME = "non_ocp_p_squared"
else:
    PROBLEM_NAME = "non_ocp_p_linear"

def export_parametric_nlp() -> AcadosOcp:

    model = AcadosModel()
    model.x = ca.SX.sym("x", 1)
    model.p_global = ca.SX.sym("p_global", 1)
    if P_SQUARED:
        model.cost_expr_ext_cost_e = (model.x - model.p_global**2)**2
    else:
        model.cost_expr_ext_cost_e = (model.x - model.p_global)**2
    model.name = "non_ocp"
    ocp = AcadosOcp()
    ocp.model = model

    ocp.constraints.lbx_e = np.array([-1.0])
    ocp.constraints.ubx_e = np.array([1.0])
    ocp.constraints.idxbx_e = np.array([0])

    ocp.cost.cost_type_e = "EXTERNAL"
    ocp.solver_options.qp_solver = "FULL_CONDENSING_HPIPM"
    ocp.solver_options.hessian_approx = "EXACT"
    ocp.solver_options.N_horizon = 0

    ocp.p_global_values = np.zeros((1,))
    ocp.solver_options.with_solution_sens_wrt_params = True
    ocp.solver_options.with_value_sens_wrt_params = True
    ocp.solver_options.nlp_solver_ext_qp_res = 1

    return ocp

def solve_and_compute_sens(p_test, tau):
    np_test = p_test.shape[0]

    ocp = export_parametric_nlp()
    ocp.solver_options.tau_min = tau
    ocp.solver_options.qp_solver_t0_init = 0
    ocp.solver_options.nlp_solver_max_iter = 2 # QP should converge in one iteration

    ocp_solver = AcadosOcpSolver(ocp, json_file="parameter_augmented_acados_ocp.json", verbose=False)

    sens_x = np.zeros(np_test)
    solution = np.zeros(np_test)

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
        out_dict = ocp_solver.eval_solution_sensitivity(0, "p_global", return_sens_x=True, return_sens_u=False)
        sens_x[i] = out_dict['sens_x'].item()

    return solution, sens_x

def main():
    p_nominal = 0.0
    delta_p = 0.002
    p_test = np.arange(p_nominal - 2, p_nominal + 2, delta_p)
    sens_list = []
    labels_list = []
    sol_list = []
    tau = 1e-6
    solution, sens_x = solve_and_compute_sens(p_test, tau)

    # Compare to numerical gradients
    sens_x_fd = np.gradient(solution, delta_p)
    test_tol = 1e-2
    median_diff = np.median(np.abs(sens_x - sens_x_fd))
    print(f"Median difference between policy gradient obtained by acados and via FD is {median_diff} should be < {test_tol}.")
    # test: check median since derivative cannot be compared at active set changes
    assert median_diff <= test_tol

    sens_list.append(sens_x)
    labels_list.append(r"$\tau = 10^{-6}$")
    sol_list.append(solution)

    tau_vals = [1e-4, 1e-3, 1e-2]
    for tau in tau_vals:
        sol_tau, sens_x_tau = solve_and_compute_sens(p_test, tau)
        sens_list.append(sens_x_tau)
        labels_list.append(r"$\tau = 10^{" + f"{int(np.log10(tau))}" + r"}$")
        # labels_list.append(r"$\tau =" + f"{tau}" + r"$")
        sol_list.append(sol_tau)

    plot_solution_sensitivities_results(p_test, sol_list, sens_list, labels_list,
                 title=None, parameter_name=r"$\theta$", fig_filename=f"solution_sens_{PROBLEM_NAME}.pdf")
    plot_solution_sensitivities_results(p_test, sol_list, sens_list, labels_list,
                 title=None, parameter_name=r"$\theta$", fig_filename=f"solution_sens_{PROBLEM_NAME}_transposed.pdf", horizontal_plot=True)

def plot_solution_sensitivities_results(p_test, sol_list, sens_list, labels_list, title=None, parameter_name="", fig_filename=None, horizontal_plot=False):
    p_min = p_test[0]
    p_max = p_test[-1]
    linestyles = ["--", "-.", "--", ":", "-.", ":"]

    nsub = 2
    if horizontal_plot:
        _, ax = plt.subplots(nrows=1, ncols=nsub, sharex=False, figsize=(12, 3.0))
    else:
        _, ax = plt.subplots(nrows=nsub, ncols=1, sharex=True, figsize=(6.5,5))

    isub = 0
    # plot analytic solution
    x_vals = np.linspace(-1, 1, 100)
    if P_SQUARED:
        ax[isub].plot([p_min, -1], [1, 1], "k-", linewidth=2, label="analytic")
        y_vals = x_vals**2
    else:
        ax[isub].plot([p_min, -1], [-1, -1], "k-", linewidth=2, label="analytic")
        y_vals = x_vals
    ax[isub].plot([1, p_max], [1, 1], "k-", linewidth=2)
    ax[isub].plot(x_vals, y_vals, "k-", linewidth=2)

    for i, sol in enumerate(sol_list):
        ax[isub].plot(p_test, sol, label=labels_list[i], linestyle=linestyles[i])
    ax[isub].set_xlim([p_test[0], p_test[-1]])
    ax[isub].set_ylabel(r"solution $x^{\star}$")
    if title is not None:
        ax[isub].set_title(title)
    ax[isub].legend()

    isub += 1

    # plot analytic sensitivity
    ax[isub].plot([p_min, -1], [0, 0], "k-", linewidth=2, label="analytic")
    ax[isub].plot([1, p_max], [0, 0], "k-", linewidth=2)
    if P_SQUARED:
        ax[isub].plot([-1, 1], [-2, 2], "k-", linewidth=2)
    else:
        ax[isub].plot([-1, 1], [1, 1], "k-", linewidth=2)

    # plot numerical sensitivities
    for i, sens_x_tau in enumerate(sens_list):
        ax[isub].plot(p_test, sens_x_tau, label=labels_list[i], color=f"C{i}", linestyle=linestyles[i])
    ax[isub].set_xlim([p_test[0], p_test[-1]])
    ax[isub].set_ylabel(r"derivative $\partial_\theta x^{\star}$")
    # ax[isub].legend(ncol=2)

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
    main()

    # to plot only analytic solution
    # plot_solution_sensitivities_results([-2, 2], [], [], [], parameter_name=r"$\theta$", fig_filename="solution_sens_non_ocp_analytic.pdf")
