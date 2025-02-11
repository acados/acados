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

def export_parametric_ocp() -> AcadosOcp:

    model = AcadosModel()
    model.x = ca.SX.sym("x", 1)
    model.p_global = ca.SX.sym("p_global", 1)
    model.disc_dyn_expr = model.x
    model.cost_expr_ext_cost = (model.x - model.p_global**2)**2
    model.cost_expr_ext_cost_e = 0
    model.name = "non_ocp"
    ocp = AcadosOcp()
    ocp.model = model

    ocp.constraints.lbx_0 = np.array([-1.0])
    ocp.constraints.ubx_0 = np.array([1.0])
    ocp.constraints.idxbx_0 = np.array([0])

    ocp.cost.cost_type = "EXTERNAL"
    ocp.solver_options.integrator_type = "DISCRETE"
    ocp.solver_options.hessian_approx = "EXACT"
    ocp.solver_options.N_horizon = 1
    ocp.solver_options.tf = 1.0

    ocp.p_global_values = np.zeros((1,))
    ocp.solver_options.with_solution_sens_wrt_params = True
    ocp.solver_options.with_value_sens_wrt_params = True
    ocp.solver_options.nlp_solver_ext_qp_res = 1

    return ocp

def solve_and_compute_sens(p_test, tau):
    np_test = p_test.shape[0]

    ocp = export_parametric_ocp()
    # ocp.solver_options.solution_sens_qp_t_lam_min = tau
    ocp.solver_options.qp_solver_tau_min = tau
    # ocp.solver_options.qp_solver_tol_comp = 1e-6
    # ocp.solver_options.nlp_solver_tol_comp = 1e1 * tau

    # solver creation arguments
    ocp_solver = AcadosOcpSolver(ocp, json_file="parameter_augmented_acados_ocp.json", verbose=False)

    sens_x = np.zeros(np_test)
    u_opt = np.zeros(np_test)

    for i, p in enumerate(p_test):
        p_val = np.array([p])

        ocp_solver.set_p_global_and_precompute_dependencies(p_val)
        status = ocp_solver.solve()
        u_opt[i] = ocp_solver.get(0, "x")

        if status not in [0]:
            print(f"OCP solver returned status {status} at {i}th p value {p}, {tau=}.")
            ocp_solver.print_statistics()
            breakpoint()

        ocp_solver.setup_qp_matrices_and_factorize()
        # Calculate the policy gradient
        out_dict = ocp_solver.eval_solution_sensitivity(0, "p_global", return_sens_x=True)
        sens_x[i] = out_dict['sens_x'].item()

    return u_opt, sens_x

def main():
    p_nominal = 0.0
    delta_p = 0.01
    p_test = np.arange(p_nominal - 2, p_nominal + 2, delta_p)
    sens_list = []
    labels_list = []
    # tau = 1e-6
    # u_opt, sens_x = solve_and_compute_sens(p_test, tau)

    # # Compare to numerical gradients
    # sens_x_fd = np.gradient(u_opt, delta_p)
    # test_tol = 1e-2
    # median_diff = np.median(np.abs(sens_x - sens_x_fd))
    # print(f"Median difference between policy gradient obtained by acados and via FD is {median_diff} should be < {test_tol}.")
    # # test: check median since derivative cannot be compared at active set changes
    # assert median_diff <= test_tol

    # sens_list = [sens_x]
    # labels_list = [r"$\tau = 10^{-6}$"]
    tau_vals = [1e-4, 1e-3, 1e-2]
    for tau in tau_vals:
        _, sens_x_tau = solve_and_compute_sens(p_test, tau)
        sens_list.append(sens_x_tau)
        labels_list.append(r"$\tau =" + f"{tau} $")

    plot_solution_sensitivities_results(p_test, u_opt, sens_list, labels_list,
                 title=None, parameter_name="p")


def plot_solution_sensitivities_results(p_test, pi, sens_list, labels_list, title=None, parameter_name=""):

    nsub = 2

    _, ax = plt.subplots(nrows=nsub, ncols=1, sharex=True, figsize=(9,9))

    isub = 0
    ax[isub].plot(p_test, pi, label='acados', color='k')
    ax[isub].set_ylabel(r"$u$")
    if title is not None:
        ax[isub].set_title(title)

    isub += 1
    linestyles = ["-", "--", "-.", ":", "-"]
    for i, sens_x_tau in enumerate(sens_list):
        ax[isub].plot(p_test, sens_x_tau, label=labels_list[i], color=f"C{i}", linestyle=linestyles[i])
    ax[isub].set_xlim([p_test[0], p_test[-1]])
    ax[isub].set_ylabel(r"$\partial_p u$")

    # isub += 1
    # ax[isub].plot(p_test, np.abs(sens_x- np_grad), "--", label='acados - finite diff')
    # ax[isub].set_ylabel(r"diff $\partial_p u$")
    # ax[isub].set_yscale("log")

    for i in range(isub+1):
        ax[i].grid()
        ax[i].legend()

    ax[-1].set_xlabel(f"{parameter_name}")

    fig_filename = f"solution_sens_{title}.pdf"
    plt.savefig(fig_filename)
    print(f"stored figure as {fig_filename}")
    plt.show()

if __name__ == "__main__":
    main()
