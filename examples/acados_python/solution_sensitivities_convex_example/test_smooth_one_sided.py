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
from non_ocp_example import plot_solution_sensitivities_results
from acados_template import AcadosOcpSolver, AcadosModel, AcadosOcp, ACADOS_INFTY


def create_parametric_nlp() -> AcadosOcp:

    model = AcadosModel()
    model.x = ca.SX.sym("x", 1)
    model.p_global = ca.SX.sym("p_global", 1)
    model.cost_expr_ext_cost_e = (model.x - model.p_global)**2
    model.name = "non_ocp"
    ocp = AcadosOcp()
    ocp.model = model

    ocp.constraints.lbx_e = np.array([-ACADOS_INFTY])
    ocp.constraints.ubx_e = np.array([1.0])
    ocp.constraints.idxbx_e = np.array([0])

    ocp.cost.cost_type_e = "EXTERNAL"
    ocp.solver_options.qp_solver = "FULL_CONDENSING_HPIPM"
    ocp.solver_options.qp_solver_ric_alg = 0
    ocp.solver_options.hessian_approx = "EXACT"
    ocp.solver_options.N_horizon = 0
    ocp.solver_options.print_level = 8

    ocp.p_global_values = np.zeros((1,))
    ocp.solver_options.with_solution_sens_wrt_params = True
    ocp.solver_options.with_value_sens_wrt_params = True
    ocp.solver_options.nlp_solver_ext_qp_res = 1
    ocp.solver_options.qp_solver_mu0 = 1e0
    # TODO: needed?
    ocp.solver_options.qp_tol = 1e-7

    return ocp


def solve_and_compute_sens(p_test, tau):
    np_test = p_test.shape[0]

    ocp = create_parametric_nlp()
    ocp.solver_options.tau_min = 0.01
    ocp.solver_options.qp_solver_t0_init = 0
    ocp.solver_options.nlp_solver_ext_qp_res = 1
    ocp.solver_options.nlp_solver_max_iter = 1 # QP should converge in one iteration
    
    ocp_solver = AcadosOcpSolver(ocp, json_file="parameter_augmented_acados_ocp.json", verbose=False)

    # p = 1.504999999999614
    # p_val = np.array([p])
    # ocp_solver.set_p_global_and_precompute_dependencies(p_val)
    # status = ocp_solver.solve()

    # breakpoint()

    sens_x = np.zeros(np_test)
    solution = np.zeros(np_test)

    for i, p in enumerate(p_test):
        p_val = np.array([p])

        ocp_solver.set_p_global_and_precompute_dependencies(p_val)
        status = ocp_solver.solve()
        solution[i] = ocp_solver.get(0, "x")[0]

        if status != 0:
            ocp_solver.print_statistics()
            # raise Exception(f"OCP solver returned status {status} at {i}th p value {p}, {tau=}.")
            print(f"OCP solver returned status {status} at {i}th p value {p}, {tau=}.")
            breakpoint()

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
    delta_p = 0.001
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
                 title=None, parameter_name=r"$\theta$", fig_filename=f"solution_sens_one_sided.pdf")
    

if __name__ == "__main__":
    main()