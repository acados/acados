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
from acados_template import AcadosOcpSolver
from sensitivity_utils import plot_results, export_parametric_ocp, evaluate_hessian_eigenvalues


def main_parametric(qp_solver_ric_alg: int, eigen_analysis=True, use_cython=False):
    """
    Evaluate policy and calculate its gradient for the pendulum on a cart with a parametric model.
    """

    if eigen_analysis and use_cython:
        raise Exception("Eigenvalue analysis is not possible with the cython interface.")

    p_nominal = 1.0
    x0 = np.array([0.0, np.pi / 2, 0.0, 0.0])
    delta_p = 0.002
    p_test = np.arange(p_nominal - 0.5, p_nominal + 0.5, delta_p)

    np_test = p_test.shape[0]
    N_horizon = 50
    T_horizon = 2.0
    Fmax = 80.0

    ocp = export_parametric_ocp(x0=x0, N_horizon=N_horizon, T_horizon=T_horizon, Fmax=Fmax, qp_solver_ric_alg=1)
    if use_cython:
        AcadosOcpSolver.generate(ocp, json_file="parameter_augmented_acados_ocp.json")
        AcadosOcpSolver.build(ocp.code_export_directory, with_cython=True)
        acados_ocp_solver = AcadosOcpSolver.create_cython_solver("parameter_augmented_acados_ocp.json")
    else:
        acados_ocp_solver = AcadosOcpSolver(ocp, json_file="parameter_augmented_acados_ocp.json")

    # create sensitivity solver
    ocp = export_parametric_ocp(x0=x0, N_horizon=N_horizon, T_horizon=T_horizon, Fmax=Fmax, hessian_approx='EXACT', qp_solver_ric_alg=qp_solver_ric_alg)
    ocp.model.name = 'sensitivity_solver'
    ocp.code_export_directory = f'c_generated_code_{ocp.model.name}'
    if use_cython:
        AcadosOcpSolver.generate(ocp, json_file=f"{ocp.model.name}.json")
        AcadosOcpSolver.build(ocp.code_export_directory, with_cython=True)
        sensitivity_solver = AcadosOcpSolver.create_cython_solver(f"{ocp.model.name}.json")
    else:
        sensitivity_solver = AcadosOcpSolver(ocp, json_file=f"{ocp.model.name}.json")

    if eigen_analysis:
        min_eig_full = np.zeros(np_test)
        min_abs_eig_full = np.zeros(np_test)
        min_abs_eig_proj_hess = np.zeros(np_test)
        min_eig_proj_hess = np.zeros(np_test)
        min_eig_P = np.zeros(np_test)
        min_abs_eig_P = np.zeros(np_test)
    else:
        min_eig_full = min_abs_eig_full = min_abs_eig_proj_hess = min_eig_proj_hess = min_eig_P = min_abs_eig_P = None

    sens_u = np.zeros(np_test)
    u_opt = np.zeros(np_test)
    for i, p in enumerate(p_test):
        p_val = np.array([p])

        for n in range(N_horizon+1):
            acados_ocp_solver.set(n, 'p', p_val)
            sensitivity_solver.set(n, 'p', p_val)
        u_opt[i] = acados_ocp_solver.solve_for_x0(x0)[0]

        acados_ocp_solver.store_iterate(filename='iterate.json', overwrite=True, verbose=False)

        sensitivity_solver.load_iterate(filename='iterate.json', verbose=False)
        sensitivity_solver.solve_for_x0(x0, fail_on_nonzero_status=False, print_stats_on_failure=False)
        # residuals = sensitivity_solver.get_stats("residuals")
        # print(f"residuals sensitivity_solver {residuals} status {sensitivity_solver.status}")

        if sensitivity_solver.get_status() not in [0, 2]:
            breakpoint()

        if eigen_analysis:
            min_eig_full[i], min_abs_eig_full[i], min_abs_eig_proj_hess[i], min_eig_proj_hess[i], min_eig_P[i], min_abs_eig_P[i] = evaluate_hessian_eigenvalues(sensitivity_solver, N_horizon)

        # Calculate the policy gradient
        _, sens_u_ = sensitivity_solver.eval_solution_sensitivity(0, "params_global")
        sens_u[i] = sens_u_.item()

    # Compare to numerical gradients
    sens_u_fd = np.gradient(u_opt, delta_p)
    u_opt_reconstructed_fd = np.cumsum(sens_u_fd) * delta_p + u_opt[0]
    u_opt_reconstructed_fd += u_opt[0] - u_opt_reconstructed_fd[0]

    u_opt_reconstructed_acados = np.cumsum(sens_u) * delta_p + u_opt[0]
    u_opt_reconstructed_acados += u_opt[0] - u_opt_reconstructed_acados[0]

    plot_results(p_test, u_opt, u_opt_reconstructed_acados, u_opt_reconstructed_fd, sens_u, sens_u_fd,
                 min_eig_full, min_eig_proj_hess, min_eig_P,
                 min_abs_eig_full, min_abs_eig_proj_hess, min_abs_eig_P,
                 eigen_analysis, qp_solver_ric_alg, parameter_name="mass")

    # test: check median since derivative cannot be compared at active set changes
    assert np.median(np.abs(sens_u - sens_u_fd)) <= 1e-2


if __name__ == "__main__":
    main_parametric(qp_solver_ric_alg=0, eigen_analysis=False, use_cython=True)
