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
from sensitivity_utils import plot_smoothed_solution_sensitivities_results, export_parametric_ocp, plot_pendulum


with_parametric_constraint = True
with_nonlinear_constraint = False

N_horizon = 50
T_horizon = 2.0
Fmax = 80.0


def solve_ocp_and_compute_sens(ocp_solver: AcadosOcpSolver, sensitivity_solver: AcadosOcpSolver, p_test, x0, tau_min, sanity_checks=True):

    ocp_solver.options_set('tau_min', tau_min)
    sensitivity_solver.options_set('tau_min', tau_min)

    np_test = p_test.shape[0]
    sens_u = np.zeros(np_test)
    u_opt = np.zeros(np_test)

    if with_parametric_constraint:
        n_lam_total = ocp_solver.get_flat('lam').shape[0]
        lambda_flat = np.zeros((np_test, n_lam_total))

    for i, p in enumerate(p_test):
        p_val = np.array([p])

        ocp_solver.set_p_global_and_precompute_dependencies(p_val)
        sensitivity_solver.set_p_global_and_precompute_dependencies(p_val)
        u_opt[i] = ocp_solver.solve_for_x0(x0, fail_on_nonzero_status=False)[0]
        status = ocp_solver.get_status()
        # ocp_solver.print_statistics()
        if status != 0:
            ocp_solver.print_statistics()
            print(f"Solver failed with status {status} for {i}th parameter value {p} and {tau_min=}.")
            breakpoint()

        iterate = ocp_solver.store_iterate_to_flat_obj()

        sensitivity_solver.load_iterate_from_flat_obj(iterate)
        sensitivity_solver.setup_qp_matrices_and_factorize()

        for j in range(1, N_horizon):
            # 1, 3 are indices of upper and lower multiplier for the parametric constraints
            lambda_flat[i, :] = ocp_solver.get_flat('lam')

        if ocp_solver.get_status() not in [0]:
            print(f"OCP solver returned status {ocp_solver.get_status()}.")
            breakpoint()
        if sensitivity_solver.get_status() not in [0, 2]:
            print(f"sensitivity solver returned status {sensitivity_solver.get_status()}.")
            # breakpoint()
        # Calculate the policy gradient
        out_dict = sensitivity_solver.eval_solution_sensitivity(0, "p_global", return_sens_x=False, sanity_checks=sanity_checks)
        sens_u[i] = out_dict['sens_u'].item()

    return u_opt, sens_u, lambda_flat

def create_solvers(x0, use_cython=False, qp_solver_ric_alg=0,
                    verbose = True, build = True, generate = True):
    ocp = export_parametric_ocp(x0=x0, N_horizon=N_horizon, T_horizon=T_horizon, Fmax=Fmax, qp_solver_ric_alg=1, with_parametric_constraint=with_parametric_constraint, with_nonlinear_constraint=with_nonlinear_constraint)

    # create nominal solver
    if use_cython:
        AcadosOcpSolver.generate(ocp, json_file="parameter_augmented_acados_ocp.json")
        AcadosOcpSolver.build(ocp.code_export_directory, with_cython=True)
        ocp_solver = AcadosOcpSolver.create_cython_solver("parameter_augmented_acados_ocp.json")
    else:
        ocp_solver = AcadosOcpSolver(ocp, build=build, generate=generate, json_file="parameter_augmented_acados_ocp.json", verbose=verbose)

    # create sensitivity solver
    ocp = export_parametric_ocp(x0=x0, N_horizon=N_horizon, T_horizon=T_horizon, Fmax=Fmax, hessian_approx='EXACT', qp_solver_ric_alg=qp_solver_ric_alg, with_parametric_constraint=with_parametric_constraint, with_nonlinear_constraint=with_nonlinear_constraint)
    # test with QP solver that does condensing: not recommended for sensitivtity solver
    ocp.solver_options.qp_solver_cond_N = int(N_horizon/4)

    ocp.model.name = 'sensitivity_solver'
    ocp.code_export_directory = f'c_generated_code_{ocp.model.name}'
    if use_cython:
        AcadosOcpSolver.generate(ocp, json_file=f"{ocp.model.name}.json")
        AcadosOcpSolver.build(ocp.code_export_directory, with_cython=True)
        sensitivity_solver = AcadosOcpSolver.create_cython_solver(f"{ocp.model.name}.json")
    else:
        sensitivity_solver = AcadosOcpSolver(ocp, build=build, generate=generate, json_file=f"{ocp.model.name}.json", verbose=verbose)

    return ocp_solver, sensitivity_solver


def main_parametric(qp_solver_ric_alg: int, use_cython=False, plot_trajectory=False):
    """
    Evaluate policy and calculate its gradient for the pendulum on a cart with a parametric model.
    """

    x0 = np.array([0.0, np.pi / 2, 0.0, 0.0])
    delta_p = 0.001
    # p_nominal = 1.0
    # p_test = np.arange(p_nominal + 0.1, p_nominal + 0.5, delta_p)
    p_test = np.arange(1.05, 1.4+delta_p, delta_p)

    ocp_solver, sensitivity_solver = create_solvers(x0, use_cython=use_cython, qp_solver_ric_alg=qp_solver_ric_alg,) # verbose=False, build=False, generate=False)
    ocp = ocp_solver.acados_ocp

    # compute policy and its gradient
    u_opt, sens_u, lambda_flat = solve_ocp_and_compute_sens(ocp_solver, sensitivity_solver, p_test, x0, tau_min=0.0)

    # Compare to numerical gradients
    sens_u_fd = np.gradient(u_opt, delta_p)
    u_opt_reconstructed_fd = np.cumsum(sens_u_fd) * delta_p + u_opt[0]
    u_opt_reconstructed_fd += u_opt[0] - u_opt_reconstructed_fd[0]

    u_opt_reconstructed_acados = np.cumsum(sens_u) * delta_p + u_opt[0]
    u_opt_reconstructed_acados += u_opt[0] - u_opt_reconstructed_acados[0]

    test_tol = 1e-2
    median_diff = np.median(np.abs(sens_u - sens_u_fd))
    print(f"Median difference between policy gradient obtained by acados and via FD is {median_diff} should be < {test_tol}.")
    # test: check median since derivative cannot be compared at active set changes
    assert median_diff <= test_tol

    # for multiplier plot
    n_lam_total = ocp_solver.get_flat('lam').shape[0]
    multipliers_bu = []
    multipliers_h = []
    nbu = ocp.dims.nbu
    nx = ocp.dims.nx
    x0_lam_idx = [*range(nbu, nx+nbu)] + [*range(2*nbu+nx, 2*nx+2*nbu)]
    n_lam_0 = ocp_solver.get(0, "lam").shape[0]
    bu_lam_idx = [*range(n_lam_0, n_lam_total, 2)]
    h_lam_idx = [*range(n_lam_0+1, n_lam_total, 2)]

    for i in range(n_lam_total):
        if np.max(np.abs(lambda_flat[:, i])) > 1e-2 and i not in x0_lam_idx:
            if i in bu_lam_idx:
                multipliers_bu += [lambda_flat[:, i]]
            elif i in h_lam_idx:
                multipliers_h += [lambda_flat[:, i]]
            else:
                print(f"found multiplier with index {i} that is not in x0_lam_idx, bu_lam_idx or h_lam_idx.")
            print(f"Multiplier {i} has non-zero values.")
    print(f"Multipliers with absolute value > 1e-2: bu {len(multipliers_bu)}, h {len(multipliers_h)}")

    # solutions to plot
    label = r'$\tau_{\mathrm{min}} = 0$'
    pi_label_pairs = []
    sens_pi_label_pairs = []

    pi_label_pairs.append((u_opt, label))
    sens_pi_label_pairs.append((sens_u, label))

    for tau_min in [1e-3, 1e-2]:
        u_opt, sens_u, _ = solve_ocp_and_compute_sens(ocp_solver, sensitivity_solver, p_test, x0, tau_min=tau_min)
        label = r'$\tau_{\mathrm{min}} = 10^{' + f"{int(np.log10(tau_min))}" + r"}$"
        pi_label_pairs.append((u_opt, label))
        sens_pi_label_pairs.append((sens_u, label))

    sens_pi_label_pairs.append((sens_u_fd, 'finite diff.'))

    # without 2-solver approach
    tau_min = 1e-6
    u_opt, sens_u, _ = solve_ocp_and_compute_sens(ocp_solver, ocp_solver, p_test, x0, tau_min=tau_min, sanity_checks=False)
    label = r"IFT approx. Hess."
    # pi_label_pairs.append((u_opt, label))
    sens_pi_label_pairs.append((sens_u, label))

    # plot
    # plot_smoothed_solution_sensitivities_results(p_test, pi_label_pairs, sens_pi_label_pairs, title=None, parameter_name=r"$\theta$",
    #              multipliers_bu=multipliers_bu, multipliers_h=multipliers_h,
    #              figsize=(7, 9),
    #              fig_filename="smoothed_solution_sensitivities.pdf",
    #              )
    plot_smoothed_solution_sensitivities_results(p_test, pi_label_pairs, sens_pi_label_pairs, title=None, parameter_name=r"$\theta$",
                multipliers_bu=multipliers_bu, multipliers_h=multipliers_h,
                figsize=(12, 3.2),
                fig_filename="smoothed_solution_sensitivities_horizontal.pdf",
                horizontal_plot=True,
                )
    if plot_trajectory:
        plot_pendulum_traj_from_ocp_iterate(ocp_solver)


def plot_pendulum_traj_from_ocp_iterate(ocp_solver: AcadosOcpSolver):
    ocp = ocp_solver.acados_ocp
    nx = ocp.dims.nx
    nu = ocp.dims.nu
    simX = np.zeros((N_horizon+1, nx))
    simU = np.zeros((N_horizon, nu))

    # get solution
    for i in range(N_horizon):
        simX[i,:] = ocp_solver.get(i, "x")
        simU[i,:] = ocp_solver.get(i, "u")
    simX[N_horizon,:] = ocp_solver.get(N_horizon, "x")

    plot_pendulum(ocp.solver_options.shooting_nodes, Fmax, simU, simX, latexify=True, time_label=ocp.model.t_label, x_labels=ocp.model.x_labels, u_labels=ocp.model.u_labels)


def main_plot_trajectories():
    x0 = np.array([0.0, np.pi / 2, 0.0, 0.0])
    ocp_solver, sensitivity_solver = create_solvers(x0, use_cython=False, qp_solver_ric_alg=0, verbose=False, build=False, generate=False)

    tau_min = 0.0
    ocp_solver.options_set('tau_min', tau_min)
    sensitivity_solver.options_set('tau_min', tau_min)

    p_test = np.array([1.4435, 1.444])
    np_test = p_test.shape[0]
    u_opt = np.zeros(np_test)

    for i, p in enumerate(p_test):
        p_val = np.array([p])

        ocp_solver.set_p_global_and_precompute_dependencies(p_val)
        sensitivity_solver.set_p_global_and_precompute_dependencies(p_val)
        u_opt[i] = ocp_solver.solve_for_x0(x0, fail_on_nonzero_status=False)[0]
        status = ocp_solver.get_status()
        ocp_solver.print_statistics()

        iterate = ocp_solver.store_iterate_to_flat_obj()
        sensitivity_solver.load_iterate_from_flat_obj(iterate)
        sensitivity_solver.setup_qp_matrices_and_factorize()
        diagnostics = sensitivity_solver.qp_diagnostics()
        print(diagnostics)

        if status != 0:
            ocp_solver.print_statistics()
            print(f"Solver failed with status {status} for {i}th parameter value {p} and {tau_min=}.")
            breakpoint()

        plot_pendulum_traj_from_ocp_iterate(ocp_solver)


if __name__ == "__main__":
    main_parametric(qp_solver_ric_alg=0, use_cython=False, plot_trajectory=True)
    # main_plot_trajectories()
