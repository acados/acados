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
from sensitivity_utils import export_parametric_ocp, plot_pendulum

TOL = 1e-6

def main(qp_solver_ric_alg: int, generate_solvers=True, plot_trajectory=False):
    """
    Evaluate policy and calculate its gradient for the pendulum on a cart with a parametric model.
    """
    p_nominal = 1.0
    x0 = np.array([0.0, np.pi / 2, 0.0, 0.0])

    nx = len(x0)

    N_horizon = 20
    T_horizon = 2.0
    p_test = p_nominal + 0.3
    Fmax = 80.0

    cost_scale_as_param = True
    with_parametric_constraint = True
    with_nonlinear_constraint = True

    ocp = export_parametric_ocp(x0=x0, N_horizon=N_horizon, T_horizon=T_horizon, Fmax=Fmax, qp_solver_ric_alg=1, cost_scale_as_param=cost_scale_as_param, with_parametric_constraint=with_parametric_constraint, with_nonlinear_constraint=with_nonlinear_constraint)
    ocp_solver = AcadosOcpSolver(ocp, json_file="parameter_augmented_acados_ocp.json", generate=generate_solvers, build=generate_solvers)

    # create sensitivity solver
    ocp = export_parametric_ocp(x0=x0, N_horizon=N_horizon, T_horizon=T_horizon, Fmax=Fmax, hessian_approx='EXACT', qp_solver_ric_alg=qp_solver_ric_alg, cost_scale_as_param=cost_scale_as_param, with_parametric_constraint=with_parametric_constraint, with_nonlinear_constraint=with_nonlinear_constraint)
    ocp.model.name = 'sensitivity_solver'
    ocp.code_export_directory = f'c_generated_code_{ocp.model.name}'
    sensitivity_solver = AcadosOcpSolver(ocp, json_file=f"{ocp.model.name}.json", generate=generate_solvers, build=generate_solvers)

    if cost_scale_as_param:
        p_vals = [
            np.array([p_nominal, 1.0]),
            np.array([p_nominal + 0.2, 1.0]),
            np.array([p_nominal + 0.3, 1.0]),
            np.array([p_nominal + 0.4, 1.0]),
        ]
    else:
        p_vals = np.array([p_test])

    for p_val in p_vals:
        print(f"Testing with p_val = {p_val}")

        # solve and compare forward and adjoint
        solve_and_compare_fwd_and_adj(ocp_solver, sensitivity_solver, x0, p_val, with_parametric_constraint)

    if plot_trajectory:
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




def solve_and_compare_fwd_and_adj(ocp_solver: AcadosOcpSolver,
                                  sensitivity_solver: AcadosOcpSolver,
                                  x0: np.ndarray,
                                  p_val: np.ndarray,
                                  with_parametric_constraint: bool):

    N_horizon = ocp_solver.N
    nx = ocp_solver.acados_ocp.dims.nx
    nu = ocp_solver.acados_ocp.dims.nu

    # set parameter value
    ocp_solver.set_p_global_and_precompute_dependencies(p_val)
    sensitivity_solver.set_p_global_and_precompute_dependencies(p_val)

    # nominal solve
    u_opt = ocp_solver.solve_for_x0(x0)[0]
    tau_iter = ocp_solver.get_stats("qp_tau_iter")
    print(f"qp tau iter: {tau_iter}\n")

    if with_parametric_constraint:
        lambdas = np.zeros((N_horizon-1, 2))
        for i in range(1, N_horizon):
            lam_ = ocp_solver.get(i, "lam")
            lambdas[i-1] = np.array([lam_[1], lam_[3]])
        # print(f"lambdas of parametric constraints: {lambdas}\n")
        max_lam = np.max(np.abs(lambdas))
        print(f"max lambda of parametric constraints: {max_lam:.2f}\n")

    # transfer iterate to sensitivity solver
    iterate = ocp_solver.store_iterate_to_obj()
    sensitivity_solver.load_iterate_from_obj(iterate)

    # setup QP matrices and factorize
    sensitivity_solver.setup_qp_matrices_and_factorize()

    qp_iter = sensitivity_solver.get_stats("qp_iter")
    print(f"qp iter: {qp_iter}\n")

    # check factorization
    # for i in range(1, N_horizon-1):
    #     P_mat = sensitivity_solver.get_from_qp_in(i, "P")
    #     K_mat = sensitivity_solver.get_from_qp_in(i, "K")
    #     Lr_mat = sensitivity_solver.get_from_qp_in(i, "Lr")
    #     print(f"stage {i} got factorization")
    #     # print(f"P_mat = {P_mat}")
    #     # print(f"K_mat = {K_mat}")
    #     print(f"Lr_mat = {Lr_mat}")

    if sensitivity_solver.get_status() not in [0, 2]:
        breakpoint()

    # adjoint direction for one stage
    seed_xstage = np.ones((nx, 1))
    seed_xstage[0, 0] = -8
    seed_ustage = np.ones((nu, 1))
    seed_ustage[0, 0] = 42

    # test different settings forward vs. adjoint
    for stages, seed_x_list, seed_u_list in [
        ([0, 1], [seed_xstage, seed_xstage], [seed_ustage, seed_ustage]),
        ([0, 3], [seed_xstage, seed_xstage], [seed_ustage, seed_ustage]),
        ([0], [seed_xstage], [seed_ustage]),
        ([5], [seed_xstage], [seed_ustage]),
    ]:
        # Calculate the policy gradient
        out_dict = sensitivity_solver.eval_solution_sensitivity(stages, "p_global")
        sens_x_forw = out_dict['sens_x']
        sens_u_forw = out_dict['sens_u']

        adj_p_ref = sum([seed_x_list[k].T @ sens_x_forw[k] + seed_u_list[k].T @ sens_u_forw[k] for k in range(len(stages))])

        # compute adjoint solution sensitivity
        adj_p = sensitivity_solver.eval_adjoint_solution_sensitivity(
                                                    seed_x=list(zip(stages, seed_x_list)),
                                                    seed_u=list(zip(stages, seed_u_list)))

        print(f"{adj_p=} {adj_p_ref=}")
        if not np.allclose(adj_p, adj_p_ref, atol=TOL):
            test_failure_message("adj_p and adj_p_ref should match.")
        else:
            print("Success: adj_p and adj_p_ref match!")

    # test with list vs. single stage API: varying seed_x
    adj_p_ref = sensitivity_solver.eval_adjoint_solution_sensitivity(seed_x=None, seed_u=[(0, seed_ustage)])
    adj_p = sensitivity_solver.eval_adjoint_solution_sensitivity(seed_x=[], seed_u=[(0, seed_ustage)])
    adj_p_zero_x_seed = sensitivity_solver.eval_adjoint_solution_sensitivity(seed_x=[(1, 0*seed_xstage)], seed_u=[(0, seed_ustage)])

    if not np.allclose(adj_p, adj_p_ref, atol=TOL):
        test_failure_message("adj_p and adj_p_ref should match.")
    if not np.allclose(adj_p, adj_p_zero_x_seed, atol=TOL):
        test_failure_message("adj_p and adj_p_zero_x_seed should match.")
    print("Success: adj_p and adj_p_ref match! Tested with None and empty list for seed_x.")

    # test with list vs. single stage API: varying seed_u
    adj_p_ref = sensitivity_solver.eval_adjoint_solution_sensitivity(seed_x=[(0, seed_xstage)], seed_u=None)
    adj_p = sensitivity_solver.eval_adjoint_solution_sensitivity(seed_x=[(0, seed_xstage)], seed_u=[])
    adj_p_zero_x_seed = sensitivity_solver.eval_adjoint_solution_sensitivity(seed_x=[(0, seed_xstage)], seed_u=[(0, 0*seed_ustage)])

    if not np.allclose(adj_p, adj_p_ref, atol=1e-7):
        test_failure_message("adj_p and adj_p_ref should match.")
    if not np.allclose(adj_p, adj_p_zero_x_seed, atol=1e-7):
        test_failure_message("adj_p and adj_p_zero_x_seed should match.")
    print("Success: adj_p and adj_p_ref match! Tested with None and empty list for seed_u.")

    # test multiple adjoint seeds at once
    seed_x_mat = np.eye(nx)
    seed_u_mat = np.zeros((nu, nx))
    adj_p_mat = sensitivity_solver.eval_adjoint_solution_sensitivity(seed_x=[(1, seed_x_mat)],
                                            seed_u=[(1, seed_u_mat)])
    print(f"{adj_p_mat=}")

    for i in range(nx):
        adj_p_vec = sensitivity_solver.eval_adjoint_solution_sensitivity(seed_x=[(1, seed_x_mat[:, [i]])],
                                            seed_u=[(1, seed_u_mat[:, [i]])])
        print(f"{adj_p_vec=} {adj_p_mat[i, :]=}")
        if not np.allclose(adj_p_vec, adj_p_mat[i, :], atol=TOL):
            test_failure_message(f"adj_p_vec and adj_p_mat[{i}, :] should match.")
        else:
            print(f"Success: adj_p_vec and adj_p_mat[{i}, :] match!")


def test_failure_message(msg):
    # print(f"ERROR: {msg}")
    raise Exception(msg)

if __name__ == "__main__":
    main(qp_solver_ric_alg=0, generate_solvers=True, plot_trajectory=False)
