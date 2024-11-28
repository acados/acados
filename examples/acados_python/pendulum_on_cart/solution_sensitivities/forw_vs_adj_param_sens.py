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

def main(qp_solver_ric_alg: int, use_cython=False, generate_solvers=True, plot_trajectory=False):
    """
    Evaluate policy and calculate its gradient for the pendulum on a cart with a parametric model.
    """
    p_nominal = 1.0
    x0 = np.array([0.0, np.pi / 2, 0.0, 0.0])
    p_test = p_nominal + 0.2

    nx = len(x0)
    nu = 1

    N_horizon = 10
    T_horizon = 1.0

    # TODO: look more into why the test fails with these settings:
    # adjoint sensitivities of first parameter dont match with forward sensitivities
    # only "small" difference: in 5th digit, might be due to numerical inaccuracies
    # N_horizon = 50
    # T_horizon = 2.0
    # p_test = p_nominal + 0.3
    Fmax = 80.0

    cost_scale_as_param = True # test with 2 parameters
    with_parametric_constraint = True
    with_nonlinear_constraint = True

    ocp = export_parametric_ocp(x0=x0, N_horizon=N_horizon, T_horizon=T_horizon, Fmax=Fmax, qp_solver_ric_alg=1, cost_scale_as_param=cost_scale_as_param, with_parametric_constraint=with_parametric_constraint, with_nonlinear_constraint=with_nonlinear_constraint)
    if use_cython:
        raise NotImplementedError()
        AcadosOcpSolver.generate(ocp, json_file="parameter_augmented_acados_ocp.json")
        AcadosOcpSolver.build(ocp.code_export_directory, with_cython=True)
        ocp_solver = AcadosOcpSolver.create_cython_solver("parameter_augmented_acados_ocp.json")
    else:
        ocp_solver = AcadosOcpSolver(ocp, json_file="parameter_augmented_acados_ocp.json", generate=generate_solvers, build=generate_solvers)

    # create sensitivity solver
    ocp = export_parametric_ocp(x0=x0, N_horizon=N_horizon, T_horizon=T_horizon, Fmax=Fmax, hessian_approx='EXACT', qp_solver_ric_alg=qp_solver_ric_alg, cost_scale_as_param=cost_scale_as_param, with_parametric_constraint=with_parametric_constraint, with_nonlinear_constraint=with_nonlinear_constraint)
    ocp.model.name = 'sensitivity_solver'
    ocp.code_export_directory = f'c_generated_code_{ocp.model.name}'
    if use_cython:
        AcadosOcpSolver.generate(ocp, json_file=f"{ocp.model.name}.json")
        AcadosOcpSolver.build(ocp.code_export_directory, with_cython=True)
        sensitivity_solver = AcadosOcpSolver.create_cython_solver(f"{ocp.model.name}.json")
    else:
        sensitivity_solver = AcadosOcpSolver(ocp, json_file=f"{ocp.model.name}.json", generate=generate_solvers, build=generate_solvers)

    # set parameter value
    if cost_scale_as_param:
        p_val = np.array([p_test, 1.0])
    else:
        p_val = np.array([p_test])

    ocp_solver.set_p_global_and_precompute_dependencies(p_val)
    sensitivity_solver.set_p_global_and_precompute_dependencies(p_val)

    u_opt = ocp_solver.solve_for_x0(x0)[0]
    iterate = ocp_solver.store_iterate_to_obj()

    if with_parametric_constraint:
        lambdas = np.zeros((N_horizon-1, 2))
        for i in range(1, N_horizon):
            lam_ = ocp_solver.get(i, "lam")
            lambdas[i-1] = np.array([lam_[1], lam_[3]])
        print(f"lambdas of parametric constraints: {lambdas}\n")
        max_lam = np.max(np.abs(lambdas))
        print(f"max lambda of parametric constraints: {max_lam:.2f}\n")

    sensitivity_solver.load_iterate_from_obj(iterate)
    sensitivity_solver.solve_for_x0(x0, fail_on_nonzero_status=False, print_stats_on_failure=False)

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
            raise Exception("adj_p and adj_p_ref should match.")
            # print("ERROR: adj_p and adj_p_ref should match.")
        else:
            print("Success: adj_p and adj_p_ref match!")

    # test with list vs. single stage API: varying seed_x
    adj_p_ref = sensitivity_solver.eval_adjoint_solution_sensitivity(seed_x=None, seed_u=[(0, seed_ustage)])
    adj_p = sensitivity_solver.eval_adjoint_solution_sensitivity(seed_x=[], seed_u=[(0, seed_ustage)])
    adj_p_zero_x_seed = sensitivity_solver.eval_adjoint_solution_sensitivity(seed_x=[(1, 0*seed_xstage)], seed_u=[(0, seed_ustage)])

    if not np.allclose(adj_p, adj_p_ref, atol=TOL):
        raise Exception("adj_p and adj_p_ref should match.")
    if not np.allclose(adj_p, adj_p_zero_x_seed, atol=TOL):
        raise Exception("adj_p and adj_p_zero_x_seed should match.")
    print("Success: adj_p and adj_p_ref match! Tested with None and empty list for seed_x.")

    # test with list vs. single stage API: varying seed_u
    adj_p_ref = sensitivity_solver.eval_adjoint_solution_sensitivity(seed_x=[(0, seed_xstage)], seed_u=None)
    adj_p = sensitivity_solver.eval_adjoint_solution_sensitivity(seed_x=[(0, seed_xstage)], seed_u=[])
    adj_p_zero_x_seed = sensitivity_solver.eval_adjoint_solution_sensitivity(seed_x=[(0, seed_xstage)], seed_u=[(0, 0*seed_ustage)])

    if not np.allclose(adj_p, adj_p_ref, atol=1e-7):
        raise Exception("adj_p and adj_p_ref should match.")
    if not np.allclose(adj_p, adj_p_zero_x_seed, atol=1e-7):
        raise Exception("adj_p and adj_p_zero_x_seed should match.")
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
            raise Exception(f"adj_p_vec and adj_p_mat[{i}, :] should match.")
        else:
            print(f"Success: adj_p_vec and adj_p_mat[{i}, :] match!")

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


if __name__ == "__main__":
    main(qp_solver_ric_alg=0, use_cython=False, generate_solvers=True, plot_trajectory=False)
