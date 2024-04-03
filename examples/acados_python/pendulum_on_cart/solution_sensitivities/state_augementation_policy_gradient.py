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

import sys
sys.path.insert(0, '../common')

from acados_template import AcadosOcp, AcadosOcpSolver, AcadosModel
import casadi as ca
import numpy as np
import scipy.linalg as scipylinalg
from sensitivity_utils import evaluate_hessian_eigenvalues, plot_results

"""
    This example computes solution sensitivities with respect to parameters using state augmentation.
    If the parameters enter only via the cost and the dynamics, this approach is not recommended.
    Instead use the solution sensitivities provided by the acados OCP solver, cf. `policy_gradient_example.py`.
"""

def export_parameter_augmented_pendulum_ode_model(param_M_as_state=True) -> AcadosModel:
    """
    Augment the normal state vector with cart mass M.

    Return:
        AcadosModel: model of the augmented state cart
        nx_original: number of states before augmentating model with parameters
    """
    model_name = "parameter_augmented_pendulum_ode"

    # constants
    m = 0.1  # mass of the ball [kg]
    g = 9.81  # gravity constant [m/s^2]
    l = 0.8  # length of the rod [m]

    # set up states
    p = ca.SX.sym("p")
    theta = ca.SX.sym("theta")
    v = ca.SX.sym("v")
    omega = ca.SX.sym("omega")

    # controls
    F = ca.SX.sym("F")
    u = ca.vertcat(F)

    # xdot
    p_dot = ca.SX.sym("p_dot")
    theta_dot = ca.SX.sym("theta_dot")
    v_dot = ca.SX.sym("v_dot")
    omega_dot = ca.SX.sym("omega_dot")
    M_dot = ca.SX.sym("M_dot")

    nx_original = 4
    if param_M_as_state:
        M = ca.SX.sym("M")  # mass of the cart [kg]
        x = ca.vertcat(p, theta, v, omega, M)
        xdot = ca.vertcat(p_dot, theta_dot, v_dot, omega_dot, M_dot)
    else:
        M = 1.0  # mass of the cart [kg]
        x = ca.vertcat(p, theta, v, omega)
        xdot = ca.vertcat(p_dot, theta_dot, v_dot, omega_dot)

    # dynamics
    cos_theta = ca.cos(theta)
    sin_theta = ca.sin(theta)
    denominator = M + m - m * cos_theta * cos_theta
    f_expl = ca.vertcat(
        v,
        omega,
        (-m * l * sin_theta * omega * omega + m * g * cos_theta * sin_theta + F) / denominator,
        (-m * l * cos_theta * sin_theta * omega * omega + F * cos_theta + (M + m) * g * sin_theta) / (l * denominator),
    )

    if param_M_as_state:
        f_expl = ca.vertcat(f_expl, 0)

    model = AcadosModel()

    model.f_impl_expr = xdot - f_expl
    model.f_expl_expr = f_expl
    model.x = x
    model.xdot = xdot
    model.u = u
    model.name = model_name

    return model, nx_original


def export_parameter_augmented_ocp(
    x0=np.array([0.0, np.pi / 6, 0.0, 0.0, 1.0]), N_horizon=50, T_horizon=2.0, Fmax=80.0,
    hessian_approx = "GAUSS_NEWTON", param_M_as_state=True, qp_solver_ric_alg=1
) -> AcadosOcp:
    """
    OCP with augmented state vector (p, theta, v, omega, M).
    """
    ocp = AcadosOcp()
    ocp.model, nx_original = export_parameter_augmented_pendulum_ode_model(param_M_as_state=param_M_as_state)

    # set dimensions
    ocp.dims.N = N_horizon
    nu = ocp.model.u.rows()
    nx = ocp.model.x.rows()

    # set cost
    Q_mat = 2 * np.diag([1e3, 1e3, 1e-2, 1e-2])
    R_mat = 2 * np.diag([1e-1])

    ocp.cost.cost_type = "NONLINEAR_LS"
    ocp.cost.cost_type_e = "NONLINEAR_LS"
    ocp.cost.W = scipylinalg.block_diag(Q_mat, R_mat)
    ocp.cost.W_e = Q_mat

    ocp.model.cost_y_expr = ca.vertcat(ocp.model.x[:nx_original], ocp.model.u)
    ocp.model.cost_y_expr_e = ocp.model.x[:nx_original]
    ocp.cost.yref = np.zeros((nx_original + nu,))
    ocp.cost.yref_e = np.zeros((nx_original,))

    # set constraints
    ocp.constraints.lbu = np.array([-Fmax])
    ocp.constraints.ubu = np.array([+Fmax])
    ocp.constraints.idxbu = np.array([0])

    ocp.constraints.x0 = x0[:nx]

    # set options
    ocp.solver_options.qp_solver = "PARTIAL_CONDENSING_HPIPM"  # FULL_CONDENSING_QPOASES
    ocp.solver_options.integrator_type = "IRK"  # "DISCRETE"
    # ocp.solver_options.print_level = 1
    ocp.solver_options.nlp_solver_type = "SQP"  # SQP_RTI, SQP
    ocp.solver_options.qp_solver_cond_N = N_horizon

    # set prediction horizon
    ocp.solver_options.tf = T_horizon

    ocp.solver_options.qp_solver_ric_alg = qp_solver_ric_alg
    ocp.solver_options.hessian_approx = hessian_approx  # 'GAUSS_NEWTON', 'EXACT'
    if hessian_approx == 'EXACT':
        ocp.solver_options.nlp_solver_step_length = 0.0
        ocp.solver_options.nlp_solver_max_iter = 1
        ocp.solver_options.qp_solver_iter_max = 200
        ocp.solver_options.tol = 1e-10
    else:
        ocp.solver_options.nlp_solver_max_iter = 400
        ocp.solver_options.tol = 1e-8

    return ocp



def main_augmented(param_M_as_state: bool, idxp: int, qp_solver_ric_alg: int, eigen_analysis=True):
    """
    Evaluate policy and calculate its gradient for the pendulum on a cart with an augmented state formulation for
    varying M.
    """

    p_nominal = 1.0
    x0 = np.array([0.0, np.pi / 2, 0.0, 0.0, p_nominal])
    if idxp == 4:
        # vary M
        parameter_name = 'M'
        delta_p = 0.005
        p_test = np.arange(p_nominal - 0.5, p_nominal + 0.5, delta_p)
        x0_augmented = [np.array(x0[:-1].tolist() + [p]) for p in p_test]

    elif idxp == 1:
        # vary theta instead
        parameter_name = r'$\theta$'
        delta_p = np.pi/100
        p_test = np.arange(0.0, np.pi/5, delta_p)
        x0_augmented = [np.array([x0[0]] + [p] + x0[2:].tolist()) for p in p_test]

    if not param_M_as_state:
        x0_augmented = [x[:-1] for x in x0_augmented]

    if idxp == 4 and not param_M_as_state:
        raise Exception(f"can only get sensitivities wrt {idxp=} if param_M_as_state==True.")

    np_test = p_test.shape[0]
    N_horizon = 50
    T_horizon = 2.0
    Fmax = 80.0

    ocp = export_parameter_augmented_ocp(x0=x0, N_horizon=N_horizon, T_horizon=T_horizon, Fmax=Fmax, param_M_as_state=param_M_as_state, qp_solver_ric_alg=1)
    acados_ocp_solver = AcadosOcpSolver(ocp, json_file="parameter_augmented_acados_ocp.json")

    # create sensitivity solver
    ocp = export_parameter_augmented_ocp(x0=x0, N_horizon=N_horizon, T_horizon=T_horizon, Fmax=Fmax, hessian_approx='EXACT', param_M_as_state=param_M_as_state, qp_solver_ric_alg=qp_solver_ric_alg)
    ocp.model.name = 'sensitivity_solver'
    sensitivity_solver = AcadosOcpSolver(ocp, json_file="sensitivity_solver.json")

    sens_u = np.zeros(np_test)

    if eigen_analysis:
        min_eig_full = np.zeros(np_test)
        min_abs_eig_full = np.zeros(np_test)
        min_abs_eig_proj_hess = np.zeros(np_test)
        min_eig_proj_hess = np.zeros(np_test)
        min_eig_P = np.zeros(np_test)
        min_abs_eig_P = np.zeros(np_test)
    else:
        min_eig_full = None
        min_abs_eig_full = None
        min_abs_eig_proj_hess = None
        min_eig_proj_hess = None
        min_eig_P = None
        min_abs_eig_P = None

    u_opt = np.zeros(np_test)
    for i, x in enumerate(x0_augmented):
        # Evaluate the policy
        u_opt[i] = acados_ocp_solver.solve_for_x0(x).item()
        # acados_ocp_solver.print_statistics()
        acados_ocp_solver.store_iterate(filename='iterate.json', overwrite=True, verbose=False)

        sensitivity_solver.load_iterate(filename='iterate.json', verbose=False)
        sensitivity_solver.solve_for_x0(x, fail_on_nonzero_status=False, print_stats_on_failure=False)
        # residuals = sensitivity_solver.get_stats("residuals")
        # print(f"residuals sensitivity_solver {residuals} status {sensitivity_solver.status}")

        if sensitivity_solver.status not in [0, 2]:
            print(f"warning")
            breakpoint()

        if eigen_analysis:
            min_eig_full[i], min_abs_eig_full[i], min_abs_eig_proj_hess[i], min_eig_proj_hess[i], min_eig_P[i], min_abs_eig_P[i] = evaluate_hessian_eigenvalues(sensitivity_solver)

        # Calculate the policy gradient
        _, sens_u_ = sensitivity_solver.eval_solution_sensitivity(0, "initial_state")
        sens_u[i] = sens_u_[:, idxp]

    # Compare to numerical gradients
    np_grad = np.gradient(u_opt, delta_p)

    u_opt_reconstructed_np_grad = np.cumsum(np_grad) * delta_p + u_opt[0]
    u_opt_reconstructed_np_grad += u_opt[0] - u_opt_reconstructed_np_grad[0]

    u_opt_reconstructed_acados = np.cumsum(sens_u) * delta_p + u_opt[0]
    u_opt_reconstructed_acados += u_opt[0] - u_opt_reconstructed_acados[0]

    plot_results(p_test, u_opt, u_opt_reconstructed_acados, u_opt_reconstructed_np_grad, sens_u, np_grad,
                 min_eig_full, min_eig_proj_hess, min_eig_P,
                 min_abs_eig_full, min_abs_eig_proj_hess, min_abs_eig_P,
                 eigen_analysis, qp_solver_ric_alg, parameter_name)


if __name__ == "__main__":
    main_augmented(param_M_as_state=True, idxp=4, qp_solver_ric_alg=0, eigen_analysis=False)
