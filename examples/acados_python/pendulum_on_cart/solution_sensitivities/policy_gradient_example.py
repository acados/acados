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

from acados_template import AcadosOcp, AcadosOcpSolver, AcadosModel, latexify_plot
from casadi import SX, vertcat, cos, sin
import numpy as np
import scipy.linalg as scipylinalg
import matplotlib.pyplot as plt
from utils import plot_pendulum



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
    p = SX.sym("p")
    theta = SX.sym("theta")
    v = SX.sym("v")
    omega = SX.sym("omega")

    # controls
    F = SX.sym("F")
    u = vertcat(F)

    # xdot
    p_dot = SX.sym("p_dot")
    theta_dot = SX.sym("theta_dot")
    v_dot = SX.sym("v_dot")
    omega_dot = SX.sym("omega_dot")
    M_dot = SX.sym("M_dot")

    nx_original = 4
    if param_M_as_state:
        M = SX.sym("M")  # mass of the cart [kg]
        x = vertcat(p, theta, v, omega, M)
        xdot = vertcat(p_dot, theta_dot, v_dot, omega_dot, M_dot)
    else:
        M = 1.0  # mass of the cart [kg]
        x = vertcat(p, theta, v, omega)
        xdot = vertcat(p_dot, theta_dot, v_dot, omega_dot)

    # dynamics
    cos_theta = cos(theta)
    sin_theta = sin(theta)
    denominator = M + m - m * cos_theta * cos_theta
    f_expl = vertcat(
        v,
        omega,
        (-m * l * sin_theta * omega * omega + m * g * cos_theta * sin_theta + F) / denominator,
        (-m * l * cos_theta * sin_theta * omega * omega + F * cos_theta + (M + m) * g * sin_theta) / (l * denominator),
    )

    if param_M_as_state:
        f_expl = vertcat(f_expl, 0)

    # algebraic variables
    # z = None

    # parameters
    param = []

    f_impl = xdot - f_expl

    model = AcadosModel()

    model.f_impl_expr = f_impl
    model.f_expl_expr = f_expl
    model.x = x
    model.xdot = xdot
    model.u = u
    # model.z = z
    model.p = param
    model.name = model_name

    return model, nx_original


def export_parameter_augmented_ocp(
    x0=np.array([0.0, np.pi / 6, 0.0, 0.0, 1.0]), N_horizon=50, T_horizon=2.0, Fmax=80.0,
    hessian_approx = "GAUSS_NEWTON", param_M_as_state=True, qp_solver_ric_alg=1
) -> AcadosOcp:
    """
    OCP with augmented state vector (p, theta, v, omega, M).
    """
    # create ocp object to formulate the OCP
    ocp = AcadosOcp()

    # set model
    ocp.model, nx_original = export_parameter_augmented_pendulum_ode_model(param_M_as_state=param_M_as_state)

    # set dimensions
    ocp.dims.N = N_horizon
    nu = ocp.model.u.size()[0]
    nx = ocp.model.x.size()[0]

    # set cost
    Q_mat = 2 * np.diag([1e3, 1e3, 1e-2, 1e-2])
    R_mat = 2 * np.diag([1e-1])

    # We use the NONLINEAR_LS cost type and GAUSS_NEWTON Hessian approximation - One can also use the external cost module to specify generic cost.
    ocp.cost.cost_type = "NONLINEAR_LS"
    ocp.cost.cost_type_e = "NONLINEAR_LS"
    ocp.cost.W = scipylinalg.block_diag(Q_mat, R_mat)
    ocp.cost.W_e = Q_mat

    ocp.model.cost_y_expr = vertcat(ocp.model.x[:nx_original], ocp.model.u)
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
    # PARTIAL_CONDENSING_HPIPM, FULL_CONDENSING_QPOASES, FULL_CONDENSING_HPIPM,
    # PARTIAL_CONDENSING_QPDUNES, PARTIAL_CONDENSING_OSQP, FULL_CONDENSING_DAQP
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
        # ocp.solver_options.hpipm_mode = 'SPEED_ABS'
    else:
        ocp.solver_options.nlp_solver_max_iter = 400
        ocp.solver_options.tol = 1e-8

    return ocp




def evaluate_hessian_eigenvalues(acados_solver: AcadosOcpSolver):
    offset = 0
    min_eigv_total = 1e12
    min_abs_eigv = 1e12
    N_horizon = acados_solver.acados_ocp.dims.N

    for i in range(N_horizon+1):
        hess_block_acados = acados_solver.get_hessian_block(i)
        nv = hess_block_acados.shape[0]
        offset += nv

        eigv = np.linalg.eigvals(hess_block_acados)
        min_eigv = np.min(eigv)
        min_eigv_total = min(min_eigv, min_eigv_total)
        min_abs_eigv = min(min_abs_eigv, np.min(np.abs(eigv)))

    # check projected Hessian
    min_abs_eig_proj_hess = 1e12
    min_eig_proj_hess = 1e12
    min_eig_P = 1e12
    min_abs_eig_P = 1e12
    for i in range(1, N_horizon):
        P_mat = acados_solver.get_from_qp_in(i, 'P')
        B_mat = acados_solver.get_from_qp_in(i-1, 'B')
        # Lr: lower triangular decomposition of R within Riccati != R in qp_in!
        Lr = acados_solver.get_from_qp_in(i-1, 'Lr')
        R_ric = Lr @ Lr.T
        proj_hess_block = R_ric + B_mat.T @ P_mat @ B_mat
        eigv = np.linalg.eigvals(proj_hess_block)
        min_eigv = np.min(eigv)
        min_eig_proj_hess = min(min_eigv, min_eig_proj_hess)
        # if min_eig_proj_hess < 0:
        #     print(f"min_eig_proj_hess < 0 at {i}")
        #     breakpoint()
        min_abs_eig_proj_hess = min(min_abs_eig_proj_hess, np.min(np.abs(eigv)))
        # P
        eigv = np.linalg.eigvals(P_mat)
        min_eig_P = min(min_eig_P, np.min(eigv))
        min_abs_eig_P = min(min_abs_eig_P, np.min(np.abs(eigv)))

    return min_eigv_total, min_abs_eigv, min_abs_eig_proj_hess, min_eig_proj_hess, min_eig_P, min_abs_eig_P


def main(param_M_as_state: bool, idxp: int, qp_solver_ric_alg: int, eigen_analysis=True):
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
    N_horizon=50
    T_horizon=2.0
    Fmax=80.0

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

    pi = np.zeros(np_test)
    for i, x in enumerate(x0_augmented):
        # Evaluate the policy
        pi[i] = acados_ocp_solver.solve_for_x0(x)[0]
        # acados_ocp_solver.print_statistics()
        acados_ocp_solver.store_iterate(filename='iterate.json', overwrite=True, verbose=False)

        sensitivity_solver.load_iterate(filename='iterate.json', verbose=False)
        sensitivity_solver.solve_for_x0(x, fail_on_nonzero_status=False, print_stats_on_failure=False)
        residuals = sensitivity_solver.get_stats("residuals")
        print(f"residuals sensitivity_solver {residuals} status {sensitivity_solver.status}")
        if sensitivity_solver.status not in [0, 2]:
            print(f"warning")
            breakpoint()

        if eigen_analysis:
            min_eig_full[i], min_abs_eig_full[i], min_abs_eig_proj_hess[i], min_eig_proj_hess[i], min_eig_P[i], min_abs_eig_P[i] = evaluate_hessian_eigenvalues(sensitivity_solver)

        # Calculate the policy gradient
        sensitivity_solver.eval_param_sens(idxp)
        sens_u[i] = sensitivity_solver.get(0, "sens_u")[0]

        # plot solution
        # if i < 0:
        #     nx = max(x0.shape)
        #     nu = 1
        #     simX = np.ndarray((N_horizon+1, nx))
        #     simU = np.ndarray((N_horizon, nu))
        #     for i in range(N_horizon):
        #         simX[i, :] = acados_ocp_solver.get(i, "x")
        #         simU[i, :] = acados_ocp_solver.get(i, "u")
        #     simX[N_horizon, :] = acados_ocp_solver.get(N_horizon, "x")
        #     plot_pendulum(np.linspace(0, T_horizon, N_horizon+1), Fmax, simU, simX, states_lables = ['$x$', r'$\theta$', '$v$', r'$\dot{\theta}$', 'mass'])

    # Compare to numerical gradients
    np_grad = np.gradient(pi, delta_p)
    central_diff = (pi[2:] - pi[:-2]) / (p_test[2:] - p_test[:-2])

    pi_reconstructed_np_grad = np.cumsum(np_grad) * delta_p + pi[0]
    pi_reconstructed_np_grad += pi[0] - pi_reconstructed_np_grad[0]

    latexify_plot()
    nsub = 3
    if eigen_analysis:
        nsub +=2

    _, ax = plt.subplots(nrows=nsub, ncols=1, sharex=True, figsize=(9,9))
    isub = 0
    ax[isub].plot(p_test, pi, label='$p$ acados')
    ax[isub].plot(p_test, pi_reconstructed_np_grad, "--", label='$p$ reconstructed from finite diff')
    ax[isub].set_ylabel(r"$u$")
    ax[isub].legend()
    ax[isub].set_title(f'qp_solver_ric_alg {qp_solver_ric_alg}')
    ax[isub].grid(True)

    isub += 1
    ax[isub].plot(p_test, sens_u, label="acados")
    ax[isub].plot(p_test, np_grad, "--", label="np.grad")
    ax[isub].legend()
    ax[isub].set_xlim([p_test[0], p_test[-1]])
    ax[isub].grid(True)
    ax[isub].set_ylabel(r"$\partial_p u$")

    isub += 1
    # ax[isub].plot(p_test, np.abs(sens_u))
    # ax[isub].plot(p_test, np.abs(np_grad), "--")
    # ax[isub].plot(p_test[1:-1], np.abs(central_diff), "-.")
    ax[isub].plot(p_test, np.abs(sens_u- np_grad), "--", label='acados - finite diff')
    ax[isub].legend()
    ax[isub].set_ylabel(r"diff $\partial_p u$")
    ax[isub].set_yscale("log")
    ax[isub].grid(True)

    if eigen_analysis:
        isub += 1
        ax[isub].plot(p_test, np.sign(min_eig_full), linestyle="-", alpha=.6, label='full Hessian')
        ax[isub].plot(p_test, np.sign(min_eig_proj_hess), "--", label='proj Hessian')
        ax[isub].plot(p_test, np.sign(min_eig_P), ":", label='$P$ Riccati')
        ax[isub].grid(True)
        ax[isub].legend()
        ax[isub].set_ylabel("sign min eig")

        isub += 1
        ax[isub].plot(p_test, min_abs_eig_full, "--", label='full Hessian')
        ax[isub].plot(p_test, min_abs_eig_proj_hess, "--", label='proj Hessian')
        ax[isub].plot(p_test, min_abs_eig_P, "--", label='$P$ Riccati')
        ax[isub].set_ylabel("min abs eig")
        ax[isub].grid(True)
        ax[isub].legend()
        ax[isub].set_yscale("log")

    ax[-1].set_xlabel(f"{parameter_name}")


    fig_filename = f"solution_sens_{qp_solver_ric_alg}.pdf"
    plt.savefig(fig_filename)
    print(f"stored figure as {fig_filename}")
    plt.show()


if __name__ == "__main__":
    # main(param_M_as_state=False, idxp=1)
    main(param_M_as_state=True, idxp=4, qp_solver_ric_alg=0)
