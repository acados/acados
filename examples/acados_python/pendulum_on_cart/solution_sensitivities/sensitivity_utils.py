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
import matplotlib.pyplot as plt
from acados_template import AcadosModel, AcadosOcp, AcadosOcpSolver, latexify_plot

latexify_plot()

def export_pendulum_ode_model_with_mass_as_p_global(dt) -> AcadosModel:

    model_name = 'pendulum_parametric'

    # constants
    m = 0.1 # mass of the ball [kg]
    g = 9.81 # gravity constant [m/s^2]
    l = 0.8 # length of the rod [m]

    # set up states & controls
    x1      = ca.SX.sym('x1')
    theta   = ca.SX.sym('theta')
    v1      = ca.SX.sym('v1')
    dtheta  = ca.SX.sym('dtheta')

    x = ca.vertcat(x1, theta, v1, dtheta)

    F = ca.SX.sym('F')
    u = F

    # xdot
    x1_dot      = ca.SX.sym('x1_dot')
    theta_dot   = ca.SX.sym('theta_dot')
    v1_dot      = ca.SX.sym('v1_dot')
    dtheta_dot  = ca.SX.sym('dtheta_dot')

    xdot = ca.vertcat(x1_dot, theta_dot, v1_dot, dtheta_dot)

    # parameters
    m_cart = ca.SX.sym('m_cart')  # mass of the cart [kg]
    p = m_cart

    # dynamics
    cos_theta = ca.cos(theta)
    sin_theta = ca.sin(theta)
    denominator = m_cart + m - m*cos_theta*cos_theta
    f_expl = ca.vertcat(v1,
                       dtheta,
                       (-m*l*sin_theta*dtheta*dtheta + m*g*cos_theta*sin_theta+F)/denominator,
                       (-m*l*cos_theta*sin_theta*dtheta*dtheta + F*cos_theta+(m_cart+m)*g*sin_theta)/(l*denominator)
                       )

    f_impl = xdot - f_expl

    # discrete dynamics via RK4
    ode = ca.Function('ode', [x, u, p], [f_expl])
    k1 = ode(x, u, p)
    k2 = ode(x + dt/2*k1, u, p)
    k3 = ode(x + dt/2*k2, u, p)
    k4 = ode(x + dt*k3, u, p)

    xf = x + dt/6 * (k1 + 2*k2 + 2*k3 + k4)

    model = AcadosModel()

    model.disc_dyn_expr = xf
    model.f_impl_expr = f_impl
    model.f_expl_expr = f_expl
    model.x = x
    model.xdot = xdot
    model.u = u
    model.p_global = p
    model.name = model_name

    # store meta information
    model.x_labels = ['$x$ [m]', r'$\theta$ [rad]', '$v$ [m]', r'$\dot{\theta}$ [rad/s]']
    model.u_labels = ['$F$']
    model.t_label = '$t$ [s]'

    return model


def export_parametric_ocp(
    x0=np.array([0.0, np.pi / 6, 0.0, 0.0]), N_horizon=50, T_horizon=2.0, Fmax=80.0,
    hessian_approx = "GAUSS_NEWTON", qp_solver_ric_alg=1,
    cost_scale_as_param=False,
    with_parametric_constraint=True,
    with_nonlinear_constraint=True
) -> AcadosOcp:

    ocp = AcadosOcp()
    dt = T_horizon/N_horizon
    ocp.model = export_pendulum_ode_model_with_mass_as_p_global(dt=dt)

    ocp.solver_options.N_horizon = N_horizon

    # set mass to one
    ocp.p_global_values = np.ones((1,))

    Q_mat = 2 * np.diag([1e3, 1e3, 1e-2, 1e-2])
    R_mat = 2 * np.diag([1e-1])

    ocp.cost.cost_type = "EXTERNAL"
    ocp.cost.cost_type_e = "EXTERNAL"

    if cost_scale_as_param:
        # add parameter to model
        cost_scale_param = ca.SX.sym('cost_scale_param')
        ocp.model.p_global = ca.vertcat(ocp.model.p_global, cost_scale_param)
        ocp.p_global_values = np.concatenate((ocp.p_global_values, np.ones((1,))))
        # add nonlinear dependency in cost
        cost_scale_factor = ca.exp(cost_scale_param)
    else:
        cost_scale_factor = 1.0

    # NOTE here we make the cost parametric
    ocp.model.cost_expr_ext_cost = cost_scale_factor * ocp.model.x.T @ Q_mat @ ocp.model.x + ocp.model.u.T @ R_mat @ ocp.model.u
    ocp.model.cost_expr_ext_cost_e = cost_scale_factor * ocp.model.x.T @ Q_mat @ ocp.model.x

    ocp.constraints.lbu = np.array([-Fmax])
    ocp.constraints.ubu = np.array([+Fmax])
    ocp.constraints.idxbu = np.array([0])

    if with_parametric_constraint:
        if with_nonlinear_constraint:
            ocp.model.con_h_expr = -ocp.model.x[0] * ocp.model.p_global[0]**2
        else:
            ocp.model.con_h_expr = -ocp.model.x[0] * ocp.model.p_global[0]
        ocp.constraints.lh = np.array([-1.5])
        ocp.constraints.uh = np.array([1.5])

    ocp.constraints.x0 = x0

    ocp.solver_options.qp_solver = "PARTIAL_CONDENSING_HPIPM"
    ocp.solver_options.integrator_type = "DISCRETE"
    ocp.solver_options.nlp_solver_type = "SQP"
    ocp.solver_options.qp_solver_cond_N = N_horizon

    ocp.solver_options.tf = T_horizon

    ocp.solver_options.qp_solver_ric_alg = qp_solver_ric_alg
    ocp.solver_options.qp_solver_cond_ric_alg = qp_solver_ric_alg
    ocp.solver_options.qp_solver_mu0 = 1e3  # makes HPIPM converge more robustly
    ocp.solver_options.hessian_approx = hessian_approx
    ocp.solver_options.nlp_solver_max_iter = 400
    ocp.solver_options.tol = 1e-8
    # ocp.solver_options.globalization = "MERIT_BACKTRACKING"

    # if hessian_approx == 'EXACT':
        # sensitivity solver settings!
    ocp.solver_options.with_solution_sens_wrt_params = True
    ocp.solver_options.with_value_sens_wrt_params = True

    return ocp

def plot_cost_gradient_results(p_test, cost_values, acados_cost_grad, np_cost_grad,
                               cost_reconstructed_np_grad, cost_reconstructed_acados=None,
                               y_scale_log=True, xlabel=None, title=None):
    _, ax = plt.subplots(nrows=4, ncols=1, sharex=True, figsize=(9,9))

    ax[0].plot(p_test, cost_values, label='cost acados', color='k')
    ax[0].plot(p_test, cost_reconstructed_np_grad, "--", label='reconstructed from finite diff')
    if cost_reconstructed_acados is not None:
        ax[0].plot(p_test, cost_reconstructed_acados, ":", label='reconstructed from acados derivatives')
    ax[0].set_ylabel(r"cost")

    ax[1].plot(p_test, np.abs(np_cost_grad), "--", label='finite diff')
    ax[1].plot(p_test, np.abs(acados_cost_grad), ":", label='acados')
    ax[1].set_ylabel(r"$|\partial_p V^*|$")

    if y_scale_log:
        ax[1].set_yscale("log")

    # plot differences
    isub = 2
    ax[isub].plot(p_test, np.abs(acados_cost_grad - np_cost_grad), "--", label='acados vs. finite diff')
    ax[isub].set_ylabel(r"abs diff $\partial_p V^*$")

    if y_scale_log:
        ax[isub].set_yscale("log")

    isub += 1
    ax[isub].plot(p_test, np.abs(acados_cost_grad - np_cost_grad) / np.abs(np_cost_grad), "--", label='acados vs. finite diff')
    ax[isub].set_ylabel(r"rel. diff. $\partial_p V^*$")

    if y_scale_log:
        ax[isub].set_yscale("log")

    for i in range(isub+1):
        ax[i].grid()
        ax[i].legend()

    if xlabel is not None:
        ax[-1].set_xlabel(xlabel)

    if title is not None:
        ax[0].set_title(title)
    ax[-1].set_xlim([p_test[0], p_test[-1]])

    fig_filename = f"cost_gradient.pdf"
    plt.savefig(fig_filename)
    print(f"stored figure as {fig_filename}")
    plt.show()


def plot_solution_sensitivities_results(p_test, pi, pi_reconstructed_acados, pi_reconstructed_np_grad, sens_u, np_grad,
                 min_eig_full=None, min_eig_proj_hess=None, min_eig_P=None,
                 min_abs_eig_full=None, min_abs_eig_proj_hess=None, min_abs_eig_P=None,
                 eigen_analysis=False, title=None, parameter_name="",
                 max_lam_parametric_constraint=None, sum_lam_parametric_constraint=None,
                 multipliers_bu=None, multipliers_h=None, plot_reconstructed=True,
                 figsize=None,
                 ):

    nsub = 5 if eigen_analysis else 3

    with_multiplier_subplot = max_lam_parametric_constraint is not None or sum_lam_parametric_constraint is not None or multipliers_bu is not None or multipliers_h is not None
    if with_multiplier_subplot:
        nsub += 1

    if figsize is None:
        figsize = (9, 9)
    _, ax = plt.subplots(nrows=nsub, ncols=1, sharex=True, figsize=figsize)

    isub = 0
    ax[isub].plot(p_test, pi, label='solution acados', color='k')
    if plot_reconstructed:
        ax[isub].plot(p_test, pi_reconstructed_acados, "--", label='reconstructed from acados solution sens.')
        ax[isub].plot(p_test, pi_reconstructed_np_grad, "--", label='reconstructed from finite diff.')
    ax[isub].set_ylabel(r"$u_0$")
    if title is not None:
        ax[isub].set_title(title)

    isub += 1
    ax[isub].plot(p_test, sens_u, label="acados")
    ax[isub].plot(p_test, np_grad, "--", label="finite diff.")
    ax[isub].set_xlim([p_test[0], p_test[-1]])
    ax[isub].set_ylabel(r"$\partial_\theta u_0$")

    isub += 1
    ax[isub].plot(p_test, np.abs(sens_u- np_grad), "--", label='acados - finite diff.')
    ax[isub].set_ylabel(r"difference $\partial_\theta u_0$")
    ax[isub].set_yscale("log")

    if with_multiplier_subplot:
        isub += 1
        isub_multipliers = isub
        if max_lam_parametric_constraint is not None:
            ax[isub].plot(p_test, max_lam_parametric_constraint, label=r'max $\lambda$ parametric constraint')
        # ax[isub].set_ylabel("max lam parametric constraint")
        # ax[isub].set_yscale("log")
        if sum_lam_parametric_constraint is not None:
            ax[isub].plot(p_test, sum_lam_parametric_constraint, label=r'sum $\lambda$ parametric constraint')

        legend_elements = []
        if multipliers_bu is not None:
            for lam in multipliers_bu:
                ax[isub].plot(p_test, lam, linestyle='--', color='C0', alpha=.6)
            legend_elements += [plt.Line2D([0], [0], color='C0', linestyle='--', label='multipliers control bounds')]
        if multipliers_h is not None:
            for lam in multipliers_h:
                ax[isub].plot(p_test, lam, linestyle='--', color='C1', alpha=.6)
            legend_elements += [plt.Line2D([0], [0], color='C1', linestyle='--', label='multipliers $h$')]
        ax[isub].legend(handles=legend_elements, ncol=2)
        ax[isub].set_ylim([0, 14])
        ax[isub].set_ylabel("multipliers")
    if eigen_analysis:
        isub += 1
        ax[isub].plot(p_test, np.sign(min_eig_full), linestyle="-", alpha=.6, label='full Hessian')
        ax[isub].plot(p_test, np.sign(min_eig_proj_hess), "--", label='proj Hessian')
        ax[isub].plot(p_test, np.sign(min_eig_P), ":", label='$P$ Riccati')
        ax[isub].set_ylabel("sign min eig")

        isub += 1
        ax[isub].plot(p_test, min_abs_eig_full, "--", label='full Hessian')
        ax[isub].plot(p_test, min_abs_eig_proj_hess, "--", label='proj Hessian')
        ax[isub].plot(p_test, min_abs_eig_P, "--", label='$P$ Riccati')
        ax[isub].set_ylabel("min abs eig")
        ax[isub].set_yscale("log")

    for isub in range(nsub):
        ax[isub].grid()
        if isub != isub_multipliers:
            ax[isub].legend()

    ax[-1].set_xlabel(f"{parameter_name}")

    plt.tight_layout()

    fig_filename = f"solution_sens_{title}.pdf"
    plt.savefig(fig_filename)
    print(f"stored figure as {fig_filename}")
    plt.show()


def plot_smoothed_solution_sensitivities_results(p_test, pi_label_pairs, sens_pi_label_pairs,
                 title=None, parameter_name="",
                 multipliers_bu=None, multipliers_h=None,
                 figsize=None,
                 fig_filename=None,
                 horizontal_plot=False,
                 ):

    nsub = 2

    with_multiplier_subplot = multipliers_bu is not None or multipliers_h is not None
    if with_multiplier_subplot:
        nsub += 1

    if figsize is None:
        figsize = (9, 9)
    if not horizontal_plot:
        _, ax = plt.subplots(nrows=nsub, ncols=1, sharex=True, figsize=figsize)
    else:
        _, ax = plt.subplots(nrows=1, ncols=nsub, sharex=False, figsize=figsize)

    linestyles = ["-", "--", "-.", ":", "-", "--", "-.", ":"]

    isub = 0
    for i, (pi, label) in enumerate(pi_label_pairs):
        ax[isub].plot(p_test, pi, label=label, linestyle=linestyles[i])
    ax[isub].set_ylabel(r"$u_0$")
    if title is not None:
        ax[isub].set_title(title)
    ax[isub].legend()
    ax[isub].legend(handlelength=1.2)

    isub += 1
    for i, (sens_pi, label) in enumerate(sens_pi_label_pairs):
        ax[isub].plot(p_test, sens_pi, label=label, linestyle=linestyles[i])
    ax[isub].set_ylabel(r"$\partial_\theta u_0$")
    if horizontal_plot:
        ax[isub].legend(loc = 'upper left', handlelength=1.2, ncol=2, columnspacing=0.5, labelspacing=0.2)
    else:
        ax[isub].legend(loc = 'upper left', handlelength=1.2)

    if with_multiplier_subplot:
        isub += 1
        legend_elements = []
        if multipliers_bu is not None:
            for lam in multipliers_bu:
                ax[isub].plot(p_test, lam, linestyle='--', color='C0', alpha=.6)
            legend_elements += [plt.Line2D([0], [0], color='C0', linestyle='--', label='multipliers control bounds')]
        if multipliers_h is not None:
            for lam in multipliers_h:
                ax[isub].plot(p_test, lam, linestyle='--', color='C1', alpha=.6)
            legend_elements += [plt.Line2D([0], [0], color='C1', linestyle='--', label='multipliers $h$')]
        if horizontal_plot:
            ax[isub].legend(handles=legend_elements, ncol=1, handlelength=1.2)
        else:
            ax[isub].legend(handles=legend_elements, ncol=2)
        ax[isub].set_ylim([0, 14])
        ax[isub].set_ylabel("multipliers")

    for isub in range(nsub):
        ax[isub].grid()
        ax[isub].set_xlim([p_test[0], p_test[-1]])

        if horizontal_plot:
            ax[isub].set_xlabel(f"{parameter_name}")

    ax[-1].set_xlabel(f"{parameter_name}")

    plt.tight_layout()
    if fig_filename is not None:
        plt.savefig(fig_filename)
        print(f"stored figure as {fig_filename}")
    plt.show()




def plot_pendulum(t, u_max, U, X_true, latexify=False, plt_show=True, time_label='$t$', x_labels=None, u_labels=None):
    """
    Params:
        t: time values of the discretization
        u_max: maximum absolute value of u
        U: arrray with shape (N_sim-1, nu) or (N_sim, nu)
        X_true: arrray with shape (N_sim, nx)
        latexify: latex style plots
    """

    if latexify:
        latexify_plot()

    nx = X_true.shape[1]
    fig, axes = plt.subplots(nx+1, 1, sharex=True)

    for i in range(nx):
        axes[i].plot(t, X_true[:, i])
        axes[i].grid()
        if x_labels is not None:
            axes[i].set_ylabel(x_labels[i])
        else:
            axes[i].set_ylabel(f'$x_{i}$')

    axes[-1].step(t, np.append([U[0]], U))

    if u_labels is not None:
        axes[-1].set_ylabel(u_labels[0])
    else:
        axes[-1].set_ylabel('$u$')

    axes[-1].hlines(u_max, t[0], t[-1], linestyles='dashed', alpha=0.7)
    axes[-1].hlines(-u_max, t[0], t[-1], linestyles='dashed', alpha=0.7)
    axes[-1].set_ylim([-1.2*u_max, 1.2*u_max])
    axes[-1].set_xlim(t[0], t[-1])
    axes[-1].set_xlabel(time_label)
    axes[-1].grid()

    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, hspace=0.4)

    fig.align_ylabels()

    if plt_show:
        plt.show()
