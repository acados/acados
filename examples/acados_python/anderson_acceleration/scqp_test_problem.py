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


import casadi as ca
import numpy as np
from acados_template import AcadosModel, AcadosOcp, ACADOS_INFTY, latexify_plot
import matplotlib.pyplot as plt


NX = 4
GOAL_POSITION_RADIUS = (5e-02)**2
LENGTH_PENDULUM = 0.8

def create_pendulum_model():
    M = 1.0
    m = 0.1
    g = 9.81
    states = ca.MX.sym("x", NX)
    states_dot = ca.MX.sym('xdot', NX)
    [p, pDot, theta, thetaDot] = ca.vertsplit(states)
    controls = ca.MX.sym("u")
    F = controls
    denominator = M + m - m*ca.cos(theta)*ca.cos(theta)
    f_x_p = pDot
    f_x_pDot = (-m*LENGTH_PENDULUM*ca.sin(theta)*thetaDot*thetaDot + m*g*ca.cos(theta)*ca.sin(theta)+F)/denominator
    f_x_theta = thetaDot
    f_x_thetaDot = (-m*LENGTH_PENDULUM*ca.cos(theta)*ca.sin(theta)*thetaDot*thetaDot + F*ca.cos(theta)+(M+m)*g*ca.sin(theta))/(LENGTH_PENDULUM*denominator)

    f_x = ca.vertcat(f_x_p, f_x_pDot, f_x_theta, f_x_thetaDot)
    f_impl = states_dot - f_x

    model = AcadosModel()
    model.f_impl_expr = f_impl
    model.f_expl_expr = f_x
    model.x = states
    model.xdot = states_dot
    model.u = controls
    model.name = 'pendulum'
    return model

def pendulum_position(p, theta):
    return ca.vertcat(p-LENGTH_PENDULUM*ca.sin(theta), LENGTH_PENDULUM*ca.cos(theta))

def pendulum_final_position_constraint(p, theta):
    pendulum_final_position = pendulum_position(p, theta)
    return ca.sumsqr(pendulum_final_position - ca.vertcat(LENGTH_PENDULUM, LENGTH_PENDULUM)) - GOAL_POSITION_RADIUS

def build_acados_test_problem(mode='GN', with_anderson_acceleration=False, globalization="FIXED_STEP", max_iter=400) -> AcadosOcp:
    print(f"Building acados test problem with mode {mode} and with_anderson_acceleration {with_anderson_acceleration}, {globalization}")

    # create ocp object to formulate the OCP
    ocp = AcadosOcp()

    # set model
    model = create_pendulum_model()

    ocp.model = model

    # Define cost
    R_mat = np.array([[1e-4]])

    ocp.cost.cost_type_0 = 'LINEAR_LS'
    ocp.cost.W_0 = R_mat
    ocp.cost.Vx_0 = np.zeros((1, NX))
    ocp.cost.Vu_0 = np.array([[1]])
    ocp.cost.yref_0 = np.zeros((1, ))

    ocp.cost.cost_type = 'LINEAR_LS'
    ocp.cost.W = R_mat
    ocp.cost.Vx = np.zeros((1, NX))
    ocp.cost.Vu = np.array([[1]])
    ocp.cost.yref  = np.zeros((1, ))

    # Define constraints
    ocp.constraints.x0 = np.array([0.0, 0.0, np.pi, 0.0])

    # Convex over Nonlinear Constraints
    if mode in ['GN', 'EXACT']:
        ocp.model.con_h_expr_e = pendulum_final_position_constraint(model.x[0], model.x[2])
        ocp.constraints.lh_e = np.array([-ACADOS_INFTY])
        ocp.constraints.uh_e = np.array([0.0])

    elif mode == 'SCQP':
        r = ca.MX.sym('r', 2, 1)
        ocp.model.con_phi_expr_e = ca.sumsqr(r)
        ocp.model.con_r_in_phi_e = r
        ocp.model.con_r_expr_e = pendulum_position(model.x[0], model.x[2]) - ca.vertcat(LENGTH_PENDULUM, LENGTH_PENDULUM)
        ocp.constraints.lphi_e = np.array([-ACADOS_INFTY])
        ocp.constraints.uphi_e = np.array([GOAL_POSITION_RADIUS])
    else:
        raise ValueError("Wrong mode name!")

    # set solver options
    N_horizon = 20
    dt = 0.05

    ocp.solver_options.integrator_type = 'ERK'
    ocp.solver_options.nlp_solver_max_iter = 400
    ocp.solver_options.qp_solver_iter_max = 1000
    ocp.solver_options.nlp_solver_type = 'SQP'
    ocp.solver_options.sim_method_num_steps = 20
    ocp.solver_options.sim_method_num_stages = 4
    ocp.solver_options.N_horizon = N_horizon
    ocp.solver_options.tf = dt*N_horizon
    ocp.solver_options.tol = 1e-12
    ocp.solver_options.qp_tol = 1e-1 * ocp.solver_options.tol
    # ocp.solver_options.cost_scaling = np.ones((N_horizon+1, ))
    ocp.solver_options.with_anderson_acceleration = with_anderson_acceleration
    ocp.solver_options.globalization = globalization
    ocp.solver_options.qp_solver_ric_alg = 0
    ocp.solver_options.qp_solver_cond_ric_alg = 0
    ocp.solver_options.nlp_solver_max_iter = max_iter
    # ocp.solver_options.qp_solver = 'FULL_CONDENSING_DAQP'
    ocp.solver_options.hessian_approx = 'EXACT' if mode == "EXACT" else 'GAUSS_NEWTON'
    if mode == "EXACT":
        ocp.solver_options.exact_hess_dyn = 1

    if max_iter < 10:
        ocp.solver_options.print_level = 4

    return ocp



def plot_pendulum(shooting_nodes, U, X_true, X_est=None, Y_measured=None, latexify=True, plt_show=True, X_true_label=None,
    time_label='$t$', x_labels=['$x$', r'$\theta$', '$v$', r'$\dot{\theta}$'],
    title = None
                  ):
    """
    Params:
        shooting_nodes: time values of the discretization
        u_max: maximum absolute value of u
        U: arrray with shape (N_sim-1, nu) or (N_sim, nu)
        X_true: arrray with shape (N_sim, nx)
        X_est: arrray with shape (N_sim-N_mhe, nx)
        Y_measured: array with shape (N_sim, ny)
        latexify: latex style plots
    """

    if latexify:
        latexify_plot()

    WITH_ESTIMATION = X_est is not None and Y_measured is not None

    N_sim = X_true.shape[0]
    nx = X_true.shape[1]

    Tf = shooting_nodes[N_sim-1]
    t = shooting_nodes

    Ts = t[1] - t[0]
    if WITH_ESTIMATION:
        N_mhe = N_sim - X_est.shape[0]
        t_mhe = np.linspace(N_mhe * Ts, Tf, N_sim-N_mhe)

    plt.figure()
    plt.subplot(nx+1, 1, 1)
    line, = plt.step(t, np.append([U[0]], U))
    if X_true_label is not None:
        line.set_label(X_true_label)
    else:
        line.set_color('r')
    if title is not None:
        plt.title(title)
    plt.ylabel('$u$')
    plt.xlabel(time_label)
    plt.grid()


    for i in range(nx):
        plt.subplot(nx+1, 1, i+2)
        line, = plt.plot(t, X_true[:, i], label='true')
        if X_true_label is not None:
            line.set_label(X_true_label)

        if WITH_ESTIMATION:
            plt.plot(t_mhe, X_est[:, i], '--', label='estimated')
            plt.plot(t, Y_measured[:, i], 'x', label='measured')

        plt.ylabel(x_labels[i])
        plt.xlabel('$t$')
        plt.grid()
        plt.legend(loc=1)

    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, hspace=0.4)

    if plt_show:
        plt.show()
    return
