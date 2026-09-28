#
# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# Licensed under the 2-Clause BSD License.

import os
import matplotlib.pyplot as plt
import numpy as np
import casadi as ca

from acados_template import latexify_plot, casadi_length, AcadosModel

def export_pendulum_ode_model() -> AcadosModel:
    model_name = 'pendulum_ode'

    # constants
    M = 1. # mass of the cart [kg] -> now estimated
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
    u = ca.vertcat(F)

    # xdot
    x1_dot      = ca.SX.sym('x1_dot')
    theta_dot   = ca.SX.sym('theta_dot')
    v1_dot      = ca.SX.sym('v1_dot')
    dtheta_dot  = ca.SX.sym('dtheta_dot')

    xdot = ca.vertcat(x1_dot, theta_dot, v1_dot, dtheta_dot)

    # dynamics
    cos_theta = ca.cos(theta)
    sin_theta = ca.sin(theta)
    denominator = M + m - m*cos_theta*cos_theta
    f_expl = ca.vertcat(v1,
                     dtheta,
                     (-m*l*sin_theta*dtheta*dtheta + m*g*cos_theta*sin_theta+F)/denominator,
                     (-m*l*cos_theta*sin_theta*dtheta*dtheta + F*cos_theta+(M+m)*g*sin_theta)/(l*denominator)
                     )

    f_impl = xdot - f_expl

    model = AcadosModel()

    model.f_impl_expr = f_impl
    model.f_expl_expr = f_expl
    model.x = x
    model.xdot = xdot
    model.u = u
    model.name = model_name

    return model


def plot_open_loop_trajectory_pwpol_u(shooting_nodes, X_traj: np.ndarray, U_fine_traj: list, plt_show=True, u_max=None, title=None,
    x_labels = ['$x$', r'$\theta$', '$v$', r'$\dot{\theta}$'],
    idxpx=None, idxpu=None
                  ):

    nx = X_traj.shape[1]
    nu = U_fine_traj[0].shape[1]
    t = shooting_nodes

    if idxpx is None:
        idxpx = list(range(nx))
    if idxpu is None:
        idxpu = list(range(nu))

    nxpx = len(idxpx)
    nxpu = len(idxpu)

    nrows = max(nxpx, nxpu)

    latexify_plot()
    fig, axes = plt.subplots(ncols=1, nrows=nxpx+nxpu, figsize=(5.5, 1.65*(nxpx+nxpu+1)), sharex=True)
    plt.subplot(nx+1, 1, 1)

    # plot U
    for i, u_interval in enumerate(U_fine_traj):
        n_pol_interval = u_interval.shape[0]
        t_grid_interval = np.linspace(shooting_nodes[i], shooting_nodes[i+1], n_pol_interval)
        for isub, iu in enumerate(idxpu):
            axes[isub].plot(t_grid_interval, u_interval[:, iu], color='C0', alpha= .5 + .5*(i%2) )
    for isub, iu in enumerate(idxpu):
        axes[isub].set_ylabel(f'$u_{iu}$')
        axes[isub].grid()

    for isub, ix in enumerate(idxpx):
        axes[isub+nxpu].plot(t, X_traj[:, ix])

        axes[isub+nxpu].set_ylabel(x_labels[ix])
        axes[isub+nxpu].grid()

    axes[-1].set_xlabel('$t$')
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, hspace=0.4)

    if title is not None:
        axes[0].set_title(title)

    # avoid plotting when running on Travis
    if os.environ.get('ACADOS_ON_CI') is None and plt_show:
        plt.show()
