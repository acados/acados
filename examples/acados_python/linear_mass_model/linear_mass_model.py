#
# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# Licensed under the 2-Clause BSD License.

from acados_template import AcadosModel, latexify_plot
import numpy as np
from casadi import SX, vertcat
import matplotlib.pyplot as plt

def export_linear_mass_model():
    model_name = 'linear_mass'

    # set up states & controls
    qx = SX.sym('qx')
    qy = SX.sym('qy')
    vx = SX.sym('vx')
    vy = SX.sym('vy')
    x = vertcat(qx, qy, vx, vy)

    ux = SX.sym('ux')
    uy = SX.sym('uy')
    u = vertcat(ux, uy)

    f_expl = vertcat(vx, vy, u)
    model = AcadosModel()

    model.f_expl_expr = f_expl
    model.x = x
    model.u = u
    model.p = []
    model.name = model_name

    return model


def plot_linear_mass_system_X_state_space(simX, latexify=False, circle=None, x_goal=None):
    """
    Params:
        simX: x trajectory
        latexify: latex style plots
    """

    if latexify:
        latexify_plot()

    fig, axs = plt.subplots(1, 1)
    if x_goal is not None:
        plt.plot(x_goal[0], x_goal[1], 'rx')
    if circle is not None:
        obs_x, obs_y, obs_rad = circle
        ts = np.linspace(0,2*np.pi,100)
        plt.plot(obs_rad * np.cos(ts)+obs_x,obs_rad * np.sin(ts)-obs_y, 'r')

    plt.grid()
    plt.plot(simX[:,0], simX[:,1], '*-b')
    plt.title('state space plot')

    axs.axis('equal')
    plt.show()


def plot_linear_mass_system_U(shooting_nodes, simU, latexify=False):
    """
    Params:
        simU: u trajectory
        latexify: latex style plots
    """

    if latexify:
        latexify_plot()

    nu = simU.shape[1]
    for i in range(nu):
        plt.subplot(nu, 1, i+1)
        line, = plt.step(shooting_nodes, np.append([simU[0,i]], simU[:,i]))
        plt.grid()
    plt.show()


def plot_linear_mass_system_X(shooting_nodes, simX, latexify=False):
    """
    Params:
        simX: x trajectory
        latexify: latex style plots
    """

    if latexify:
        latexify_plot()

    nx = simX.shape[1]
    for i in range(nx):
        plt.subplot(nx, 1, i+1)
        line, = plt.plot(shooting_nodes, simX[:,i])
        plt.grid()
    plt.show()


