#
# Copyright 2019 Gianluca Frison, Dimitris Kouzoupis, Robin Verschueren,
# Andrea Zanelli, Niels van Duijkeren, Jonathan Frey, Tommaso Sartor,
# Branimir Novoselnik, Rien Quirynen, Rezart Qelibari, Dang Doan,
# Jonas Koenemann, Yutao Chen, Tobias Schöls, Jonas Schlagenhauf, Moritz Diehl
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

from acados_template import AcadosModel
import numpy as np
from casadi import SX, vertcat
import matplotlib
import matplotlib.pyplot as plt

def export_linear_mass_model():
    model_name = 'linear_mass'

    # set up states & controls
    qx = SX.sym('qx')
    qy = SX.sym('qy')
    vx = SX.sym('vx')
    vy = SX.sym('vy')
    x = vertcat(qx, qy, vx, vy)

    # u = SX.sym('u', 2)
    ux = SX.sym('ux')
    uy = SX.sym('uy')
    u = vertcat(ux, uy)

    f_expl = vertcat(vx, vy, u)
    model = AcadosModel()

    # model.f_impl_expr = f_impl
    model.f_expl_expr = f_expl
    # model.disc_dyn_expr
    model.x = x
    model.u = u
    # model.xdot = xdot
    # model.z = z
    model.p = []
    model.name = model_name

    return model


def plot_linear_mass_system_X_state_space(simX, latexify=False, circle=None, x_goal=None):
    """
    Params:
        simX: x trajectory
        latexify: latex style plots
    """

    # latexify plot
    if latexify:
        params = {'backend': 'ps',
                'text.latex.preamble': r"\usepackage{gensymb} \usepackage{amsmath}",
                'axes.labelsize': 10,
                'axes.titlesize': 10,
                'legend.fontsize': 10,
                'xtick.labelsize': 10,
                'ytick.labelsize': 10,
                'text.usetex': True,
                'font.family': 'serif'
        }

        matplotlib.rcParams.update(params)

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
    return

def plot_linear_mass_system_U(shooting_nodes, simU, latexify=False):
    """
    Params:
        simU: u trajectory
        latexify: latex style plots
    """
    nu = simU.shape[1]
    for i in range(nu):
        plt.subplot(nu, 1, i+1)
        line, = plt.step(shooting_nodes, np.append([simU[0,i]], simU[:,i]))
        plt.grid()
    plt.show()

    return

def plot_linear_mass_system_X(shooting_nodes, simX, latexify=False):
    """
    Params:
        simX: x trajectory
        latexify: latex style plots
    """
    nx = simX.shape[1]
    for i in range(nx):
        plt.subplot(nx, 1, i+1)
        line, = plt.plot(shooting_nodes, simX[:,i])
        plt.grid()
    plt.show()

    return


