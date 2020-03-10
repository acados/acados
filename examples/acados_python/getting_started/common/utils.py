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

import os
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

def plot_pendulum(h, u_max, U, X_true, X_est=None, Y_measured=None, latexify=False):
    """
    Params:
        h: time step
        u_max: maximum absolute value of u
        U: arrray with shape (N_sim-1, nu) or (N_sim, nu)
        X_true: arrray with shape (N_sim, nx)
        X_est: arrray with shape (N_sim-N_mhe, nx)
        Y_measured: array with shape (N_sim, ny)
        latexify: latex style plots
    """

    # latexify plot
    if latexify:
        params = {'backend': 'ps',
                'text.latex.preamble': [r"\usepackage{gensymb} \usepackage{amsmath}"],
                'axes.labelsize': 10,
                'axes.titlesize': 10,
                'legend.fontsize': 10,
                'xtick.labelsize': 10,
                'ytick.labelsize': 10,
                'text.usetex': True,
                'font.family': 'serif'
        }

        matplotlib.rcParams.update(params)

    WITH_ESTIMATION = X_est is not None and Y_measured is not None

    N_sim = X_true.shape[0]
    nx = X_true.shape[1]

    Tf = N_sim*h
    t = np.linspace(0.0, Tf, N_sim)

    if WITH_ESTIMATION:
        N_mhe = N_sim - X_est.shape[0]
        t_mhe = np.linspace(N_mhe, Tf, N_sim)

    plt.subplot(nx+1, 1, 1)
    plt.step(t[:U.shape[0]], U, color='r')
    plt.title('closed-loop simulation')
    plt.ylabel('$u$')
    plt.xlabel('$t$')
    plt.hlines(u_max, t[0], t[-2], linestyles='dashed', alpha=0.7)
    plt.hlines(-u_max, t[0], t[-2], linestyles='dashed', alpha=0.7 )
    plt.ylim([1.2*u_max, -1.2*u_max])
    plt.grid()

    states_lables = ['$x$', r'$\theta$', '$v$', r'$\dot{\theta}$']

    for i in range(nx):
        plt.subplot(nx+1, 1, i+2)
        plt.plot(t, X_true[:,i], label='true')

        if WITH_ESTIMATION:
            plt.plot(t_mhe, X_est[:, i], label='estimated')
            plt.plot(t, Y_measured[:, i], 'x', label='measured')

        plt.ylabel(states_lables[i])
        plt.xlabel('$t$')
        plt.grid()
        plt.legend(loc=1)

    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, hspace=0.4)

    # avoid plotting when running on Travis
    if os.environ.get('ACADOS_ON_TRAVIS') is None:
        plt.show()
