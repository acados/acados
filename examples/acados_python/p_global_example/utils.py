#
# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# Licensed under the 2-Clause BSD License.

import matplotlib.pyplot as plt
import numpy as np
from acados_template import latexify_plot

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
