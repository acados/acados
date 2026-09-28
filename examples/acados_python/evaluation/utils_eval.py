#
# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# Licensed under the 2-Clause BSD License.

from typing import Dict
import matplotlib.pyplot as plt
import numpy as np
from acados_template import latexify_plot


def plot_pendulum_eval(t, U, X_true, eval: Dict, latexify=False, time_label='$t$', x_labels=None, u_labels=None):
    """
    Plots pendulum with evaluation of cost function.
    Params:
        t: time values of the discretization
        U: arrray with shape (N_sim-1, nu) or (N_sim, nu)
        X_true: arrray with shape (N_sim, nx)
        eval: evaluation dictionary
        latexify: latex style plots
    """

    if latexify:
        latexify_plot()

    nx = X_true.shape[1]
    fig, axes = plt.subplots(nx + 2, 1, sharex=True)

    for i in range(nx):
        axes[i].plot(t, X_true[:, i])
        axes[i].grid()
        if x_labels is not None:
            axes[i].set_ylabel(x_labels[i])
        else:
            axes[i].set_ylabel(f'$x_{i}$')

    axes[-2].plot(t[:-1], eval['cost_without_slacks'], label='cost w/o slacks')
    axes[-2].plot(t[:-1], eval['cost'], label='cost with slacks')
    axes[-2].grid()
    axes[-2].legend()
    axes[-2].set_ylabel('cost')

    axes[-1].step(t, np.append([U[0]], U))

    if u_labels is not None:
        axes[-1].set_ylabel(u_labels[0])
    else:
        axes[-1].set_ylabel('$u$')

    axes[-1].set_xlim(t[0], t[-1])
    axes[-1].set_xlabel(time_label)
    axes[-1].grid()

    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, hspace=0.4)

    fig.align_ylabels()

    return fig, axes
