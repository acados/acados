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


import shutil
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from typing import Optional, List

def latexify_plot() -> None:
    text_usetex = True if shutil.which('latex') else False
    params = {
            'text.latex.preamble': r"\usepackage{gensymb} \usepackage{amsmath}",
            'axes.labelsize': 12,
            'axes.titlesize': 12,
            'legend.fontsize': 12,
            'xtick.labelsize': 12,
            'ytick.labelsize': 12,
            'text.usetex': text_usetex,
            'font.family': 'serif'
    }

    matplotlib.rcParams.update(params)
    return


def plot_convergence(residuals: list,
                     list_labels: list,
                     xlim: Optional[float] = None,
                     ylim: tuple = None,
                     fig_filename: str = None):
    latexify_plot()

    assert len(residuals) == len(list_labels), f"Lists of data and labels do not have the same length, got {len(residuals)} and {len(list_labels)}"

    plt.figure(figsize=(4.5, 3.0))
    for i in range(len(residuals)):
        iters = np.arange(0, len(residuals[i]))
        data = np.array(residuals[i]).squeeze()
        plt.semilogy(iters, data, label=list_labels[i])
    plt.legend(loc='best')
    plt.xlabel("iteration number")
    plt.ylabel("KKT residual norm")
    if ylim is not None:
        plt.ylim(ylim)
    if xlim is not None:
        plt.xlim(0, xlim)
    else:
        plt.xlim(0, max([len(data) for data in residuals]))
    plt.tight_layout()
    plt.grid()
    if fig_filename is not None:
        plt.savefig(fig_filename, dpi=300, bbox_inches='tight', pad_inches=0.01)
    plt.show()

def plot_contraction_rates(rates_list: list,
                          labels: list,
                          fig_filename: str = None):
    latexify_plot()
    plt.figure(figsize=(4.5, 3.0))
    for rates, label in zip(rates_list, labels):
        iters = np.arange(0, len(rates))
        plt.plot(iters, rates, label=label)
    plt.legend(loc='best')
    plt.xlabel("iteration number")
    plt.ylabel("empirical contraction rate")
    plt.xlim(0, max([len(data) for data in rates_list]))
    plt.ylim(0, 1.1)
    plt.tight_layout()
    plt.grid()
    if fig_filename is not None:
        plt.savefig(fig_filename, dpi=300, bbox_inches='tight', pad_inches=0.01)
    plt.show()


def plot_trajectories(
    x_traj_list: List[np.array],
    u_traj_list: List[np.array],
    labels_list: List[str],
    time_traj_list: List[np.array],
    x_labels=None,
    u_labels=None,
    idxbu=[],
    lbu=None,
    ubu=None,
    X_ref=None,
    U_ref=None,
    fig_filename=None,
    x_min=None,
    x_max=None,
    title=None,
    idxpx=None,
    idxpu=None,
    color_list=None,
    linestyle_list=None,
    single_column = False,
    alpha_list = None,
    time_label = None,
    idx_xlogy = None,
    show_legend = True,
    bbox_to_anchor = None,
    ncol_legend = 2,
    figsize=None,
    show_plot: bool = True,
    latexify: bool = True,
):
    if latexify:
        latexify_plot()

    nx = x_traj_list[0].shape[1]
    nu = u_traj_list[0].shape[1]
    Ntraj = len(x_traj_list)

    if idxpx is None:
        idxpx = list(range(nx))
    if idxpu is None:
        idxpu = list(range(nu))

    if color_list is None:
        color_list = [f"C{i}" for i in range(Ntraj)]
    if linestyle_list is None:
        linestyle_list = Ntraj * ['-']
    if alpha_list is None:
        alpha_list = Ntraj * [0.8]

    if idx_xlogy is None:
        idx_xlogy = []

    if time_label is None:
        time_label = "$t$"

    if x_labels is None:
        x_labels = [f"$x_{i}$" for i in range(nx)]
    if u_labels is None:
        u_labels = [f"$u_{i}$" for i in range(nu)]

    nxpx = len(idxpx)
    nxpu = len(idxpu)
    nrows = max(nxpx, nxpu)

    if figsize is None:
        if single_column:
            figsize = (6.0, 2*(nxpx+nxpu+1))
        else:
            figsize = (10, (nxpx+nxpu))

    if single_column:
        fig, axes = plt.subplots(ncols=1, nrows=nxpx+nxpu, figsize=figsize, sharex=True)
    else:
        fig, axes = plt.subplots(ncols=2, nrows=nrows, figsize=figsize, sharex=True)
        axes = np.ravel(axes, order='F')

    if title is not None:
        axes[0].set_title(title)

    for i in idxpx:
        isubplot = idxpx.index(i)
        for x_traj, time_traj, label, color, linestyle, alpha in zip(x_traj_list, time_traj_list, labels_list, color_list, linestyle_list, alpha_list):
            axes[isubplot].plot(time_traj, x_traj[:, i], label=label, alpha=alpha, color=color, linestyle=linestyle)

        if X_ref is not None:
            axes[isubplot].step(
                time_traj_list[0],
                X_ref[:, i],
                alpha=0.8,
                where="post",
                label="reference",
                linestyle="dotted",
                color="k",
            )
        axes[isubplot].set_ylabel(x_labels[i])
        axes[isubplot].grid()
        axes[isubplot].set_xlim(time_traj_list[0][0], time_traj_list[0][-1])

        if i in idx_xlogy:
            axes[isubplot].set_yscale('log')

        if x_min is not None:
            axes[isubplot].set_ylim(bottom=x_min[i])

        if x_max is not None:
            axes[isubplot].set_ylim(top=x_max[i])

    for i in idxpu:
        for u_traj, time_traj, label, color, linestyle, alpha in zip(u_traj_list, time_traj_list, labels_list, color_list, linestyle_list, alpha_list):
            vals = u_traj[:, i]
            axes[i+nrows].step(time_traj, np.append([vals[0]], vals), label=label, alpha=alpha, color=color, linestyle=linestyle)

        if U_ref is not None:
            axes[i+nrows].step(time_traj, np.append([U_ref[0, i]], U_ref[:, i]), alpha=0.8,
                               label="reference", linestyle="dotted", color="k")

        axes[i+nrows].set_ylabel(u_labels[i])
        axes[i+nrows].grid()

        if i in idxbu:
            axes[i+nrows].hlines(
                ubu[i], time_traj[0], time_traj[-1], linestyles="dashed", alpha=0.4, color="k"
            )
            axes[i+nrows].hlines(
                lbu[i], time_traj[0], time_traj[-1], linestyles="dashed", alpha=0.4, color="k"
            )
            axes[i+nrows].set_xlim(time_traj[0], time_traj[-1])
            bound_margin = 0.05
            u_lower = (1-bound_margin) * lbu[i] if lbu[i] > 0 else (1+bound_margin) * lbu[i]
            axes[i+nrows].set_ylim(bottom=u_lower, top=(1+bound_margin) * ubu[i])

    axes[nxpx+nxpu-1].set_xlabel(time_label)
    if not single_column:
        axes[nxpx-1].set_xlabel(time_label)

    if bbox_to_anchor is None and single_column:
        bbox_to_anchor=(0.5, -0.75)
    elif bbox_to_anchor is None:
        bbox_to_anchor=(0.5, -1.5)

    if show_legend:
        axes[nxpx+nxpu-1].legend(loc="lower center", ncol=ncol_legend, bbox_to_anchor=bbox_to_anchor)

    fig.align_ylabels()
    # fig.tight_layout()

    if not single_column:
        for i in range(nxpu, nxpx):
            fig.delaxes(axes[i+nrows])

    if fig_filename is not None:
        plt.savefig(fig_filename, bbox_inches="tight", transparent=True, pad_inches=0.05)
        print(f"\nstored figure in {fig_filename}")

    if show_plot:
        plt.show()
