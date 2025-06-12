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

from typing import Optional

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
