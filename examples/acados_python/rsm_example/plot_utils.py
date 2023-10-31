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

import matplotlib.pyplot as plt
import numpy as np
from acados_template import latexify_plot

latexify_plot()

def plot_rsm_trajectories(simX, simU, simY, Ts):
    Nsim = simU.shape[0]
    t = np.linspace(0.0, Ts*Nsim, Nsim)

    fig, axes = plt.subplots(ncols=1, nrows=4, sharex=True)

    for i in range(2):
        axes[i].step(t, simU[:,i], color='C1')
        axes[i].step(t, simY[:, 2+i], '--', color='k')
        axes[i].grid(True)

    for i in range(2):
        axes[2+i].plot(t, simX[:,i], color='C2')
        axes[2+i].step(t, simY[:, i], '--', color='k')
        axes[2+i].grid(True)

    axes[0].set_title('closed-loop simulation')
    axes[0].set_ylabel(r'$u^q$')
    axes[1].set_ylabel(r'$u^d$')
    axes[2].set_ylabel(r'$\psi^d$')
    axes[3].set_ylabel(r'$\psi^q$')
    axes[-1].set_xlabel(r'$t$')
    axes[0].set_xlim(t[0], t[-1])


def plot_hexagon(simU, u_max, Uref):

    Nsim = simU.shape[0]
    r = u_max

    x1 = r
    y1 = 0
    x2 = r*np.cos(np.pi/3)
    y2 = r*np.sin(np.pi/3)

    q1 = -(y2 - y1/x1*x2)/(1-x2/x1)
    m1 = -(y1 + q1)/x1

    # q1 <= uq + m1*ud <= -q1
    # q1 <= uq - m1*ud <= -q1

    # box constraints
    q2 = r*np.sin(np.pi/3)

    ud = np.linspace(-1.5*u_max, 1.5*u_max, 100)

    ms = [0., m1, -m1, 0., m1, -m1]
    qs = [q2, -q1, q1, -q2,  q1, -q1]

    crossings_ud = [(qs[i]-qs[i-1])/(ms[i-1]-ms[i]) for i in range(len(ms))]
    crossings_uq = [ms[i]*crossings_ud[i] + qs[i] for i in range(len(ms))]

    fig, ax = plt.subplots(ncols=1, nrows=1)
    ax.fill(crossings_ud, crossings_uq, color='k', alpha=0.15)

    for m, q in zip(ms, qs):
        ax.plot(ud, m*ud + q, color='k', alpha=0.3)

    circle = plt.Circle((0, 0), u_max*np.sqrt(3)/2, color='k', fill=False, alpha=0.8)

    cmap = plt.get_cmap('viridis_r')
    ax.scatter(simU[:,0], simU[:,1], color=[cmap(tau) for tau in np.linspace(0.1, 1., Nsim)])
    ax.scatter(Uref[:,0], Uref[:,1], color='k', alpha=0.9, marker='x')
    ax.grid(True)
    ax.set_ylabel('$u^q$')
    ax.set_xlabel('$u^d$')
    ax.set_xlim([-1.5*u_max, 1.5*u_max])
    ax.set_ylim([-1.5*u_max, 1.5*u_max])
    ax.set_aspect('equal', 'box')
    ax.add_artist(circle)
