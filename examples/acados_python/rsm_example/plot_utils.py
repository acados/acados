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

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from acados_template import latexify_plot

def plot_rsm_trajectories(simX, simU, ocp, Ts):
    latexify_plot()
    Nsim = simU.shape[0]
    t = np.linspace(0.0, Ts*Nsim, Nsim)
    plt.subplot(4, 1, 1)
    plt.step(t, simU[:,0], color='r')
    plt.plot([0, Ts*Nsim], [ocp.cost.yref[2], ocp.cost.yref[2]], '--')
    plt.title('closed-loop simulation')
    plt.ylabel('u_d')
    plt.xlabel('t')
    plt.grid(True)
    plt.subplot(4, 1, 2)
    plt.step(t, simU[:,1], color='r')
    plt.plot([0, Ts*Nsim], [ocp.cost.yref[3], ocp.cost.yref[3]], '--')
    plt.ylabel('u_q')
    plt.xlabel('t')
    plt.grid(True)
    plt.subplot(4, 1, 3)
    plt.plot(t, simX[:,0])
    plt.plot([0, Ts*Nsim], [ocp.cost.yref[0], ocp.cost.yref[0]], '--')
    plt.ylabel('psi_d')
    plt.xlabel('t')
    plt.grid(True)
    plt.subplot(4, 1, 4)
    plt.plot(t, simX[:,1])
    plt.plot([0, Ts*Nsim], [ocp.cost.yref[1], ocp.cost.yref[1]], '--')
    plt.ylabel('psi_q')
    plt.xlabel('t')
    plt.grid(True)

def plot_hexagon(simU, u_max):
    latexify_plot()

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
    # -q2 <= uq  <= q2

    plt.figure()
    plt.plot(simU[:,0], simU[:,1], 'o')
    plt.xlabel('ud')
    plt.ylabel('uq')
    ud = np.linspace(-1.5*u_max, 1.5*u_max, 100)
    plt.plot(ud, -m1*ud -q1)
    plt.plot(ud, -m1*ud +q1)
    plt.plot(ud, +m1*ud -q1)
    plt.plot(ud, +m1*ud +q1)
    plt.plot(ud, -q2*np.ones((100, 1)))
    plt.plot(ud, q2*np.ones((100, 1)))
    plt.grid(True)
    ax = plt.gca()
    ax.set_xlim([-1.5*u_max, 1.5*u_max])
    ax.set_ylim([-1.5*u_max, 1.5*u_max])
    circle = plt.Circle((0, 0), u_max*np.sqrt(3)/2, color='red', fill=False)
    ax.add_artist(circle)