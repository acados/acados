# -*- coding: future_fstrings -*-
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

import sys
sys.path.insert(0, '../common')

from pendulum_model import export_pendulum_ode_model
from export_mhe_ode_model_with_noisy_param import export_mhe_ode_model_with_noisy_param

from export_ocp_solver import export_ocp_solver
from export_mhe_solver_with_param import export_mhe_solver_with_param

import numpy as np
from utils import plot_pendulum


# general
Tf = 1.0
N = 20
h = Tf/N
Fmax = 80

# NOTE: hard coded in export_pendulum_ode_model;
l_true = 0.8

# ocp model and solver
model = export_pendulum_ode_model()

nx = model.x.size()[0]
nu = model.u.size()[0]

Q_ocp = np.diag([1e3, 1e3, 1e-2, 1e-2])
R_ocp = 1e-2 *np.eye(1)

acados_solver_ocp = export_ocp_solver(model, N, h, Q_ocp, R_ocp, Fmax, use_cython=True)

# mhe model and solver
model_mhe = export_mhe_ode_model_with_noisy_param()

nx_augmented = model_mhe.x.size()[0]
nw = model_mhe.u.size()[0]
ny = nx

Q0_mhe = np.diag([0.1, 0.1, 0.1, 0.1, 1])
Q_mhe  = 10.*np.diag([0.2, 0.2, 2, 2, 0.1])
R_mhe  = 2*np.diag([0.1, 0.1, 0.1, 0.1])

acados_solver_mhe = export_mhe_solver_with_param(model_mhe, N, h, Q_mhe, Q0_mhe, R_mhe, use_cython=True)

# simulation
v_stds = [0.2, 0.5, 1, 1]
v_stds = [0, 0, 0, 0]

simX = np.ndarray((N+1, nx))
simU = np.ndarray((N, nu))
simY = np.ndarray((N+1, nx))

simXest = np.zeros((N+1, nx))
simWest = np.zeros((N, nx_augmented))
sim_l_est = np.zeros((N+1, 1))

# arrival cost mean (with wrong guess for l)
# x0_bar = np.array([0.0, np.pi, 0.0, 0.0, 1])
x0_bar = np.array([0.0, 0, 0.0, 0.0, 0.2])

# solve ocp problem
status = acados_solver_ocp.solve()

if status != 0:
    raise Exception(f'acados returned status {status}.')

# get solution
for i in range(N):
    simX[i,:] = acados_solver_ocp.get(i, "x")
    simU[i,:] = acados_solver_ocp.get(i, "u")
    simY[i,:] = simX[i,:] + np.transpose(np.diag(v_stds) @ np.random.standard_normal((nx, 1)))

simX[N,:] = acados_solver_ocp.get(N, "x")
simY[N,:] = simX[N,:] + np.transpose(np.diag(v_stds) @ np.random.standard_normal((nx, 1)))

# set measurements and controls
yref_0 = np.zeros((nx + 2*nx_augmented, ))
yref_0[:nx] = simY[0, :]
yref_0[nx+nx_augmented:] = x0_bar
acados_solver_mhe.set(0, "yref", yref_0)
acados_solver_mhe.set(0, "p", simU[0,:])

# set initial guess to x0_bar
acados_solver_mhe.set(0, "x", x0_bar)

yref = np.zeros((2*nx, ))
for j in range(1, N):
    yref = np.zeros((nx+nx_augmented, ))
    yref[:nx] = simY[j, :]
    acados_solver_mhe.set(j, "yref", yref)
    acados_solver_mhe.set(j, "p", simU[j,:])

    # set initial guess to x0_bar
    acados_solver_mhe.set(j, "x", x0_bar)

# set initial guess to x0_bar
acados_solver_mhe.set(N, "x", x0_bar)

# solve mhe problem
status = acados_solver_mhe.solve()

if status != 0:
    raise Exception(f'acados returned status {status}.')

# get solution
for i in range(N):
    x_augmented = acados_solver_mhe.get(i, "x")

    simXest[i,:] = x_augmented[0:nx]
    sim_l_est[i,:] = x_augmented[nx]
    simWest[i,:] = acados_solver_mhe.get(i, "u")

x_augmented = acados_solver_mhe.get(N, "x")
simXest[N,:] = x_augmented[0:nx]
sim_l_est[N,:] = x_augmented[nx]

print('difference |x0_est - x0_bar|', np.linalg.norm(x0_bar[0:nx] - simXest[0, :]))
print('difference |x_est - x_true|', np.linalg.norm(simXest - simX))

ts = np.linspace(0, Tf, N+1)
plot_pendulum(ts, Fmax, simU, simX, simXest, simY, latexify=False)


import matplotlib.pyplot as plt

plt.figure()
plt.plot(ts, l_true*np.ones((N+1, 1)), '-')
plt.plot(ts, sim_l_est, '.-')
plt.grid()
plt.ylabel('l')
plt.legend(['true l', 'estimated l'])

plt.show()