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
from export_mhe_ode_model import export_mhe_ode_model

from export_ocp_solver import export_ocp_solver
from export_mhe_solver import export_mhe_solver

import numpy as np
from utils import plot_pendulum

# general

Tf = 1.0
N = 20
h = Tf/N
Fmax = 80

# ocp model and solver
model = export_pendulum_ode_model()

nx = model.x.size()[0]
nu = model.u.size()[0]

Q_ocp = np.diag([1e3, 1e3, 1e-2, 1e-2])
R_ocp = 1e-2 *np.eye(1)

acados_solver_ocp = export_ocp_solver(model, N, h, Q_ocp, R_ocp, Fmax)

# mhe model and solver
model_mhe = export_mhe_ode_model()

nx = model_mhe.x.size()[0]
nw = model_mhe.u.size()[0]
ny = nx

Q0_mhe = 100*np.eye((nx))
Q_mhe  = 0.1*np.eye(nx)
R_mhe  = 0.1*np.eye(nx)
# Q_mhe = np.zeros((nx, nx))
# Q0_mhe = np.diag([0.01, 1, 1, 1])
# R_mhe  = np.diag([0.1, 10, 10, 10])
acados_solver_mhe = export_mhe_solver(model_mhe, N, h, Q_mhe, Q0_mhe, R_mhe)

# simulation
v_stds = [0.1, 0.01, 0.01, 0.01]
v_stds = [0, 0, 0, 0]

simX = np.ndarray((N+1, nx))
simU = np.ndarray((N, nu))
simY = np.ndarray((N+1, nx))

simXest = np.zeros((N+1, nx))
simWest = np.zeros((N, nx))

# arrival cost mean
x0_bar = np.array([0.0, np.pi, 0.0, 0.0])

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
yref_0 = np.zeros((3*nx, ))
yref_0[:nx] = simY[0, :]
yref_0[2*nx:] = x0_bar
acados_solver_mhe.set(0, "yref", yref_0)
acados_solver_mhe.set(0, "p", simU[0,:])
#acados_solver_mhe.set(0, "x", simX[0,:])

yref = np.zeros((2*nx, ))
for j in range(1,N):
    yref[:nx] = simY[j, :]
    acados_solver_mhe.set(j, "yref", yref)
    acados_solver_mhe.set(j, "p", simU[j,:])
    # acados_solver_mhe.set(j, "x", simX[j,:])


# solve mhe problem
status = acados_solver_mhe.solve()

if status != 0 and status != 2:
    raise Exception(f'acados returned status {status}.')

# get solution
for i in range(N):
    simXest[i,:] = acados_solver_mhe.get(i, "x")
    simWest[i,:] = acados_solver_mhe.get(i, "u")

simXest[N, :] = acados_solver_mhe.get(N, "x")

print('difference |x0_est - x0_bar|', np.linalg.norm(x0_bar - simXest[0, :]))
print('difference |x_est - x_true|', np.linalg.norm(simXest - simX))

plot_pendulum(np.linspace(0, Tf, N+1), Fmax, simU, simX, simXest, simY, latexify=False)
