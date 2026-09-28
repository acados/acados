#
# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# Licensed under the 2-Clause BSD License.

import sys
sys.path.insert(0, '../common')

from pendulum_model import export_pendulum_ode_model
from export_mhe_ode_model_with_param import export_mhe_ode_model_with_param

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
v_stds = np.array([0.4, 0.2, 0.2, 0.2])

# ocp model and solver
model = export_pendulum_ode_model()

nx = model.x.rows()
nu = model.u.rows()

Q_ocp = np.diag([1e3, 1e3, 1e-2, 1e-2])
R_ocp = 1e-2 *np.eye(1)

acados_solver_ocp = export_ocp_solver(model, N, h, Q_ocp, R_ocp, Fmax)

# mhe model and solver
model_mhe = export_mhe_ode_model_with_param()

nx_augmented = model_mhe.x.rows()
nw = model_mhe.u.rows()
ny = nx


Q0_mhe = np.diag([0.1, 0.01, 0.01, 0.01, 10])
Q_mhe  = 10.*np.diag([0.1, 0.01, 0.01, 0.01])
R_mhe  = 0.1*np.diag(1./v_stds**2)

acados_solver_mhe = export_mhe_solver_with_param(model_mhe, N, h, Q_mhe, Q0_mhe, R_mhe)

# simulation

simX = np.zeros((N+1, nx))
simU = np.zeros((N, nu))
simY = np.zeros((N+1, nx))

simXest = np.zeros((N+1, nx))
simWest = np.zeros((N, nx))
sim_l_est = np.zeros((N+1, 1))

# arrival cost mean (with wrong guess for l)
x0_bar = np.array([0.0, np.pi, 0.0, 0.0, 1])

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
yref_0 = np.zeros((2*nx + nx_augmented, ))
yref_0[:nx] = simY[0, :]
yref_0[2*nx:] = x0_bar
acados_solver_mhe.set(0, "yref", yref_0)
acados_solver_mhe.set(0, "p", simU[0,:])

# set initial guess to x0_bar
acados_solver_mhe.set(0, "x", x0_bar)

yref = np.zeros((2*nx, ))
for j in range(1, N):
    # set measurements and controls
    yref[:nx] = simY[j, :]
    acados_solver_mhe.set(j, "yref", yref)
    acados_solver_mhe.set(j, "p", simU[j,:])

    # set initial guess to x0_bar
    acados_solver_mhe.set(j, "x", x0_bar)

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
print('difference |l_est - l_true|', np.abs(sim_l_est[0] - l_true))

ts = np.linspace(0, Tf, N+1)
plot_pendulum(ts, Fmax, simU, simX, simXest, simY, latexify=False)
