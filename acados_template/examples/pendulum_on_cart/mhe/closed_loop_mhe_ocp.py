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

# authors: Katrin Baumgaertner, Jonathan Frey

import sys
sys.path.insert(0, '../common')

from pendulum_model import export_pendulum_ode_model
from export_mhe_ode_model import export_mhe_ode_model

from export_ocp_solver import export_ocp_solver
from export_mhe_solver import export_mhe_solver

from export_ode_mhe_integrator import export_ode_mhe_integrator
from utils import plot_pendulum
import numpy as np
from scipy.linalg import block_diag


Tf_ocp = 1.0
N_ocp = 20

Ts = Tf_ocp/N_ocp # time step

Tf_mhe = 0.5*Tf_ocp
N_mhe = int(Tf_mhe/Ts)

u_max = 80

# state and measurement noise
v_stds_mhe = np.array([0.1, 0.1, .5, 0.3]) # measurement noise stds
w_stds_mhe = np.array([0.01, 0.001, 0.001, 0.001]) # state noise stds

v_stds_plant = .8 * np.array([0.1, 0.1, .5, 0.3])
w_stds_plant = np.zeros((4,))

V = np.diag(v_stds_plant)
W = np.diag(w_stds_plant)

# ocp model and solver
model = export_pendulum_ode_model()

nx = model.x.size()[0]
nu = model.u.size()[0]

Q_ocp = np.diag([1e3, 1e3, 1e-2, 1e-2])
R_ocp = 1 * np.eye(1)

acados_solver_ocp = export_ocp_solver(model, N_ocp, Ts, Q_ocp, R_ocp, Fmax=u_max)

# mhe model and solver
model_mhe = export_mhe_ode_model()

nw = model_mhe.u.size()[0]
ny = nx

# inverse covariances, R_mhe has to be scaled with h
Q_mhe = np.diag(1/w_stds_mhe)
R_mhe = 1/Ts*np.diag(1/v_stds_mhe)

# arrival cost weighting
Q0_mhe = 0.01*Q_mhe

acados_solver_mhe = export_mhe_solver(model_mhe, N_mhe, Ts, Q_mhe, Q0_mhe, R_mhe)

# integrator/plant
plant = export_ode_mhe_integrator(model_mhe, Ts)

# simulation
Nsim = 100

simX = np.zeros((Nsim+1, nx))
simU = np.zeros((Nsim, nu))
simY = np.zeros((Nsim+1, nx))

simXest = np.zeros((Nsim+1, nx))
simWest = np.zeros((Nsim+1, nx))

# arrival cost mean & initial state
x0_plant = np.array([0.1, np.pi + 0.5, -0.05, 0.05])
x0_bar = np.array([0.0, np.pi, 0.0, 0.0])

u0 = np.zeros((nu,))

# initial state
simX[0,:] = x0_plant

# init solvers
for i in range(N_mhe):
    acados_solver_mhe.set(i, "x", x0_bar)

# simulate for N_mhe steps with zero input
for i in range(N_mhe):
    
    # measurement
    simY[i,:] = simX[i,:] + (V @ np.random.standard_normal((nx,1))).T

    # simulate one step 
    w = W @ np.random.standard_normal((nx,))
    plant.set("u", w)
    plant.set("p", u0)
    plant.set("x", simX[i,:])

    # solve
    status = plant.solve()

    if status != 0:
        raise Exception('integrator returned status {} in step {}. Exiting.'.format(status, i))

    # get solution
    simX[i+1,:] = plant.get("x")
   

# reference for mhe
yref = np.zeros((2*nx, ))
yref_0 = np.zeros((3*nx, ))

# closed loop
for i in range(N_mhe, Nsim):

    ### estimation ###
    k = i - N_mhe
    
    # set measurements
    yref_0[:nx] = simY[k, :]
    yref_0[2*nx:] = x0_bar

    acados_solver_mhe.set(0, "yref", yref_0)
    # set controls
    acados_solver_mhe.set(0, "p", simU[k,:])

    for j in range(1, N_mhe):
        # set measurements
        yref[:nx] = simY[k+j, :]
        acados_solver_mhe.set(j, "yref", yref)
        # set controls
        acados_solver_mhe.set(j, "p", simU[k+j,:])

    status = acados_solver_mhe.solve()

    if status != 0:
        raise Exception('estimator returned status {} in step {}. Exiting.'.format(status, i))
    simXest[i,:] = acados_solver_mhe.get(N_mhe, "x")

    # update arrival cost
    x0_bar = acados_solver_mhe.get(1, "x")

    ### control ###
    # update initial condition of ocp solver
    acados_solver_ocp.set(0, "lbx", simXest[i, :])
    acados_solver_ocp.set(0, "ubx", simXest[i, :])

    status = acados_solver_ocp.solve()
    # acados_solver_ocp.print_statistics()
    if status != 0:
        raise Exception('controller returned status {} in step {}. Exiting.'.format(status, i))

    simU[i:, ] = acados_solver_ocp.get(0, "u")

    ### simulation ###
    # measurement
    simY[i,:] = simX[i, :] + (V @ np.random.standard_normal((nx,1))).T

    w = W @ np.random.standard_normal((nx,))
    plant.set("u", w)
    plant.set("p", simU[i,:])
    plant.set("x", simX[i,:])

    status = plant.solve()
    if status != 0:
        raise Exception('integrator returned status {} in step {}. Exiting.'.format(status, i))

    simX[i+1,:] = plant.get("x")

# plot
print('estimation error p', np.linalg.norm(simX[:, 0] - simXest[:, 0]))
print('estimation error theta', np.linalg.norm(simX[:, 1] - simXest[:, 1]))
print('estimation error v', np.linalg.norm(simX[:, 2] - simXest[:, 2]))
print('estimation error dtheta', np.linalg.norm(simX[:, 3] - simXest[:, 3]))

plot_pendulum(np.linspace(0, Ts*Nsim, Nsim+1), u_max, simU, simX, simXest[N_mhe:, :], simY)
