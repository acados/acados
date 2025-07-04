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


def main():

    np.random.seed(42)

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

    V_mat = np.diag(v_stds_plant)
    W_mat = np.diag(w_stds_plant)

    # ocp model and solver
    model = export_pendulum_ode_model()

    nx = model.x.rows()
    nu = model.u.rows()

    Q_ocp = np.diag([1e3, 1e3, 1e-2, 1e-2])
    R_ocp = 1 * np.eye(1)

    acados_ocp_solver = export_ocp_solver(model, N_ocp, Ts, Q_ocp, R_ocp, Fmax=u_max)

    # mhe model and solver
    model_mhe = export_mhe_ode_model()

    # inverse covariances
    Q_mhe = np.diag(1/w_stds_mhe**2)
    R_mhe = np.diag(1/v_stds_mhe**2)

    # arrival cost weighting
    Q0_mhe = 0.005*Q_mhe

    acados_mhe_solver = export_mhe_solver(model_mhe, N_mhe, Ts, Q_mhe, Q0_mhe, R_mhe)

    # integrator/plant
    plant = export_ode_mhe_integrator(model_mhe, Ts)

    # simulation
    Nsim = 100

    simX = np.zeros((Nsim+1, nx))
    simU = np.zeros((Nsim, nu))
    simY = np.zeros((Nsim+1, nx))

    simXest = np.zeros((Nsim+1, nx))

    # arrival cost mean & initial state
    x0_plant = np.array([0.1, np.pi + 0.5, -0.05, 0.05])
    x0_bar = np.array([0.0, np.pi, 0.0, 0.0])

    u0 = np.zeros((nu,))

    # initial state
    simX[0,:] = x0_plant

    # initialize MHE solver
    for i in range(N_mhe):
        acados_mhe_solver.set(i, "x", x0_bar)

    # simulate for N_mhe steps with zero input
    for i in range(N_mhe):

        # measurement
        simY[i,:] = simX[i,:] + (V_mat @ np.random.standard_normal((nx,1))).T

        # simulate
        w = W_mat @ np.random.standard_normal((nx,))
        simX[i+1,:] = plant.simulate(x=simX[i,:], u=w, p=u0)

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

        acados_mhe_solver.set(0, "yref", yref_0)
        # set controls
        acados_mhe_solver.set(0, "p", simU[k,:])

        for j in range(1, N_mhe):
            # set measurements
            yref[:nx] = simY[k+j, :]
            acados_mhe_solver.set(j, "yref", yref)
            # set controls
            acados_mhe_solver.set(j, "p", simU[k+j,:])

        status = acados_mhe_solver.solve()

        if status != 0:
            raise Exception(f'estimator returned status {status} in step {i}.')
        simXest[i,:] = acados_mhe_solver.get(N_mhe, "x")

        # update arrival cost
        x0_bar = acados_mhe_solver.get(1, "x")

        ### control ###
        simU[i:, ] = acados_ocp_solver.solve_for_x0(simXest[i, :])

        ### simulation ###
        # measurement
        simY[i,:] = simX[i, :] + (V_mat @ np.random.standard_normal((nx,1))).T
        w = W_mat @ np.random.standard_normal((nx,))
        simX[i+1,:] = plant.simulate(x=simX[i,:], u=w, p=simU[i,:])

    # plot
    print('estimation error p', np.linalg.norm(simX[:, 0] - simXest[:, 0]))
    print('estimation error theta', np.linalg.norm(simX[:, 1] - simXest[:, 1]))
    print('estimation error v', np.linalg.norm(simX[:, 2] - simXest[:, 2]))
    print('estimation error dtheta', np.linalg.norm(simX[:, 3] - simXest[:, 3]))

    plot_pendulum(np.linspace(0, Ts*Nsim, Nsim+1), u_max, simU, simX, simXest[N_mhe:, :], simY)

if __name__ == '__main__':
    main()
