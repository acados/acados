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

from ekf import build_ekf
from setup_ocp_solver import setup_ocp_solver, setup_integrator
from utils import load_trajectory, wind_profile
from visualization import plot_3d_animated, plot_3d_static, plot_controls, plot_estimations, plot_states

# %% ### Setup ###
np.random.seed(42)

RTI = True
dt = 0.1
N_horizon = 40
ekf_substeps = 5   # simulator and estimator steps per control interval
Nsim_oribts = 2.0  # number of orbits to simulate
wind_phases = [np.array([3.0, 1.0, 0.0]), np.array([1.0, 3.0, 0.0]), np.array([3.0, 1.0, 0.0])]
trajectory_file = 'data/trajectory_tilt20deg_r50m.csv'

# noise covariances and initial estimation covariance
W = np.diag([0.2]*3 + [0.1]*3 + [0.1]*9 + [0.1]*3 + [0.01]*2 + [0.01]*2 + [0.5, 0.5, 0.1]) ** 2 * (dt / ekf_substeps)
V = np.diag([1]*3 + [0.6]*3 + [0.02]*3 + [0.5] + [0.05]*3) ** 2
P_0 = np.diag([1.]*3 + [0.5]*3 + [0.01]*9 + [0.05]*3 + [0.1]*2 + [0.05]*2 + [1, 1, .1]) ** 2

x_ref, u_ref, orbit_duration = load_trajectory(trajectory_file, dt=dt, number_of_orbits=Nsim_oribts+1)
x0 = x_ref[0]

Tf = dt * N_horizon

ocp, ocp_solver = setup_ocp_solver(x0=x0, N_horizon=N_horizon, Tf=Tf, RTI=RTI, soft_constraints=True)
integrator = setup_integrator(ocp)
integrator.set("T", dt / ekf_substeps)

ekf, measure = build_ekf(ocp.model, integrator, W=W, V=V)

nx = ocp_solver.acados_ocp.dims.nx
nu = ocp_solver.acados_ocp.dims.nu

Nsim = int(orbit_duration / dt * Nsim_oribts)
simX = np.zeros((Nsim+1, nx))
simU = np.zeros((Nsim, nu))
simX_pred = np.zeros((Nsim + 1, nx, N_horizon + 1))
simU_pred = np.zeros((Nsim + 1, nu, N_horizon)) 

simX_est = np.zeros((Nsim+1, nx))               # estimated states
simWind_est = np.zeros((Nsim+1, 3))             # estimated wind
simP_est = np.zeros((Nsim+1, *P_0.shape))       # covariance of [x, wind]
simX_est[0] = x0
simP_est[0] = P_0

simX[0, :] = x0

t_preparation = np.zeros((Nsim))
t_feedback = np.zeros((Nsim))

# wind disturbance: starts at rest, 1-cosine gust into each phase
t_wind = dt * np.arange(Nsim + N_horizon)
wind = wind_profile(wind_phases, gust_duration=3, t=t_wind)

# %% ### closed loop ###
for i in range(Nsim):
    if i%10 == 0 or i == Nsim-1:
        print(f"{i}/{Nsim-1}", end='\r')

    ### preparation phase ###
    # set reference and parameter values
    for k in range(N_horizon):
        yref_k = np.concatenate([x_ref[i+k], u_ref[i+k]])
        ocp_solver.cost_set(k, "yref", yref_k)
        ocp_solver.set(k, "p", simWind_est[i, :])
    ocp_solver.cost_set(N_horizon, "yref", x_ref[i+N_horizon])
    ocp_solver.set(N_horizon, "p", simWind_est[i, :])

    if RTI:
        ocp_solver.options_set('rti_phase', 1)
        status = ocp_solver.solve()
        t_preparation[i] = ocp_solver.get_stats('time_tot')

    # set initial state
    ocp_solver.set(0, "lbx", simX_est[i, :])
    ocp_solver.set(0, "ubx", simX_est[i, :])

    ### feedback phase ###
    if RTI: ocp_solver.options_set('rti_phase', 2)
    status = ocp_solver.solve()
    t_feedback[i] = ocp_solver.get_stats('time_tot')

    simU[i, :] = ocp_solver.get(0, "u")

    # simulation and estimation loop
    simX[i+1] = simX[i]
    simX_est[i+1], simWind_est[i+1], simP_est[i+1] = simX_est[i], simWind_est[i], simP_est[i]
    for _ in range(ekf_substeps):
        simX[i+1, :] = integrator.simulate(x=simX[i+1, :], u=simU[i, :], p=wind[i])
        meassurement = measure(simX[i+1, :], wind[i])
        simX_est[i+1], simWind_est[i+1], simP_est[i+1] = ekf(simX_est[i+1], simWind_est[i+1], simP_est[i+1], meassurement, simU[i])

    simX_pred[i] = np.array([ocp_solver.get(k, "x") for k in range(N_horizon + 1)]).T
    simU_pred[i] = np.array([ocp_solver.get(k, "u") for k in range(N_horizon)]).T

# evaluate timings
print(f"\nTiming statistics (ms):")
print(f"  Preparation:  min={1e3*np.min(t_preparation):7.3f}   median={1e3*np.median(t_preparation):7.3f}   max={1e3*np.max(t_preparation):7.3f}")
print(f"  Feedback:     min={1e3*np.min(t_feedback):7.3f}   median={1e3*np.median(t_feedback):7.3f}   max={1e3*np.max(t_feedback):7.3f}")

# %% ### visualization
t = dt * np.arange(simX.shape[0])
plot_3d_static(X=simX, x_ref=x_ref, wind=wind, num_planes=9)
plot_controls(U=simU, X=simX, t=t, u_ref=u_ref, x_ref=x_ref)
plot_states(X=simX, t=t, x_ref=x_ref, wind=wind)
plot_estimations(X=simX, X_est=simX_est, t=t, P=simP_est, wind=wind, wind_est=simWind_est)
plt.show()

# %% additional 3d animation
anim, fig_anim = plot_3d_animated(X=simX, U=simU, t=t, X_pred=simX_pred, U_pred=simU_pred, x_ref=x_ref, u_ref=u_ref, wind=wind, wind_est=simWind_est, P=simP_est)
plt.show()
