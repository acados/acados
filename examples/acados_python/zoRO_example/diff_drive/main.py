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

import os

import numpy as np
import matplotlib.pyplot as plt

import casadi

from diff_drive_zoro_mpc import ZoroMPCSolver
from mpc_parameters import MPCParam, PathTrackingParam
from diff_drive_utils import plot_timings, plot_trajectory, compute_min_dis
from nominal_path_tracking_casadi import NominalPathTrackingSolver
from track_spline import TrackSpline

N_EXEC = 5
N_SIM = 175

def main():
    cfg_zo = MPCParam()
    zoroMPC = ZoroMPCSolver(cfg_zo)

    # Differential equation of the model
    x_int = casadi.SX.sym('x_int', cfg_zo.nx)  # x, y, theta, v, omega
    u_int = casadi.SX.sym('u_int', cfg_zo.nu)  # a, alpha
    f_x = casadi.vertcat(x_int[3]*casadi.cos(x_int[2]),
                        x_int[3]*casadi.sin(x_int[2]),
                        x_int[4],
                        u_int[0],
                        u_int[1])
    # Create an integrator
    dae = {'x': x_int, 'p': u_int, 'ode': f_x}
    opts = {'tf': cfg_zo.delta_t} # interval length
    I = casadi.integrator('I', 'rk', dae, opts)

    # Process Noise
    np.random.seed(1)
    process_noise = np.random.multivariate_normal(np.zeros((cfg_zo.nw,)), cfg_zo.W_mat, N_SIM)

    # Reference trajectory
    spline_control_points = np.array([[-2.0, 2.0],
                               [0.0, 0.0],
                               [2.0, 0.0],
                               [4.0, 0.0],
                               [6.0, 0.0],
                               [8.0, 2.0],
                               [8.0, 4.0],
                               [8.0, 6.0],
                               [8.0, 8.0]])
    spline_seg_length = np.array([0.0, np.pi, 2.0, 2.0, 2.0, np.pi, 2.0, 2.0, 2.0])
    track_spline = TrackSpline(spline_control_points, spline_seg_length)
    cfg_path = PathTrackingParam()
    path_tracking_solver = NominalPathTrackingSolver(track_spline, cfg_path=cfg_path, cfg_traj=cfg_zo)
    x_init = np.array([0.0, cfg_path._v_s_0])
    x_e    = np.array([1.0, cfg_path._v_s_e])
    path_tracking_solver.solve(x_init=x_init, x_e=x_e)

    time_prep = []
    time_prop = []
    time_feedback = []
    time_sim = []
    time_qp = []

    for i_exec in range(N_EXEC):
        zoroMPC.initialized = False
        zoroMPC.acados_ocp_solver.reset()
        # closed loop mpc
        traj_zo = np.zeros((N_SIM+1, cfg_zo.nx))
        traj_zo[0,:] = path_tracking_solver.x_robot_ref[0,:]
        for i_sim in range(N_SIM):
            x_ref_interp, u_ref_interp = path_tracking_solver.interpolate_reference_trajectory(robot_state=traj_zo[i_sim,:])
            u_opt, status = zoroMPC.solve(x_current=traj_zo[i_sim, :], y_ref = np.hstack((x_ref_interp, u_ref_interp)), \
                obs_position=cfg_zo.obs_pos.flatten(), obs_radius=cfg_zo.obs_radius)

            # collect timings
            if i_exec == 0:
                time_prep.append(zoroMPC.rti_phase1_t)
                time_feedback.append(zoroMPC.rti_phase2_t)
                time_prop.append(zoroMPC.propagation_t)
                time_sim.append(zoroMPC.acados_integrator_time)
                time_qp.append(zoroMPC.acados_qp_time)
            else:
                time_prep[i_sim] = min(time_prep[i_sim], zoroMPC.rti_phase1_t)
                time_feedback[i_sim] = min(time_feedback[i_sim], zoroMPC.rti_phase2_t)
                time_prop[i_sim] = min(time_prop[i_sim], zoroMPC.propagation_t)
                time_sim[i_sim] = min(time_sim[i_sim], zoroMPC.acados_integrator_time)
                time_qp[i_sim] = min(time_qp[i_sim], zoroMPC.acados_qp_time)

            if status != 0:
                print('error status=',status,'Reset Solver')
                zoroMPC.initialized = False
                zoroMPC.acados_ocp_solver.reset()
                u_opt, status = zoroMPC.solve(x_current=traj_zo[i_sim, :], \
                    y_ref = np.hstack((x_ref_interp, u_ref_interp)), \
                    obs_position=cfg_zo.obs_pos.flatten(), obs_radius=cfg_zo.obs_radius)

            print(i_sim, u_opt, traj_zo[i_sim,:2])
            traj_zo[i_sim+1,:] = I(x0=traj_zo[i_sim, :], p=u_opt)['xf'].full().flatten()
            traj_zo[i_sim+1,:] += process_noise[i_sim,:]
            min_dist = compute_min_dis(cfg=cfg_zo, s=traj_zo[i_sim+1,:])
            if min_dist < 1e-8:
                print("collision take place")
                return False


    total_time = [time_prep[i] + time_feedback[i] + time_prop[i] for i in range(len(time_prep))]
    timing_dict = {
                   "integrator": 1e3*np.array(time_sim),
                   "preparation": 1e3*np.array(time_prep),
                   "QP": 1e3*np.array(time_qp),
                   "feedback": 1e3*np.array(time_feedback),
                   "propagation": 1e3*np.array(time_prop),
                   "total": 1e3*np.array(total_time)
                }
    plot_timings(timing_dict)
    plot_trajectory(cfg_zo, path_tracking_solver.x_robot_ref, traj_zo)


if __name__ == "__main__":
    main()