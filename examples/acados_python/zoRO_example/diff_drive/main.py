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

import numpy as np
import casadi

from diff_drive_zoro_mpc import ZoroMPCSolver
from mpc_parameters import MPCParam, PathTrackingParam
from diff_drive_utils import plot_timings, plot_trajectory, compute_min_dis, get_results_filename, store_results, load_results, plot_timing_comparison
from nominal_path_tracking_casadi import NominalPathTrackingSolver
from track_spline import TrackSpline

N_SIM = 175

def run_closed_loop_simulation(use_custom_update: bool, n_executions: int = 1):
    cfg_zo = MPCParam()
    cfg_zo.use_custom_update = use_custom_update
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

    for i_exec in range(n_executions):
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

            # print(i_sim, u_opt, traj_zo[i_sim,:2])
            traj_zo[i_sim+1,:] = I(x0=traj_zo[i_sim, :], p=u_opt)['xf'].full().flatten()
            traj_zo[i_sim+1,:] += process_noise[i_sim,:]
            min_dist = compute_min_dis(cfg=cfg_zo, s=traj_zo[i_sim+1,:])
            if min_dist < 1e-8:
                print("collision takes place")
                return False

    total_time = [time_prep[i] + time_feedback[i] + time_prop[i] for i in range(len(time_prep))]
    timings = {
                   "preparation": 1e3*np.array(time_prep),
                   "integrator": 1e3*np.array(time_sim),
                   "propagation": 1e3*np.array(time_prop),
                   "feedback": 1e3*np.array(time_feedback),
                   "QP": 1e3*np.array(time_qp),
                   "total": 1e3*np.array(total_time)
                }

    results = {
        "timings": timings,
        "trajectory": traj_zo,
        "ref_trajectory": path_tracking_solver.x_robot_ref,
        "cfg_zo": cfg_zo
    }


    results_filename = get_results_filename(use_custom_update, n_executions)
    store_results(results_filename, results)
    del zoroMPC


def plot_result_trajectory(n_executions: int, use_custom_update=True):
    results_filename = get_results_filename(use_custom_update, n_executions)
    results = load_results(results_filename)
    plot_trajectory(results['cfg_zo'], results['ref_trajectory'], results['trajectory'])

def plot_result_timings(n_executions: int, use_custom_update=True):
    results_filename = get_results_filename(use_custom_update, n_executions)
    results = load_results(results_filename)
    plot_timings(results['timings'], use_custom_update)

def compare_results(n_executions: int):
    results1 = load_results(get_results_filename(use_custom_update=True, n_executions=n_executions))
    results2 = load_results(get_results_filename(use_custom_update=False, n_executions=n_executions))
    traj_diff = results1['trajectory'] - results2['trajectory']
    error = np.max(np.abs(traj_diff))
    print(f"trajectory diff after closed loop simulation {error:.2e}")
    tol = 1e-5
    if error > tol:
        raise Exception(f"zoRO implementations differ too much, error = {error:.2e} > tol = {tol:.2e}")


def plot_result_timing_comparison(n_executions: int):
    fast_timings = load_results(get_results_filename(use_custom_update=True, n_executions=n_executions))['timings']
    slow_timings = load_results(get_results_filename(use_custom_update=False, n_executions=n_executions))['timings']
    plot_timing_comparison([fast_timings, slow_timings], ['zoRO-24', 'zoRO-21'])

def timing_comparison():
    n_executions = 50
    # run_closed_loop_simulation(use_custom_update=False, n_executions=n_executions)
    # run_closed_loop_simulation(use_custom_update=True, n_executions=n_executions)
    plot_result_timings(n_executions=n_executions, use_custom_update=True)
    plot_result_timings(n_executions=n_executions, use_custom_update=False)

    plot_result_timing_comparison(n_executions=n_executions)


if __name__ == "__main__":
    n_executions = 2
    run_closed_loop_simulation(use_custom_update=True, n_executions=n_executions)
    run_closed_loop_simulation(use_custom_update=False, n_executions=n_executions)
    compare_results(n_executions=n_executions)

    # plot_result_trajectory(n_executions=n_executions, use_custom_update=True)
    # timing_comparison()
