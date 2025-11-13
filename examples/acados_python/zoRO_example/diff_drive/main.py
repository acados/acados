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
import matplotlib.pyplot as plt

from diff_drive_zoro_mpc import ZoroMPCSolver
from mpc_parameters import MPCParam, PathTrackingParam
from diff_drive_utils import plot_timings, plot_trajectory, compute_min_dis, get_results_filename, store_results, load_results, plot_timing_comparison, plot_multiple_trajectories
from nominal_path_tracking_casadi import NominalPathTrackingSolver
from track_spline import TrackSpline

N_SIM = 175

def run_closed_loop_simulation(use_custom_update: bool, feedback_optimization_mode: str, n_executions: int = 1):
    cfg_zo = MPCParam()
    cfg_zo.use_custom_update = use_custom_update
    cfg_zo.feedback_optimization_mode = feedback_optimization_mode
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
    x_init = np.array([0.0, cfg_path.v_s_0])
    x_e    = np.array([1.0, cfg_path.v_s_e])
    path_tracking_solver.solve(x_init=x_init, x_e=x_e)

    time_prep = []
    time_riccati = []
    time_zoro = []
    time_feedback = []
    time_sim = []
    time_qp = []

    for i_exec in range(n_executions):
        zoroMPC.initialized = False
        zoroMPC.acados_ocp_solver.reset()
        # closed loop mpc
        traj_zo = np.zeros((N_SIM+1, cfg_zo.nx))
        traj_zo[0,:] = path_tracking_solver.x_robot_ref[0,:]
        traj_u = np.zeros((N_SIM+1, cfg_zo.nu))
        for i_sim in range(N_SIM):
            x_ref_interp, u_ref_interp = path_tracking_solver.interpolate_reference_trajectory(robot_state=traj_zo[i_sim,:])
            u_opt, status = zoroMPC.solve(x_current=traj_zo[i_sim, :], y_ref = np.hstack((x_ref_interp, u_ref_interp)), \
                obs_position=cfg_zo.obs_pos.flatten(), obs_radius=cfg_zo.obs_radius)

            # collect timings
            if i_exec == 0:
                time_prep.append(zoroMPC.rti_phase1_t)
                time_feedback.append(zoroMPC.rti_phase2_t)
                time_riccati.append(zoroMPC.riccati_t)
                time_zoro.append(zoroMPC.zoro_t)
                time_sim.append(zoroMPC.acados_integrator_time)
                time_qp.append(zoroMPC.acados_qp_time)
            else:
                time_prep[i_sim] = min(time_prep[i_sim], zoroMPC.rti_phase1_t)
                time_feedback[i_sim] = min(time_feedback[i_sim], zoroMPC.rti_phase2_t)
                time_riccati[i_sim] = min(time_riccati[i_sim], zoroMPC.riccati_t)
                time_zoro[i_sim] = min(time_zoro[i_sim], zoroMPC.zoro_t)
                time_sim[i_sim] = min(time_sim[i_sim], zoroMPC.acados_integrator_time)
                time_qp[i_sim] = min(time_qp[i_sim], zoroMPC.acados_qp_time)

            if status != 0:
                print(f"error status={status}. Reset Solver.")
                zoroMPC.initialized = False
                zoroMPC.acados_ocp_solver.reset()
                u_opt, status = zoroMPC.solve(x_current=traj_zo[i_sim, :], \
                    y_ref = np.hstack((x_ref_interp, u_ref_interp)), \
                    obs_position=cfg_zo.obs_pos.flatten(), obs_radius=cfg_zo.obs_radius)
                if status != 0:
                    print("Failure when resolving OCP.")

            # print(i_sim, u_opt, traj_zo[i_sim,:2])
            traj_u[i_sim, :] = u_opt
            traj_zo[i_sim+1,:] = I(x0=traj_zo[i_sim, :], p=u_opt)['xf'].full().flatten()
            traj_zo[i_sim+1,:] += process_noise[i_sim,:]
            min_dist = compute_min_dis(cfg=cfg_zo, s=traj_zo[i_sim+1,:])
            if min_dist < 1e-8:
                print("collision takes place")
                return False

    total_time = [time_prep[i] + time_feedback[i] + time_riccati[i] + time_zoro[i] for i in range(len(time_prep))]
    timings = {
                   "preparation": 1e3*np.array(time_prep),
                   "integrator": 1e3*np.array(time_sim),
                   "riccati": 1e3*np.array(time_riccati),
                   "zoRO": 1e3*np.array(time_zoro),
                   "feedback": 1e3*np.array(time_feedback),
                   "QP": 1e3*np.array(time_qp),
                   "total": 1e3*np.array(total_time)
                }

    results = {
        "timings": timings,
        "trajectory": traj_zo,
        "trajectory_input": traj_u,
        "ref_trajectory": path_tracking_solver.x_robot_ref,
        "cfg_zo": cfg_zo
    }


    results_filename = get_results_filename(use_custom_update, feedback_optimization_mode, n_executions)
    store_results(results_filename, results)
    del zoroMPC



def solve_single_zoro_problem_visualize_uncertainty(feedback_optimization_mode: str="CONSTANT_FEEDBACK", converg_thr:float=1e-7, num_nominal4init:int=1):
    cfg_zo = MPCParam()
    cfg_zo.use_custom_update = True
    cfg_zo.feedback_optimization_mode = feedback_optimization_mode
    cfg_zo.zoRO_iter = 50
    zoroMPC = ZoroMPCSolver(cfg_zo, output_P_matrices=True)

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
    x_init = np.array([0.0, cfg_path.v_s_0])
    x_e = np.array([1.0, cfg_path.v_s_e])
    path_tracking_solver.solve(x_init=x_init, x_e=x_e)

    # zoro solution
    x0 = path_tracking_solver.x_robot_ref[70,:]
    x_ref_interp, u_ref_interp = path_tracking_solver.interpolate_reference_trajectory(robot_state=x0)
    u_opt, status = zoroMPC.solve(x_current=x0, y_ref = np.hstack((x_ref_interp, u_ref_interp)),
                    obs_position=cfg_zo.obs_pos.flatten(), obs_radius=cfg_zo.obs_radius, converg_thr=converg_thr, num_nominal4init=num_nominal4init)

    # get solution
    x_opt = np.zeros((cfg_zo.n_hrzn+1, cfg_zo.nx))
    for i in range(cfg_zo.n_hrzn+1):
        x_opt[i, :] = zoroMPC.acados_ocp_solver.get(i, "x")

    print(f"x_opt = {x_opt}")
    print(f"status = {status}")
    plot_trajectory(cfg_zo, x_ref_interp, x_opt,
                    P_matrices=zoroMPC.ocp.zoro_description.backoff_scaling_gamma**2 * zoroMPC.P_mats, closed_loop=False, fig_name_concat=f"OCP_{feedback_optimization_mode}")


def plot_result_trajectory(n_executions: int, use_custom_update=True, feedback_optimization_mode: str="CONSTANT_FEEDBACK"):
    results_filename = get_results_filename(use_custom_update, feedback_optimization_mode, n_executions)
    results = load_results(results_filename)
    plot_trajectory(results['cfg_zo'], results['ref_trajectory'], results['trajectory'], fig_name_concat=f"MPC_{feedback_optimization_mode}")


def closed_loop_trajectories_comparison(n_executions: int, list_feedback_optimization_mode:list):
    list_traj_label_tuple = []
    traj_ref = None
    cfg_zo = None
    for feedback_optimization_mode in list_feedback_optimization_mode:
        results_filename = get_results_filename(use_custom_update=True, feedback_optimization_mode=feedback_optimization_mode, n_executions=n_executions)
        results = load_results(results_filename)
        if traj_ref is None:
            traj_ref = results['ref_trajectory'].copy()
            cfg_zo = results['cfg_zo']
        else:
            if not np.allclose(traj_ref, results['ref_trajectory']):
                raise Exception("The reference trajectories are different.")
        
        if feedback_optimization_mode == "CONSTANT_FEEDBACK":
            label = 'ZORO, const. feedback'
        elif feedback_optimization_mode == "RICCATI_CONSTANT_COST":
            label = 'RZORO, constant Hess.'
        elif feedback_optimization_mode == "RICCATI_BARRIER_1":
            label = 'Riccati-ZORO, adaptive Hess.'
        else:
            label = feedback_optimization_mode
        list_traj_label_tuple.append((label, results['trajectory'], results['trajectory_input']))
    plot_multiple_trajectories(cfg_zo, traj_ref, list_traj_label_tuple, closed_loop=True)


def compare_results(n_executions: int, feedback_optimization_mode:int=0):
    results1 = load_results(get_results_filename(use_custom_update=True, feedback_optimization_mode=feedback_optimization_mode, n_executions=n_executions))
    results2 = load_results(get_results_filename(use_custom_update=False, feedback_optimization_mode=feedback_optimization_mode, n_executions=n_executions))
    traj_diff = results1['trajectory'] - results2['trajectory']
    error = np.max(np.abs(traj_diff))
    print(f"trajectory diff after closed loop simulation {error:.2e}")
    tol = 1.5*1e-5
    if error > tol:
        raise Exception(f"zoRO implementations differ too much, error = {error:.2e} > tol = {tol:.2e}")

def timing_comparison(n_executions: int):
    dict_results = {}

    for _, tuple in enumerate(zip([True, True, True, True, False, False], ["CONSTANT_FEEDBACK", "RICCATI_CONSTANT_COST", "RICCATI_BARRIER_1", "RICCATI_BARRIER_2", "RICCATI_CONSTANT_COST", "CONSTANT_FEEDBACK"])):
        results_filename = get_results_filename(use_custom_update=tuple[0], feedback_optimization_mode=tuple[1], n_executions=n_executions)
        results = load_results(results_filename)
        plot_timings(results['timings'], use_custom_update=tuple[0], fig_name_concat=tuple[1])
        temp_key = "zoRO-24" if tuple[0] else "zoRO-21"
        if tuple[1] != "CONSTANT_FEEDBACK":
            temp_key += "-riccati"
        dict_results[temp_key] = results

    # Compare zoro-riccati with and without custom update
    plot_timing_comparison([dict_results["zoRO-24-riccati"]['timings'], dict_results["zoRO-21-riccati"]['timings']], ['zoRO-24-riccati', 'zoRO-21-riccati'], fig_name_concat="_riccati")
    # Compare zoro with and without custom update
    plot_timing_comparison([dict_results["zoRO-24"]['timings'], dict_results["zoRO-21"]['timings']], ['zoRO-24', 'zoRO-21'], fig_name_concat="_fixedK")
    # Compare zoro-riccati to zoro with custom update
    plot_timing_comparison([dict_results["zoRO-24-riccati"]['timings'], dict_results["zoRO-24"]['timings']], ['zoRO-24-riccati', 'zoRO-24'], fig_name_concat="_zoRO-24")
    # Compare zoro-riccati to zoro without custom update
    plot_timing_comparison([dict_results["zoRO-21-riccati"]['timings'], dict_results["zoRO-21"]['timings']], ['zoRO-21-riccati', 'zoRO-21'], fig_name_concat="_zoRO-21")


if __name__ == "__main__":
    n_executions = 2
    for feedback_optimization_mode in ["CONSTANT_FEEDBACK", "RICCATI_CONSTANT_COST"]:
        run_closed_loop_simulation(use_custom_update=True, feedback_optimization_mode=feedback_optimization_mode, n_executions=n_executions)
        run_closed_loop_simulation(use_custom_update=False, feedback_optimization_mode=feedback_optimization_mode, n_executions=n_executions)
        plot_result_trajectory(n_executions=n_executions, use_custom_update=True, feedback_optimization_mode=feedback_optimization_mode)
        compare_results(n_executions=n_executions, feedback_optimization_mode=feedback_optimization_mode)

    # Feedback gain computed using riccati with sum of constant cost matrices and Hessian of tightened constraints weighted by 1/h**2
    run_closed_loop_simulation(use_custom_update=True, feedback_optimization_mode="RICCATI_BARRIER_1", n_executions=n_executions)
    plot_result_trajectory(n_executions=n_executions, use_custom_update=True, feedback_optimization_mode="RICCATI_BARRIER_1")
    # Feedback gain computed using riccati with sum of constant cost matrices and Hessian of tightened constraints weighted by -1/(h*backoff*2)
    run_closed_loop_simulation(use_custom_update=True, feedback_optimization_mode="RICCATI_BARRIER_2", n_executions=n_executions)
    plot_result_trajectory(n_executions=n_executions, use_custom_update=True, feedback_optimization_mode="RICCATI_BARRIER_2")

    closed_loop_trajectories_comparison(n_executions=n_executions, list_feedback_optimization_mode=["CONSTANT_FEEDBACK", "RICCATI_BARRIER_1"])
    closed_loop_trajectories_comparison(n_executions=n_executions, list_feedback_optimization_mode=["CONSTANT_FEEDBACK", "RICCATI_CONSTANT_COST", "RICCATI_BARRIER_1", "RICCATI_BARRIER_2"])

    timing_comparison(n_executions=n_executions)

    solve_single_zoro_problem_visualize_uncertainty(feedback_optimization_mode="CONSTANT_FEEDBACK")
    solve_single_zoro_problem_visualize_uncertainty(feedback_optimization_mode="RICCATI_CONSTANT_COST")
    solve_single_zoro_problem_visualize_uncertainty(feedback_optimization_mode="RICCATI_BARRIER_1")
    solve_single_zoro_problem_visualize_uncertainty(feedback_optimization_mode="RICCATI_BARRIER_2")
