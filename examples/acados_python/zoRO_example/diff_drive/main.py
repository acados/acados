import sys
import os
local_path = os.path.dirname(os.path.abspath(__file__))
mpc_source_dir = os.path.join(local_path, '..')
sys.path.append(mpc_source_dir)

import numpy as np
import matplotlib.pyplot as plt
from time import process_time

import casadi

from diff_drive_zoro_mpc import ZoroMPCSolver
from mpc_parameters import MPCParam, PathTrackingParam
from diff_drive_utils import plot_timings, compute_min_dis
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

    # plot trajectory
    fig = plt.figure(1)
    plt.rcParams['font.size'] = '16'
    ax = fig.add_subplot(1,1,1)
    for idx_obs in range(cfg_zo.num_obs):
        circ = plt.Circle(cfg_zo.obs_pos[idx_obs,:], cfg_zo.obs_radius[idx_obs], \
            edgecolor="red", facecolor=(1,0,0,.5))
        ax.add_artist(circ)
    plt.plot(path_tracking_solver.x_robot_ref[:, 0], path_tracking_solver.x_robot_ref[:, 1], c='m', label='ref')
    plt.plot(traj_zo[:, 0], traj_zo[:, 1], c='b', label='opt sqp')
    plt.legend()
    plt.show()


if __name__ == "__main__":
    main()