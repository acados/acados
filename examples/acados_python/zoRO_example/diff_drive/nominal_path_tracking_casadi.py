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

import sys
import os
local_path = os.path.dirname(os.path.abspath(__file__))
mpc_source_dir = os.path.join(local_path, '..')
sys.path.append(mpc_source_dir)

import numpy as np
import casadi
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

from track_spline import TrackSpline
from mpc_parameters import PathTrackingParam, MPCParam

class NominalPathTrackingSolver:
    def __init__(self, track_spline:TrackSpline, cfg_path: PathTrackingParam, cfg_traj: MPCParam) -> None:
        self._track_spline = track_spline
        self._path_cfg = cfg_path
        self._traj_cfg = cfg_traj
        self._dt_sol = 0.1
        self._x_sol  = np.zeros((cfg_path.nx, cfg_path.n_hrzn+1))
        self._u_sol  = np.zeros((cfg_path.nu, cfg_path.n_hrzn  ))
        self._x_robot_ref  = np.zeros((cfg_path.n_hrzn+1, cfg_traj.nx, ))
        self._u_robot_ref  = np.zeros((cfg_path.n_hrzn,   cfg_traj.nu  ))
        self._spline_timestamps   = np.full((cfg_path.n_hrzn+1, ), np.nan)
        self._setup_MPC()

    def _setup_MPC(self) -> None:
        self._g    = []
        self._lb_g = []
        self._ub_g = []
        self._J = 0.
        self._define_opt_variables()
        self._setup_system_dynamics_constraints()
        self._setup_state_input_constraints()
        self._setup_time_optimal_cost()

        ocp = {
            'f': self._J,
            'x': casadi.vertcat(casadi.vec(self._x_traj), casadi.vec(self._u_traj), self._dt_sym),
            'p': casadi.vertcat(self._param_x_0, self._param_x_e),
            'g': casadi.vertcat(*self._g),
        }

        options = {
            "print_time": False,
            "error_on_fail": True,
            "ipopt.check_derivatives_for_naninf": "yes",
            "ipopt.print_level": 3,
        }
        self._ocp_solver = casadi.nlpsol("solver", "ipopt", ocp, options)

    def _define_opt_variables(self) -> None:
        self._x_traj = casadi.MX.sym('x_sym', self._path_cfg.nx, self._path_cfg.n_hrzn+1)
        self._u_traj = casadi.MX.sym('u_sym', self._path_cfg.nu, self._path_cfg.n_hrzn)
        self._dt_sym = casadi.MX.sym('dt_sym')
        self._param_x_0 = casadi.MX.sym('param_x_0', self._path_cfg.nx)
        self._param_x_e = casadi.MX.sym('param_x_e', self._path_cfg.nx)

    def _setup_system_dynamics_constraints(self) -> None:
        # Initial state
        self._g.append(self._x_traj[:,0]- self._param_x_0)
        self._lb_g.append(np.zeros((self._path_cfg.nx, )))
        self._ub_g.append(np.zeros((self._path_cfg.nx, )))
        # Middle states
        for i_hrzn in range(self._path_cfg.n_hrzn):
            temp_s_next = casadi.vertcat(self._x_traj[0, i_hrzn] + self._x_traj[1, i_hrzn] * self._dt_sym \
                                             + self._u_traj[0, i_hrzn] * self._dt_sym**2 / 2, \
                                         self._x_traj[1, i_hrzn] + self._u_traj[0, i_hrzn] * self._dt_sym)
            self._g.append(self._x_traj[:, i_hrzn+1] - temp_s_next)
            self._lb_g.append(np.zeros((self._path_cfg.nx, )))
            self._ub_g.append(np.zeros((self._path_cfg.nx, )))
        # Terminal State
        self._g.append(self._x_traj[:,-1]- self._param_x_e)
        self._lb_g.append([0., 0.])
        self._ub_g.append([0., 0.])

    def _setup_state_input_constraints(self) -> None:
        for i_hrzn in range(0, self._path_cfg.n_hrzn):
            v_i     = self._track_spline.ca_fun_v_t    (self._x_traj[0, i_hrzn], self._x_traj[1, i_hrzn])
            a_i     = self._track_spline.ca_fun_a_t    (self._x_traj[0, i_hrzn], self._x_traj[1, i_hrzn], self._u_traj[0, i_hrzn])
            omega_i = self._track_spline.ca_fun_omega_t(self._x_traj[0, i_hrzn], self._x_traj[1, i_hrzn])
            alpha_i = self._track_spline.ca_fun_alpha_t(self._x_traj[0, i_hrzn], self._x_traj[1, i_hrzn], self._u_traj[0, i_hrzn])
            self._g.append(casadi.vertcat(v_i, a_i, omega_i, alpha_i))
            self._lb_g.append([self._traj_cfg.min_forward_velocity, self._traj_cfg.min_forward_acceleration, -self._traj_cfg.max_angular_velocity, -self._traj_cfg.max_angular_acceleration])
            self._ub_g.append([self._traj_cfg.max_forward_velocity, self._traj_cfg.max_forward_acceleration,  self._traj_cfg.max_angular_velocity,  self._traj_cfg.max_angular_acceleration])

    def _setup_time_optimal_cost(self) -> None:
        self._J = self._dt_sym

    def solve(self, x_init:np.ndarray, x_e:np.ndarray) -> None:
        self._x_sol[0,:] = x_init[0] + (x_e[0]- x_init[0])*np.linspace(0, 1, self._path_cfg.n_hrzn+1, endpoint=True)
        self._x_sol[1,:] = (x_e[0]- x_init[0]) / self._track_spline.spline_length * self._traj_cfg.max_forward_velocity
        self._u_sol[0, :] = 0.0
        self._dt_sol = self._traj_cfg.max_forward_velocity / ((x_e[0]- x_init[0]) * self._track_spline.spline_length) / self._path_cfg.n_hrzn
        solution = self._ocp_solver(
            x0 = np.vstack(( self._x_sol.reshape((-1, 1), order='F'), self._u_sol.reshape((-1, 1), order='F'), self._dt_sol)),
            p = np.concatenate((x_init, x_e)),
            lbg = np.concatenate(self._lb_g),
            ubg = np.concatenate(self._ub_g),
        )
        primal_sol = solution["x"].full().flatten()
        self._x_sol = primal_sol[:(self._path_cfg.n_hrzn+1)*self._path_cfg.nx].reshape((self._path_cfg.nx, -1), order='F')
        self._u_sol = primal_sol[(self._path_cfg.n_hrzn+1)*self._path_cfg.nx:-1].reshape((self._path_cfg.nu, -1), order='F')
        self._dt_sol = primal_sol[-1]
        self._compute_whole_reference_trajectory()

    def plot_results(self) -> None:
        self._track_spline.plot_trajectories_over_time(s=self._x_sol[0, :], v_s=self._x_sol[1,:], a_s=self._u_sol[0,:], timestamps=self._spline_timestamps)

    def _compute_whole_reference_trajectory(self) -> None:
        for i_hrzn in range(0, self._path_cfg.n_hrzn+1):
            self._x_robot_ref[i_hrzn, 0] = self._track_spline.ca_fun_pos_x(self._x_sol[0, i_hrzn]).full()
            self._x_robot_ref[i_hrzn, 1] = self._track_spline.ca_fun_pos_y(self._x_sol[0, i_hrzn]).full()
            self._x_robot_ref[i_hrzn, 2] = self._track_spline.ca_fun_heading(self._x_sol[0, i_hrzn]).full()
            self._x_robot_ref[i_hrzn, 3] = self._track_spline.ca_fun_v_t(self._x_sol[0, i_hrzn], self._x_sol[1, i_hrzn]).full()
            self._x_robot_ref[i_hrzn, 4] = self._track_spline.ca_fun_omega_t(self._x_sol[0, i_hrzn], self._x_sol[1, i_hrzn]).full()
            if i_hrzn < self._path_cfg.n_hrzn:
                self._u_robot_ref[i_hrzn, 0] = self._track_spline.ca_fun_a_t(self._x_sol[0, i_hrzn], self._x_sol[1, i_hrzn], self._u_sol[0, i_hrzn]).full()
                self._u_robot_ref[i_hrzn, 1] = self._track_spline.ca_fun_alpha_t(self._x_sol[0, i_hrzn], self._x_sol[1, i_hrzn], self._u_sol[0, i_hrzn]).full()
        self._spline_timestamps = np.linspace(0, self._path_cfg.n_hrzn*self._dt_sol, self._path_cfg.n_hrzn+1, endpoint=True)
        self.f_x_ref = interp1d(self._spline_timestamps, self._x_robot_ref, axis=0)
        self.f_u_ref = interp1d(self._spline_timestamps[:-1], self._u_robot_ref, axis=0)

    def interpolate_reference_trajectory(self, robot_state:np.ndarray):
        dist2spline = np.linalg.norm(robot_state[:2] - self._x_robot_ref[:,:2], axis=1)
        idx = np.argmin(dist2spline)
        ref_traj_timestamps = idx*self._dt_sol + np.linspace(0, self._traj_cfg.n_hrzn*self._traj_cfg.delta_t, self._traj_cfg.n_hrzn+1, endpoint=True)
        x_ref_interp = self.f_x_ref(ref_traj_timestamps)
        u_ref_interp = self.f_u_ref(ref_traj_timestamps)
        return x_ref_interp, u_ref_interp

    @property
    def x_robot_ref(self) -> np.ndarray:
        return self._x_robot_ref

    @property
    def u_robot_ref(self) -> np.ndarray:
        return self._u_robot_ref

    @property
    def x_sol(self) -> np.ndarray:
        return self._x_sol

    @property
    def u_sol(self) -> np.ndarray:
        return self._u_sol

    @property
    def spline_timestamps(self) -> np.ndarray:
        return self._spline_timestamps

if __name__ == '__main__':
    control_points = np.array([[-2.0, 2.0],
                               [0.0, 0.0],
                               [2.0, 0.0],
                               [4.0, 0.0],
                               [6.0, 0.0],
                               [8.0, 2.0],
                               [8.0, 4.0],
                               [8.0, 6.0],
                               [8.0, 8.0]])
    spline_seg_length = np.array([0.0, np.pi, 2.0, 2.0, 2.0, np.pi, 2.0, 2.0, 2.0])
    track_spline = TrackSpline(control_points, spline_seg_length)
    cfg_path = PathTrackingParam()
    cfg_traj = MPCParam()
    path_tracking_solver = NominalPathTrackingSolver(track_spline, cfg_path=cfg_path, cfg_traj=cfg_traj)
    x_init = np.array([0.0, cfg_path.v_s_0])
    x_e    = np.array([1.0, cfg_path.v_s_e])
    path_tracking_solver.solve(x_init=x_init, x_e=x_e)
    path_tracking_solver.plot_results()
    plt.show()