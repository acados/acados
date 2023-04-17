import sys
import os
local_path = os.path.dirname(os.path.abspath(__file__))
mpc_source_dir = os.path.join(local_path, '..')
sys.path.append(mpc_source_dir)

import numpy as np
import casadi
from acados_template import AcadosOcp, AcadosOcpSolver, AcadosModel
import matplotlib.pyplot as plt

from track_spline import TrackSpline
from mpc_parameters import PathTrackingParam, MPCParam

class NominalPathTrackingSolver:
    def __init__(self, track_spline:TrackSpline, path_cfg: PathTrackingParam, traj_cfg: MPCParam) -> None:
        self._track_spline = track_spline
        self._path_cfg = path_cfg
        self._traj_cfg = traj_cfg
        self._dt_sol = 0.1
        self._x_sol  = np.zeros((path_cfg.nx, path_cfg.n_hrzn+1))
        self._u_sol  = np.zeros((path_cfg.nu, path_cfg.n_hrzn  ))
        self._x_ref  = np.zeros((traj_cfg.n_hrzn+1, traj_cfg.nx))
        self._spline_timestamps   = np.full((path_cfg.n_hrzn+1, path_cfg.nx ), np.nan)
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
        for i_hrzn in range(1, self._path_cfg.n_hrzn):
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

    def plot_results(self) -> None:
        self._spline_timestamps = np.linspace(0, self._path_cfg.n_hrzn*self._dt_sol, self._path_cfg.n_hrzn+1, endpoint=True)
        self._track_spline.plot_trajectories_over_time(s=self._x_sol[0, :], v_s=self._x_sol[1,:], a_s=self._u_sol[0,:], timestamps=self._spline_timestamps)
        plt.show()

    def compute_reference_trajectory(self, robot_state:np.ndarray) -> None:
        dist2spline = np.linalg.norm(robot_state[:2] - self._x_sol.T, axis=1)
        idx = np.argmin(dist2spline)
        ref_traj_timestamps = idx*self._dt_sol + np.linspace(0, self._traj_cfg.n_hrzn*self._traj_cfg.delta_t, self._path_cfg.n_hrzn+1, endpoint=True)
        s_ref = np.interp(ref_traj_timestamps, self._spline_timestamps, self._x_sol[0, :])
        v_ref = np.interp(ref_traj_timestamps, self._spline_timestamps, self._x_sol[1, :])
        self._x_ref[0, :] = robot_state.copy()
        for i_hrzn in range(1, self._traj_cfg.n_hrzn+1):
            self._x_ref[i_hrzn, 0] = self._track_spline.ca_fun_pos_x(s_ref[i_hrzn])
            self._x_ref[i_hrzn, 1] = self._track_spline.ca_fun_pos_y(s_ref[i_hrzn])
            self._x_ref[i_hrzn, 2] = self._track_spline.ca_fun_heading(s_ref[i_hrzn])
            self._x_ref[i_hrzn, 3] = self._track_spline.ca_fun_v_t(s_ref[i_hrzn], v_ref[i_hrzn])
            self._x_ref[i_hrzn, 4] = self._track_spline.ca_fun_omega_t(s_ref[i_hrzn], v_ref[i_hrzn])

    @property
    def x_ref(self) -> np.ndarray:
        return self._x_ref

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
    path_cfg = PathTrackingParam()
    traj_cfg = MPCParam()
    path_tracking_solver = NominalPathTrackingSolver(track_spline, path_cfg=path_cfg, traj_cfg=traj_cfg)
    x_init = np.array([0.0, path_cfg._v_s_0])
    x_e    = np.array([1.0, path_cfg._v_s_e])
    path_tracking_solver.solve(x_init=x_init, x_e=x_e)
    path_tracking_solver.plot_results()