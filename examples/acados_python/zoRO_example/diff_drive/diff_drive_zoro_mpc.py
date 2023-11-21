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
from scipy.linalg import block_diag
from time import process_time

from acados_template import AcadosOcp, AcadosOcpSolver, ZoroDescription

from diff_drive_model import export_diff_drive_model, RobotState
from mpc_parameters import MPCParam

class ZoroMPCSolver:
    def __init__(self, cfg: MPCParam) -> None:
        # import model
        self.cfg = cfg
        self.model = export_diff_drive_model(cfg)

        # create ocp object to formulate the OCP
        self.ocp = AcadosOcp()
        self.ocp.model = self.model

        # set prediction horizon
        self.ocp.dims.N = cfg.n_hrzn
        self.ocp.solver_options.tf = cfg.n_hrzn * cfg.delta_t

        # set stage cost = ||u(k)-u_r(k)||_R^2 + ||x(k)-x_r(k)||_Q^2
        self.ocp.cost.cost_type = 'LINEAR_LS'
        cost_QR = block_diag(cfg.Q, cfg.R)
        self.ocp.cost.W = cost_QR
        self.ocp.cost.Vx = np.zeros((cfg.nu+cfg.nx, cfg.nx))
        self.ocp.cost.Vx[:cfg.nx, :cfg.nx] = np.eye(cfg.nx)
        self.ocp.cost.Vu = np.zeros((cfg.nu+cfg.nx, cfg.nu))
        self.ocp.cost.Vu[cfg.nx:, :] = np.eye(cfg.nu)
        self.ocp.cost.yref  = np.zeros((cfg.nu+cfg.nx, ))

        # set terminal cost = ||x_e - x_ref_e||_{Q_e}^2
        self.ocp.cost.cost_type_e = 'LINEAR_LS'
        self.ocp.cost.W_e = cfg.Q_e
        self.ocp.cost.Vx_e = np.eye(cfg.nx)
        self.ocp.cost.yref_e = np.zeros((cfg.nx, ))

        # input constraints
        self.lba = cfg.min_forward_acceleration
        self.uba = cfg.max_forward_acceleration
        self.ub_ang_a = cfg.max_angular_acceleration
        self.ocp.constraints.lbu = np.array([self.lba, -self.ub_ang_a])
        self.ocp.constraints.ubu = np.array([self.uba,  self.ub_ang_a])
        self.ocp.constraints.idxbu = np.array(range(cfg.nu))

        # state constraints
        self.lbv = cfg.min_forward_velocity
        self.ubv = cfg.max_forward_velocity
        self.ubw = cfg.max_angular_velocity
        self.ocp.constraints.idxbx = np.array([RobotState.VEL.value, RobotState.OMEGA.value])
        self.ocp.constraints.lbx = np.array([self.lbv, -self.ubw])
        self.ocp.constraints.ubx = np.array([self.ubv,  self.ubw])
        # initial state
        self.ocp.constraints.x0 = np.zeros((cfg.nx, ))
        # terminal state
        v_e = cfg.term_forward_velocity
        omega_e = cfg.term_angular_velocity
        self.ocp.constraints.idxbx_e = np.array([RobotState.VEL.value, RobotState.OMEGA.value])
        self.ocp.constraints.lbx_e = np.array([-v_e, -omega_e])
        self.ocp.constraints.ubx_e = np.array([v_e, omega_e])

        # Set collision avoidance constraints
        num_obs = cfg.num_obs
        self.ocp.dims.np = num_obs*2
        self.ocp.parameter_values = np.zeros((num_obs*2, ))
        self.ocp.constraints.lh = cfg.obs_radius
        self.ocp.constraints.uh = 1e3 * np.ones((num_obs, ))
        self.ocp.constraints.lh_e = cfg.obs_radius
        self.ocp.constraints.uh_e = 1e3 * np.ones((num_obs, ))

        # custom update: disturbance propagation
        self.ocp.solver_options.custom_update_filename = 'custom_update_function.c'
        self.ocp.solver_options.custom_update_header_filename = 'custom_update_function.h'

        self.ocp.solver_options.custom_update_copy = False
        self.ocp.solver_options.custom_templates = [
            ('custom_update_function_zoro_template.in.c', 'custom_update_function.c'),
            ('custom_update_function_zoro_template.in.h', 'custom_update_function.h'),
        ]

        # solver options
        self.ocp.solver_options.qp_solver = 'FULL_CONDENSING_DAQP' #'FULL_CONDENSING_QPOASES', 'FULL_CONDENSING_DAQP'
        self.ocp.solver_options.hessian_approx = 'GAUSS_NEWTON'
        self.ocp.solver_options.nlp_solver_type = 'SQP_RTI' # SQP, SQP_RTI
        self.ocp.solver_options.qp_tol = 1e-4

        # integrator options
        self.ocp.solver_options.integrator_type = 'IRK'
        self.ocp.solver_options.sim_method_num_stages = 2
        self.ocp.solver_options.sim_method_jac_reuse = 1

        # zoro stuff
        zoro_description = ZoroDescription()
        zoro_description.backoff_scaling_gamma = 3.0
        # uncertainty propagation: P_{k+1} = (A_k+B_kK) @ P_k @ (A_k+B_kK)^T + G @ W @ G^T
        # G.shape = (nx, nw), W.shape = (nw, nw)
        zoro_description.fdbk_K_mat = cfg.fdbk_K_mat
        zoro_description.unc_jac_G_mat = cfg.unc_jac_G_mat # G above
        zoro_description.P0_mat = cfg.P0_mat
        zoro_description.W_mat = cfg.W_mat
        zoro_description.idx_lbx_t = [1]
        zoro_description.idx_ubx_t = [0, 1]
        zoro_description.idx_lh_t = list(range(0, cfg.num_obs))
        zoro_description.idx_uh_t = []
        zoro_description.idx_lh_e_t = list(range(0, cfg.num_obs))
        zoro_description.idx_uh_e_t = []

        ## dummy linear constraints for testing
        # self.ocp.constraints.C = np.array([[1., 0., 0., 0., 0.0], [0., 1., 0., 0., 0.]])
        # self.ocp.constraints.D = np.array([[-0., 0.], [-0., 0.]])
        # self.ocp.constraints.C_e = np.array([[1., 0., 0., 0., 0.0], [0., 1., 0., 0., 0.]])
        # self.ocp.constraints.lg = np.array([-1e3, -1e3])
        # self.ocp.constraints.ug = np.array([1e3, 1e3])
        # self.ocp.constraints.lg_e = np.array([-1e3, -1e3])
        # self.ocp.constraints.ug_e = np.array([1e3, 1e3])
        # zoro_description.idx_lg_t = [0,1]
        # zoro_description.idx_ug_t = [0,1]
        # zoro_description.idx_lg_e_t = [0,1]
        # zoro_description.idx_ug_e_t = [0,1]

        self.ocp.zoro_description = zoro_description

        self.acados_ocp_solver = AcadosOcpSolver(self.ocp, json_file = 'acados_ocp_' + self.model.name + '.json')
        # AcadosOcpSolver.generate(self.ocp, json_file='acados_ocp_' + self.model.name + '.json')
        # AcadosOcpSolver.build(self.ocp.code_export_directory, with_cython=True)
        # self.acados_ocp_solver = AcadosOcpSolver.create_cython_solver('acados_ocp_' + self.model.name + '.json')

        self.initialized = False

        # timers
        self.rti_phase1_t = 0.
        self.rti_phase2_t = 0.
        self.propagation_t = 0.
        self.acados_integrator_time = 0.
        self.acados_qp_time = 0.

    def solve(self, x_current, y_ref, obs_position, obs_radius, p0_mat=None):
        """
        x_current: np.ndarray (nx,)
        y_ref: np.ndarray, (n_hrzn+1, nx + nu)
        obs_position: np.ndarray, (num_obs * 2, ), [x0, y0, x1, y1, ...]
        obs_radius: np.ndarray, (num_obs, ), [r0, r1, ...]
        """
        if p0_mat is None:
            p0_mat = np.zeros((self.cfg.nx, self.cfg.nx))

        self.rti_phase1_t = 0.
        self.rti_phase2_t = 0.
        self.propagation_t = 0.
        self.acados_integrator_time = 0.
        self.acados_qp_time = 0.

        if not self.initialized:
            # initialize solver
            for i in range(0, self.cfg.n_hrzn):
                self.acados_ocp_solver.set(i, 'x', np.array(x_current))
                self.acados_ocp_solver.set(i, 'u', np.zeros((self.cfg.nu,)))
            self.acados_ocp_solver.set(self.cfg.n_hrzn, 'x', np.array(x_current))
            self.x_temp_sol = np.tile(x_current, (self.cfg.n_hrzn+1, 1))
            self.u_temp_sol = np.zeros((self.cfg.n_hrzn, self.cfg.nu))

            # Initialize P matrix
            self.P_mats = np.zeros((self.cfg.n_hrzn+1, self.cfg.nx, self.cfg.nx))
            self.P_mats[0,:,:] = p0_mat.copy()
            for i in range(1, self.cfg.n_hrzn+1):
                self.P_mats[i,:,:] = self.cfg.W_mat
            self.initialized = True

        # set the current state
        self.acados_ocp_solver.set(0, "lbx", x_current)
        self.acados_ocp_solver.set(0, "ubx", x_current)
        # set reference trajectory
        for i_stage in range(self.cfg.n_hrzn):
            self.acados_ocp_solver.cost_set(i_stage, 'yref', y_ref[i_stage,:])
        self.acados_ocp_solver.cost_set(self.cfg.n_hrzn, 'yref', y_ref[self.cfg.n_hrzn,:self.cfg.nx])

        for i_stage in range(self.cfg.n_hrzn+1):
            self.acados_ocp_solver.set(i_stage,"p", obs_position)

        for i_sqp in range(self.cfg.zoRO_iter):
            # preparation rti_phase
            self.acados_ocp_solver.options_set('rti_phase', 1)
            status = self.acados_ocp_solver.solve()
            self.rti_phase1_t += self.acados_ocp_solver.get_stats("time_tot")
            self.acados_integrator_time += self.acados_ocp_solver.get_stats("time_sim")

            if self.cfg.use_custom_update:
                t_start = process_time()
                self.acados_ocp_solver.custom_update([self.cfg.P0_mat.flatten()])
                self.propagation_t += process_time() - t_start
            else:
                t_start = process_time()
                self.propagate_and_update(obs_position=obs_position, obs_radius=obs_radius, p0_mat=self.cfg.P0_mat)
                self.propagation_t += process_time() - t_start

            # feedback rti_phase
            self.acados_ocp_solver.options_set('rti_phase', 2)
            status = self.acados_ocp_solver.solve()
            self.acados_qp_time += self.acados_ocp_solver.get_stats("time_qp")
            self.rti_phase2_t += self.acados_ocp_solver.get_stats("time_tot")

            # Get solution
            for i_stage in range(self.cfg.n_hrzn):
                self.x_temp_sol[i_stage,:] = self.acados_ocp_solver.get(i_stage, "x")
                self.u_temp_sol[i_stage,:] = self.acados_ocp_solver.get(i_stage, "u")
            self.x_temp_sol[self.cfg.n_hrzn,:] = self.acados_ocp_solver.get(self.cfg.n_hrzn, "x")

            residuals = self.acados_ocp_solver.get_residuals()

        u_opt = self.acados_ocp_solver.get(0, "u")

        return u_opt, status


    def check_inequality_bound_feasibility(self) -> bool:
        temp_lam = self.acados_ocp_solver.get(self.cfg.n_hrzn, "lam")
        collision_cstr_active = not (np.isclose(\
            temp_lam[self.cfg.num_state_cstr:self.cfg.num_state_cstr+self.cfg.num_obs], \
                b=0., atol=1e-7).all())
        cstr_index = list(range(self.cfg.nu+self.cfg.num_state_cstr, \
            self.cfg.nu+self.cfg.num_state_cstr+self.cfg.num_obs))
        i_stage = 1
        while (not collision_cstr_active) and (i_stage < self.cfg.n_hrzn):
            temp_lam = self.acados_ocp_solver.get(i_stage, "lam")
            collision_cstr_active =  not (np.isclose(temp_lam[cstr_index], b=0., atol=1e-7).all())
            i_stage += 1
        return collision_cstr_active


    def get_state_input_trajectories(self):
        return self.x_temp_sol[:self.cfg.n_hrzn+1,:]


    def get_coll_cstr_status(self) -> bool:
        temp_lam = self.acados_ocp_solver.get(self.cfg.n_hrzn, "lam")
        collision_cstr_active = not (np.isclose(\
            temp_lam[self.cfg.num_state_cstr:self.cfg.num_state_cstr+self.cfg.num_obs], \
                b=0., atol=1e-7).all())
        cstr_index = list(range(self.cfg.nu+self.cfg.num_state_cstr, \
            self.cfg.nu+self.cfg.num_state_cstr+self.cfg.num_obs))
        i_stage = 1
        while (not collision_cstr_active) and (i_stage < self.cfg.n_hrzn):
            temp_lam = self.acados_ocp_solver.get(i_stage, "lam")
            collision_cstr_active =  not (np.isclose(temp_lam[cstr_index], b=0., atol=1e-7).all())
            i_stage += 1
        return collision_cstr_active


    def propagate_and_update(self, obs_position, obs_radius, p0_mat):
        # debug_list = []
        lbx_tightened = np.zeros((self.cfg.num_state_cstr, ))
        ubx_tightened = np.zeros((self.cfg.num_state_cstr, ))
        dist_guess = np.zeros((self.cfg.num_obs, 2))
        Pj_diag = np.zeros((self.cfg.nx,))
        temp_P_mat = self.P_mats[0,:,:] = p0_mat.copy()

        i_mpc_stage = 0
        i_diff_ctrl_stage = 0

        backoff_scaling_gamma = self.ocp.zoro_description.backoff_scaling_gamma
        while i_mpc_stage < self.cfg.n_hrzn:
            # get the A matrix
            temp_A = self.acados_ocp_solver.get_from_qp_in(i_mpc_stage, "A")
            temp_B = self.acados_ocp_solver.get_from_qp_in(i_mpc_stage, "B")
            temp_AK = temp_A - temp_B @ self.cfg.fdbk_K_mat
            temp_P_mat = temp_AK @ temp_P_mat @ temp_AK.T + self.cfg.W_mat
            i_mpc_stage += 1
            self.P_mats[i_mpc_stage,:,:] = temp_P_mat.copy()

            """ Compute backoff using P, set bounds with backoff
            v - ubv + sqrt(P(3, 3)) <= 0,
            w - ubw + sqrt(P(4, 4)) <= 0,
            a - uba + k1*sqrt(P(3, 3)) <= 0,
            alpha - ubalpha + k2*sqrt(P(4, 4)) <= 0,
            """
            Pj = self.P_mats[i_mpc_stage,:,:]
            Pj_diag = np.diag(Pj)
            # dont tighten x constraints at last node.
            if i_mpc_stage < self.cfg.n_hrzn:
                # v, cstr_idx[0]
                lbx_tightened[0] = self.lbv
                ubx_tightened[0] = self.ubv - backoff_scaling_gamma * np.sqrt(Pj_diag[3] + self.cfg.backoff_eps)
                # w, cstr_idx[1]
                lbx_tightened[1] = -self.ubw + backoff_scaling_gamma * np.sqrt(Pj_diag[4] + self.cfg.backoff_eps)
                ubx_tightened[1] = self.ubw - backoff_scaling_gamma * np.sqrt(Pj_diag[4] + self.cfg.backoff_eps)
                self.acados_ocp_solver.constraints_set(i_mpc_stage, "lbx", lbx_tightened)
                self.acados_ocp_solver.constraints_set(i_mpc_stage, "ubx", ubx_tightened)

            # obstacle constraint not enforced at node 0
            if i_mpc_stage > 0:
                x_guess_i = self.x_temp_sol[i_mpc_stage, 0:2]
                for i_obs in range(self.cfg.num_obs):
                    dist_guess[i_obs,:] = x_guess_i - obs_position[2*i_obs:2*(i_obs+1)]
                dist_norm_mat = dist_guess @ Pj[:2,:2] @ dist_guess.T
                dist_norm = np.diag(dist_norm_mat)
                sqr_dist_val = dist_guess[:,0]**2 + dist_guess[:,1]**2
                h_backoff = backoff_scaling_gamma * np.sqrt(dist_norm / sqr_dist_val + self.cfg.backoff_eps)
                self.acados_ocp_solver.constraints_set(i_mpc_stage, "lh", obs_radius + h_backoff)
        return
