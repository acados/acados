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

from dataclasses import dataclass, field
import numpy as np

@dataclass
class MPCParam():
    # dimensions
    nx: int=5
    nu: int=2
    nw: int=5
    delta_t: float=0.1
    n_hrzn: int=20

    # matrix of the cost function
    Q: np.ndarray=field(default_factory=lambda: np.zeros(0))
    R: np.ndarray=field(default_factory=lambda: np.zeros(0))
    Q_e: np.ndarray=field(default_factory=lambda: np.zeros(0))

    # constraints
    num_state_cstr: int=2
    min_forward_velocity: float=0.
    max_forward_velocity: float=1.0
    max_angular_velocity: float=1.0
    min_forward_acceleration: float=-1.0
    max_forward_acceleration: float=0.3
    max_angular_acceleration: float=2.84
    term_forward_velocity: float=0.05
    term_angular_velocity: float=0.05

    # feedback gain scalar parameter (full structured matrix defined below)
    fdbk_k: float=6.0
    feedback_optimization_mode: str = "CONSTANT_FEEDBACK"
    # -1: Pre-computed Feedback
    #  0: Feedback gain computed using riccati with constant cost matrices
    #  1: Feedback gain computed using riccati with sum of constant cost matrices and Hessian of tightened constraints weighted by 1/h**2

    # uncertainty / distrubance
    unc_jac_G_mat: np.ndarray=field(default_factory=lambda: np.zeros(0))
    W_mat: np.ndarray=field(default_factory=lambda: np.zeros(0))
    P0_mat: np.ndarray=field(default_factory=lambda: np.zeros(0))

    # obstacles
    num_obs: int=3
    _obs_radius: np.ndarray=field(default_factory=lambda: np.array([1.0, 0.7, 0.55]))
    _obs_pos: np.ndarray=field(default_factory=lambda: np.array([[0.0, 1.0], [3.0, 0.68], [7.0, 1.2]]))

    # zoRO
    backoff_eps: float=1e-8
    zoRO_iter: int=2
    use_custom_update: bool=True

    def __post_init__(self):
        self.Q = np.eye(self.nx)
        self.R = np.eye(self.nu) * 1e-1
        self.Q_e = np.eye(self.nx)
        self.unc_jac_G_mat = np.eye(self.nx)
        self.W_mat = np.diag([2.0e-06, 2.0e-06, 4.0e-06, 1.5e-03, 7.0e-03]) * 4.0
        self.P0_mat = np.diag([2.0e-06, 2.0e-06, 4.0e-06, 1.5e-03, 7.0e-03])

    @property
    def fdbk_K_mat(self)->np.ndarray:
        return np.array([[0., 0., 0., self.fdbk_k, 0.], \
                         [0., 0., 0., 0., self.fdbk_k]])

    @property
    def obs_radius(self)->np.ndarray:
        return self._obs_radius

    @obs_radius.setter
    def obs_radius(self, radius: np.ndarray):
        assert radius.size == self.num_obs
        self._obs_radius = radius

    @property
    def obs_pos(self)->np.ndarray:
        return self._obs_pos

    @obs_pos.setter
    def obs_pos(self, pos: np.ndarray):
        assert pos.size[0] == self.num_obs and  pos.size[1] == 2
        self._obs_pos = pos


@dataclass
class PathTrackingParam:
    nx: int=2
    nu: int=1
    nu_wT: int=2
    n_hrzn: int=500
    v_s_0: float=0.001
    v_s_e: float=0.001
