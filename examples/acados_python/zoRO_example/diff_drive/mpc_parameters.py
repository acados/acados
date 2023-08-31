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
    _nx: int=5
    _nu: int=2
    _nw: int=5
    _delta_t: float=0.1
    _n_hrzn: int=20

    # matrix of the cost function
    _Q: np.ndarray=np.zeros(0)
    _R: np.ndarray=np.zeros(0)
    _Q_e: np.ndarray=np.zeros(0)

    # constraints
    _num_state_cstr: int=2
    _min_forward_velocity: float=0.
    _max_forward_velocity: float=1.0
    _max_angular_velocity: float=1.0
    _min_forward_acceleration: float=-1.0
    _max_forward_acceleration: float=0.3
    _max_angular_acceleration: float=2.84
    _term_forward_velocity: float=0.05
    _term_angular_velocity: float=0.05

    # feedback matrix
    _fdbk_k: float=6.0
    _fdbk_K_mat: np.ndarray=np.zeros(0)

    # uncertainty / distrubance
    _unc_jac_G_mat: np.ndarray=np.zeros(0)
    _W_mat: np.ndarray=np.zeros(0)
    _P0_mat: np.ndarray=np.zeros(0)

    # obstacles
    _num_obs: int=3
    _obs_radius: np.ndarray=np.array([1.0, 0.7, 0.55])
    _obs_pos: np.ndarray=np.array([[0.0, 1.0], [3.0, 0.68], [7.0, 1.2]])

    # zoRO
    _backoff_eps: float=1e-8
    _zoRO_iter: int=2
    _use_custom_update: bool=True

    def __post_init__(self):
        self._Q: np.eye(self._nx)
        self._R: np.eye(self._nu) * 1e-1
        self._Q_e: np.eye(self._nx)

        self._fdbk_K_mat = np.array([[0., 0., 0., self._fdbk_k, 0.], \
            [0., 0., 0., 0., self._fdbk_k]])
        self._unc_jac_G_mat = np.eye(self._nx)
        self._W_mat = np.diag([2.0e-06, 2.0e-06, 4.0e-06, 1.5e-03, 7.0e-03])
        self._P0_mat = np.diag([2.0e-06, 2.0e-06, 4.0e-06, 1.5e-03, 7.0e-03])

    @property
    def nx(self)->int:
        return self._nx

    @property
    def nu(self)->int:
        return self._nu

    @property
    def nw(self)->int:
        return self._nw

    @property
    def delta_t(self)->float:
        return self._delta_t

    @property
    def n_hrzn(self)->int:
        return self._n_hrzn

    @property
    def Q(self)->np.ndarray:
        return np.eye(self._nx)

    @property
    def R(self)->np.ndarray:
        return np.eye(self._nu) * 1e-1

    @property
    def Q_e(self)->np.ndarray:
        return np.eye(self._nx)

    @property
    def num_state_cstr(self)->int:
        return self._num_state_cstr

    @property
    def min_forward_velocity(self)->float:
        return self._min_forward_velocity

    @property
    def max_forward_velocity(self)->float:
        return self._max_forward_velocity

    @property
    def max_angular_velocity(self)->float:
        return self._max_angular_velocity

    @property
    def min_forward_acceleration(self)->float:
        return self._min_forward_acceleration

    @property
    def max_forward_acceleration(self)->float:
        return self._max_forward_acceleration

    @property
    def max_angular_acceleration(self)->float:
        return self._max_angular_acceleration

    @property
    def term_forward_velocity(self)->float:
        return self._term_forward_velocity

    @property
    def term_angular_velocity(self)->float:
        return self._term_angular_velocity

    @property
    def fdbk_k(self)->float:
        return self._fdbk_k

    @property
    def fdbk_K_mat(self)->np.ndarray:
        return np.array([[0., 0., 0., self._fdbk_k, 0.], \
                         [0., 0., 0., 0., self._fdbk_k]])

    @property
    def unc_jac_G_mat(self)->np.ndarray:
        return np.eye(self._nw)

    @property
    def W_mat(self)->np.ndarray:
        return self._W_mat

    @property
    def P0_mat(self)->np.ndarray:
        return self._P0_mat

    @property
    def num_obs(self)->int:
        return self._num_obs

    @num_obs.setter
    def num_obs(self, n:int):
        self._num_obs = n

    @property
    def obs_radius(self)->np.ndarray:
        return self._obs_radius

    @obs_radius.setter
    def obs_radius(self, radius: np.ndarray):
        assert radius.size == self._num_obs
        self._obs_radius = radius

    @property
    def obs_pos(self)->np.ndarray:
        return self._obs_pos

    @obs_pos.setter
    def obs_pos(self, pos: np.ndarray):
        assert pos.size[0] == self._num_obs and  pos.size[1] == 2
        self._obs_pos = pos

    @property
    def backoff_eps(self)->float:
        return self._backoff_eps

    @property
    def zoRO_iter(self)->int:
        return self._zoRO_iter

    @zoRO_iter.setter
    def zoRO_iter(self, n: int):
        self._zoRO_iter = n

    @property
    def use_custom_update(self)->bool:
        return self._use_custom_update

    @use_custom_update.setter
    def use_custom_update(self, val: bool):
        self._use_custom_update = val

@dataclass
class PathTrackingParam:
    _nx: int=2
    _nu: int=1
    _nu_wT:  int=2
    _n_hrzn: int=500
    _v_s_0:  float=0.001
    _v_s_e: float=0.001


    @property
    def nx(self)->int:
        return self._nx

    @property
    def nu(self)->int:
        return self._nu

    @property
    def n_hrzn(self)->int:
        return self._n_hrzn

    @property
    def nu_wT(self)->int:
        return self._nu_wT

    @property
    def v_s_0(self) -> float:
        return self._v_s_0

    @property
    def v_s_e(self) -> float:
        return self._v_s_e