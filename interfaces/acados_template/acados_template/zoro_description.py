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
class ZoroDescription:
    """
    Zero-Order Robust Optimization (zoRO) scheme.

    The uncertainty propagation is performed by:
    $$P_{k+1} = (A_k + B_kK)P_k(A_k + B_kK)^\top + GWG^\top$$

    For advanced users.
    """
    backoff_scaling_gamma: float = 1.0
    """backoff scaling factor, for stochastic MPC"""
    fdbk_K_mat: np.ndarray = None
    """constant feedback gain matrix K"""
    unc_jac_G_mat: np.ndarray = None    # default: an identity matrix
    """matrix G, describes how noise enters the dynamics"""
    P0_mat: np.ndarray = None
    """initial uncertainty matrix $\bar{P}_0$"""
    W_mat: np.ndarray = None
    """matrix W, covariance of noise in stochastic setting, defines uncertainty ellipsoids in robust setting"""
    idx_lbx_t: list = field(default_factory=list)
    """Indices of constraints to be tightened within the lower bounds on x for intermediate shooting nodes 1,...,N-1"""
    idx_ubx_t: list = field(default_factory=list)
    """Indices of constraints to be tightened within the upper bounds on x for intermediate shooting nodes 1,...,N-1"""
    idx_lbx_e_t: list = field(default_factory=list)
    """Indices of constraints to be tightened within the lower bounds on x for terminal shooting node"""
    idx_ubx_e_t: list = field(default_factory=list)
    """Indices of constraints to be tightened within the upper bounds on x for terminal shooting node"""
    idx_lbu_t: list = field(default_factory=list)
    """Indices of constraints to be tightened within the lower bounds on u for intermediate shooting nodes 1,...,N-1"""
    idx_ubu_t: list = field(default_factory=list)
    """Indices of constraints to be tightened within the upper bounds on u for intermediate shooting nodes 1,...,N-1"""
    idx_lg_t: list = field(default_factory=list)
    idx_ug_t: list = field(default_factory=list)
    idx_lg_e_t: list = field(default_factory=list)
    idx_ug_e_t: list = field(default_factory=list)
    idx_lh_t: list = field(default_factory=list)
    idx_uh_t: list = field(default_factory=list)
    idx_lh_e_t: list = field(default_factory=list)
    idx_uh_e_t: list = field(default_factory=list)

def process_zoro_description(zoro_description: ZoroDescription):
    zoro_description.nw, _ = zoro_description.W_mat.shape
    if zoro_description.unc_jac_G_mat is None:
        zoro_description.unc_jac_G_mat = np.eye(zoro_description.nw)
    zoro_description.nlbx_t = len(zoro_description.idx_lbx_t)
    zoro_description.nubx_t = len(zoro_description.idx_ubx_t)
    zoro_description.nlbx_e_t = len(zoro_description.idx_lbx_e_t)
    zoro_description.nubx_e_t = len(zoro_description.idx_ubx_e_t)
    zoro_description.nlbu_t = len(zoro_description.idx_lbu_t)
    zoro_description.nubu_t = len(zoro_description.idx_ubu_t)
    zoro_description.nlg_t = len(zoro_description.idx_lg_t)
    zoro_description.nug_t = len(zoro_description.idx_ug_t)
    zoro_description.nlg_e_t = len(zoro_description.idx_lg_e_t)
    zoro_description.nug_e_t = len(zoro_description.idx_ug_e_t)
    zoro_description.nlh_t = len(zoro_description.idx_lh_t)
    zoro_description.nuh_t = len(zoro_description.idx_uh_t)
    zoro_description.nlh_e_t = len(zoro_description.idx_lh_e_t)
    zoro_description.nuh_e_t = len(zoro_description.idx_uh_e_t)
    return zoro_description
