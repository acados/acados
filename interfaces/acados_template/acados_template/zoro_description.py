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
from .acados_dims import AcadosOcpDims


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
    # Inputs:
    input_P0_diag: bool = False
    """Determines if diag(P0) is an input to the custom update function"""
    input_P0: bool = True
    """Determines if P0 is an input to the custom update function, specified in column-major format"""

    input_W_diag: bool = False
    """Determines if diag(W) is an input to the custom update function"""
    input_W_add_diag: bool = False
    """
    Determines if the concatenation of diag(W_{add}^k) is an input to the custom update function

    In case this is used W_k = W + W_{add}^k.
    """

    # Outputs:
    output_P_matrices: bool = False
    """Determines if the matrices P_k are outputs of the custom update function"""

    data_size: int = 0
    """size of data vector when calling custom update, computed automatically"""


    def make_consistent(self, dims: AcadosOcpDims) -> None:
        self.nw, _ = self.W_mat.shape
        if self.unc_jac_G_mat is None:
            self.unc_jac_G_mat = np.eye(self.nw)
        self.nlbx_t = len(self.idx_lbx_t)
        self.nubx_t = len(self.idx_ubx_t)
        self.nlbx_e_t = len(self.idx_lbx_e_t)
        self.nubx_e_t = len(self.idx_ubx_e_t)
        self.nlbu_t = len(self.idx_lbu_t)
        self.nubu_t = len(self.idx_ubu_t)
        self.nlg_t = len(self.idx_lg_t)
        self.nug_t = len(self.idx_ug_t)
        self.nlg_e_t = len(self.idx_lg_e_t)
        self.nug_e_t = len(self.idx_ug_e_t)
        self.nlh_t = len(self.idx_lh_t)
        self.nuh_t = len(self.idx_uh_t)
        self.nlh_e_t = len(self.idx_lh_e_t)
        self.nuh_e_t = len(self.idx_uh_e_t)

        if self.input_P0_diag and self.input_P0:
            raise Exception("Only one of input_P0_diag and input_P0 can be True")

        # Print input note:
        print(f"\nThe data of the generated custom update function consists of the concatenation of:")
        i_component = 1
        data_size = 0
        if self.input_P0_diag:
            size_i = dims.nx
            print(f"{i_component}) input: diag(P0), size: [nx] = {size_i}")
            i_component += 1
            data_size += size_i
        if self.input_P0:
            size_i = dims.nx ** 2
            print(f"{i_component}) input: P0; full matrix in column-major format, size: [nx*nx] = {size_i}")
            i_component += 1
            data_size += size_i
        if self.input_W_diag:
            size_i = self.nw
            print(f"{i_component}) input: diag(W), size: [nw] = {size_i}")
            i_component += 1
            data_size += size_i
        if self.input_W_add_diag:
            size_i = dims.N * self.nw
            print(f"{i_component}) input: concatenation of diag(W_gp^k) for i=0,...,N-1, size: [N * nw] = {size_i}")
            i_component += 1
            data_size += size_i
        if self.output_P_matrices:
            size_i = dims.nx * dims.nx * (dims.N+1)
            print(f"{i_component}) output: concatenation of colmaj(P^k) for i=0,...,N, size: [nx*nx*(N+1)] = {size_i}")
            i_component += 1
            data_size += size_i

        self.data_size = data_size
        print("\n")
