# -*- coding: future_fstrings -*-
#
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


from dataclasses import dataclass
from typing import List
import numpy as np

@dataclass
class AcadosOcpIterate:

    x_traj: List[np.ndarray]
    u_traj: List[np.ndarray]
    z_traj: List[np.ndarray]
    sl_traj: List[np.ndarray]
    su_traj: List[np.ndarray]
    pi_traj: List[np.ndarray]
    lam_traj: List[np.ndarray]

@dataclass
class AcadosOcpIterates:

    iterate_list: List[AcadosOcpIterate]
    __iterate_fields = ["x", "u", "z", "sl", "su", "pi", "lam"]

    def as_array(self, field: str, ) -> np.ndarray:
        """
        Returns the iterate given by field as a numpy array of shape (nlp_iter+1, N_horizon (+ 1), n_field).
        This will fail if the dimension of value `field` is varying stagewise.
        """

        if field not in self.__iterate_fields:
            raise Exception(f"Invalid field: got {field}, expected value in {self.__iterate_fields}")

        attr = f"{field}_traj"
        traj_ = [getattr(iterate, attr) for iterate in self.iterate_list]

        try:
            traj = np.array(traj_, dtype=float)
        except ValueError:
            raise Exception(f"Stage-wise dimensions are not the same for {field} trajectory.")

        return traj
