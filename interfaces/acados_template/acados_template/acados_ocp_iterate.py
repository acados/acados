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
class AcadosOcpFlattenedIterate:
    """
    This class is used to store the primal-dual iterate of an optimal control problem in a flattened form.
    This allows faster interactions with the solver.
    """
    x: np.ndarray
    u: np.ndarray
    z: np.ndarray
    sl: np.ndarray
    su: np.ndarray
    pi: np.ndarray
    lam: np.ndarray

    def allclose(self, other, rtol=1e-5, atol=1e-6) -> bool:
        if not isinstance(other, AcadosOcpFlattenedIterate):
            raise TypeError(f"Expected AcadosOcpFlattenedIterate, got {type(other)}")
        return (
            np.allclose(self.x, other.x, rtol=rtol, atol=atol) and
            np.allclose(self.u, other.u, rtol=rtol, atol=atol) and
            np.allclose(self.z, other.z, rtol=rtol, atol=atol) and
            np.allclose(self.sl, other.sl, rtol=rtol, atol=atol) and
            np.allclose(self.su, other.su, rtol=rtol, atol=atol) and
            np.allclose(self.pi, other.pi, rtol=rtol, atol=atol) and
            np.allclose(self.lam, other.lam, rtol=rtol, atol=atol)
        )

    def __add__(self, other):
        if not isinstance(other, AcadosOcpFlattenedIterate):
            raise TypeError(f"Expected AcadosOcpFlattenedIterate, got {type(other)}")
        return AcadosOcpFlattenedIterate(
            x=self.x + other.x,
            u=self.u + other.u,
            z=self.z + other.z,
            sl=self.sl + other.sl,
            su=self.su + other.su,
            pi=self.pi + other.pi,
            lam=self.lam + other.lam
        )

    def __sub__(self, other):
        if not isinstance(other, AcadosOcpFlattenedIterate):
            raise TypeError(f"Expected AcadosOcpFlattenedIterate, got {type(other)}")
        return AcadosOcpFlattenedIterate(
            x=self.x - other.x,
            u=self.u - other.u,
            z=self.z - other.z,
            sl=self.sl - other.sl,
            su=self.su - other.su,
            pi=self.pi - other.pi,
            lam=self.lam - other.lam
        )

    def inf_norm(self) -> float:
        """
        Returns the infinity norm of the iterate, which is the maximum absolute value of its elements.
        """
        return np.max(np.abs(np.concatenate((self.x, self.u, self.z, self.sl, self.su, self.pi, self.lam))))


@dataclass
class AcadosOcpFlattenedBatchIterate:
    """
    Similar to :py:class:`~acados_template.acados_ocp.AcadosOcpFlattenedIterate` but with an additional field N_batch.
    The fields x, u, z, sl, su, pi, lam are of shape (N_batch, nx_total), (N_batch, nu_total), etc.
    """
    x: np.ndarray # shape (N_batch, nx_total)
    u: np.ndarray # shape (N_batch, nu_total)
    z: np.ndarray
    sl: np.ndarray
    su: np.ndarray
    pi: np.ndarray
    lam: np.ndarray
    N_batch: int


@dataclass
class AcadosOcpIterate:
    """
    This class is used to store the primal-dual iterate of an optimal control problem.
    """

    x_traj: List[np.ndarray]
    u_traj: List[np.ndarray]
    z_traj: List[np.ndarray]
    sl_traj: List[np.ndarray]
    su_traj: List[np.ndarray]
    pi_traj: List[np.ndarray]
    lam_traj: List[np.ndarray]

    def flatten(self) -> AcadosOcpFlattenedIterate:
        return AcadosOcpFlattenedIterate(
            x=np.concatenate(self.x_traj),
            u=np.concatenate(self.u_traj) if len(self.u_traj) > 0 else np.array([]),
            z=np.concatenate(self.z_traj) if len(self.z_traj) > 0 else np.array([]),
            sl=np.concatenate(self.sl_traj),
            su=np.concatenate(self.su_traj),
            pi=np.concatenate(self.pi_traj) if len(self.pi_traj) > 0 else np.array([]),
            lam=np.concatenate(self.lam_traj),
        )

    def allclose(self, other, rtol=1e-5, atol=1e-8) -> bool:
        if not isinstance(other, AcadosOcpIterate):
            raise TypeError(f"Expected AcadosOcpIterate, got {type(other)}")
        s = self.flatten()
        o = other.flatten()
        return s.allclose(o, rtol=rtol, atol=atol)

@dataclass
class AcadosOcpIterates:
    """
    This class is used to store a list of :py:class:`~acados_template.acados_ocp.AcadosOcpIterate` objects.
    """

    iterate_list: List[AcadosOcpIterate]
    __iterate_fields = ["x", "u", "z", "sl", "su", "pi", "lam"]

    def as_array(self, field: str, ) -> np.ndarray:
        """
        Returns the iterate given by field as a numpy array of shape (nlp_iter+1, N_horizon (+ 1), n_field).
        This will fail if the dimension of value `field` is varying stagewise.
        """

        if field not in self.__iterate_fields:
            raise ValueError(f"Invalid field: got {field}, expected value in {self.__iterate_fields}")

        attr = f"{field}_traj"
        traj_ = [getattr(iterate, attr) for iterate in self.iterate_list]

        try:
            traj = np.array(traj_, dtype=float)
        except ValueError:
            raise ValueError(f"Stage-wise dimensions are not the same for {field} trajectory.")

        return traj
