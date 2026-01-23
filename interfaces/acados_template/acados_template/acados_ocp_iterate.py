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
from deprecated.sphinx import deprecated
import warnings


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
    
    def __getitem__(self, key) -> AcadosOcpFlattenedIterate:
        if isinstance(key, slice):
            ValueError("Slicing of batch iterates not supported.")
        return AcadosOcpFlattenedIterate(
            x=self.x[key],
            u=self.u[key],
            z=self.z[key],
            sl=self.sl[key],
            su=self.su[key],
            pi=self.pi[key],
            lam=self.lam[key],
        )

    def __setitem__(self, idx, value: AcadosOcpFlattenedIterate):
        self.x[idx] = value.x
        self.u[idx] = value.u
        self.z[idx] = value.z
        self.sl[idx] = value.sl
        self.su[idx] = value.su
        self.pi[idx] = value.pi
        self.lam[idx] = value.lam


@dataclass
class AcadosOcpIterate:
    """
    This class is used to store the primal-dual iterate of an optimal control problem.
    """

    x: List[np.ndarray]
    u: List[np.ndarray]
    z: List[np.ndarray]
    sl: List[np.ndarray]
    su: List[np.ndarray]
    pi: List[np.ndarray]
    lam: List[np.ndarray]

    def __init__(self, x=None, u=None, z=None, sl=None, su=None, pi=None, lam=None, x_traj=None, u_traj=None, z_traj=None, sl_traj=None, su_traj=None, pi_traj=None, lam_traj=None):

        self.x = x
        self.u = u
        self.z = z
        self.sl = sl
        self.su = su
        self.pi = pi
        self.lam = lam

        if x_traj is not None:
            self.x = x_traj
            warnings.warn("Parameter 'x_traj' is deprecated, use 'x' instead.", DeprecationWarning, stacklevel=2)
        if u_traj is not None:
            self.u = u_traj
            warnings.warn("Parameter 'u_traj' is deprecated, use 'u' instead.", DeprecationWarning, stacklevel=2)
        if z_traj is not None:
            self.z = z_traj
            warnings.warn("Parameter 'z_traj' is deprecated, use 'z' instead.", DeprecationWarning, stacklevel=2)
        if sl_traj is not None:
            self.sl = sl_traj
            warnings.warn("Parameter 'sl_traj' is deprecated, use 'sl' instead.", DeprecationWarning, stacklevel=2)
        if su_traj is not None:
            self.su = su_traj
            warnings.warn("Parameter 'su_traj' is deprecated, use 'su' instead.", DeprecationWarning, stacklevel=2)
        if pi_traj is not None:
            self.pi = pi_traj
            warnings.warn("Parameter 'pi_traj' is deprecated, use 'pi' instead.", DeprecationWarning, stacklevel=2)
        if lam_traj is not None:
            self.lam = lam_traj
            warnings.warn("Parameter 'lam_traj' is deprecated, use 'lam' instead.", DeprecationWarning, stacklevel=2)

        # TODO this is only required as long as the deprecated parameters are supported
        # remove and make all arguments required when deprecated arguments are removed

        # Validate that all required fields are provided
        all_provided = all(val is not None for val in self.__dict__.values())

        if not all_provided:
            raise ValueError("AcadosOcpIterate requires all of (x, u, z, sl, su, pi, lam) to be provided.")


    @property
    @deprecated(version="0.5.4", reason="Property 'x_traj' is deprecated, use 'x' instead.")
    def x_traj(self) -> List[np.ndarray]:
        return self.x

    @property
    @deprecated(version="0.5.4", reason="Property 'u_traj' is deprecated, use 'u' instead.")
    def u_traj(self) -> List[np.ndarray]:
        return self.u

    @property
    @deprecated(version="0.5.4", reason="Property 'z_traj' is deprecated, use 'z' instead.")
    def z_traj(self) -> List[np.ndarray]:
        return self.z

    @property
    @deprecated(version="0.5.4", reason="Property 'sl_traj' is deprecated, use 'sl' instead.")
    def sl_traj(self) -> List[np.ndarray]:
        return self.sl

    @property
    @deprecated(version="0.5.4", reason="Property 'su_traj' is deprecated, use 'su' instead.")
    def su_traj(self) -> List[np.ndarray]:
        return self.su

    @property
    @deprecated(version="0.5.4", reason="Property 'pi_traj' is deprecated, use 'pi' instead.")
    def pi_traj(self) -> List[np.ndarray]:
        return self.pi

    @property
    @deprecated(version="0.5.4", reason="Property 'lam_traj' is deprecated, use 'lam' instead.")
    def lam_traj(self) -> List[np.ndarray]:
        return self.lam

    def flatten(self) -> AcadosOcpFlattenedIterate:
        return AcadosOcpFlattenedIterate(
            x=np.concatenate(self.x),
            u=np.concatenate(self.u) if len(self.u) > 0 else np.array([]),
            z=np.concatenate(self.z) if len(self.z) > 0 else np.array([]),
            sl=np.concatenate(self.sl),
            su=np.concatenate(self.su),
            pi=np.concatenate(self.pi) if len(self.pi) > 0 else np.array([]),
            lam=np.concatenate(self.lam),
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
    __iterate_fields = ("x", "u", "z", "sl", "su", "pi", "lam")

    def as_array(self, field: str, ) -> np.ndarray:
        """
        Returns the iterate given by field as a numpy array of shape (nlp_iter+1, N_horizon (+ 1), n_field).
        This will fail if the dimension of value `field` is varying stagewise.
        """

        if field not in self.__iterate_fields:
            raise ValueError(f"Invalid field: got {field}, expected value in {self.__iterate_fields}")

        traj_ = [getattr(iterate, field) for iterate in self.iterate_list]

        try:
            traj = np.array(traj_, dtype=float)
        except ValueError:
            raise ValueError(f"Stage-wise dimensions are not the same for {field} trajectory.")

        return traj
