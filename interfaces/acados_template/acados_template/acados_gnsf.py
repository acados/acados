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
#


from acados_template import AcadosModel, AcadosOcpDims, AcadosSimDims
import numpy as np
import casadi as ca
from typing import Union


class AcadosGnsfDims():

    def __init__(self, nx1: int, nz1: int, nuhat: int, ny: int, nout: int):

        self.__nx1 = nx1
        self.__nz1 = nz1
        self.__nuhat = nuhat
        self.__ny = ny
        self.__nout = nout

    @property
    def nx1(self):
        """GNSF: Dimension nx1."""
        return self.__nx1

    @property
    def nz1(self):
        """GNSF: Dimension nz1."""
        return self.__nz1

    @property
    def nuhat(self):
        """GNSF: Dimension nuhat."""
        return self.__nuhat

    @property
    def ny(self):
        """GNSF: Dimension ny."""
        return self.__ny

    @property
    def nout(self):
        """GNSF: Dimension nout."""
        return self.__nout



class AcadosGnsfModel():

    def __init__(self,
                 y: Union[ca.SX, ca.MX],
                 uhat: Union[ca.SX, ca.MX],
                 phi: Union[ca.SX, ca.MX],
                 f_LO: Union[ca.SX, ca.MX],
                 A: np.ndarray,
                 B: np.ndarray,
                 C: np.ndarray,
                 E: np.ndarray,
                 L_x: np.ndarray,
                 L_xdot: np.ndarray,
                 L_u: np.ndarray,
                 L_z: np.ndarray,
                 A_LO: np.ndarray,
                 c: np.ndarray,
                 E_LO: np.ndarray,
                 B_LO: np.ndarray,
                 c_LO: np.ndarray,
                 ipiv_x: np.ndarray,
                 ipiv_z: np.ndarray,
                 purely_linear: bool,
                 nontrivial_f_LO: bool,
                 ):

        # symbolics and expressions
        self.__y = y
        self.__uhat = uhat
        self.__phi = phi
        self.__f_LO = f_LO

        # matrices and vectors
        self.__A = A
        self.__B = B
        self.__C = C
        self.__E = E
        self.__L_x = L_x
        self.__L_xdot = L_xdot
        self.__L_u = L_u
        self.__L_z = L_z
        self.__A_LO = A_LO
        self.__c = c
        self.__E_LO = E_LO
        self.__B_LO = B_LO
        self.__c_LO = c_LO

        self.__ipiv_x = ipiv_x
        self.__ipiv_z = ipiv_z

        self.__purely_linear = purely_linear
        self.__nontrivial_f_LO = nontrivial_f_LO

        self._dims = None

        self._detect_dims()
        self._make_consistent()

    @property
    def y(self):
        """GNSF: Symbolic expression for y in GNSF formulation."""
        return self.__y

    @property
    def uhat(self):
        """GNSF: Symbolic expression for uhat in GNSF formulation."""
        return self.__uhat

    @property
    def phi(self):
        """GNSF: Symbolic expression for phi in GNSF formulation."""
        return self.__phi

    @property
    def f_LO(self):
        """GNSF: Symbolic expression for f_LO (linear output function) in GNSF formulation."""
        return self.__f_LO

    @property
    def A(self):
        """GNSF: Matrix A in GNSF formulation."""
        return self.__A

    @property
    def B(self):
        """GNSF: Matrix B in GNSF formulation."""
        return self.__B

    @property
    def C(self):
        """GNSF: Matrix C in GNSF formulation."""
        return self.__C

    @property
    def E(self):
        """GNSF: Matrix E in GNSF formulation."""
        return self.__E

    @property
    def L_x(self):
        """GNSF: Matrix L_x in GNSF formulation."""
        return self.__L_x

    @property
    def L_xdot(self):
        """GNSF: Matrix L_xdot in GNSF formulation."""
        return self.__L_xdot

    @property
    def L_u(self):
        """GNSF: Matrix L_u in GNSF formulation."""
        return self.__L_u

    @property
    def L_z(self):
        """GNSF: Matrix L_z in GNSF formulation."""
        return self.__L_z

    @property
    def A_LO(self):
        """GNSF: Matrix A_LO (linear output) in GNSF formulation."""
        return self.__A_LO

    @property
    def c(self):
        """GNSF: Vector c in GNSF formulation."""
        return self.__c

    @property
    def E_LO(self):
        """GNSF: Matrix E_LO (linear output) in GNSF formulation."""
        return self.__E_LO

    @property
    def B_LO(self):
        """GNSF: Matrix B_LO (linear output) in GNSF formulation."""
        return self.__B_LO

    @property
    def nontrivial_f_LO(self):
        """GNSF: Flag indicating whether GNSF structure has nontrivial f_LO."""
        return self.__nontrivial_f_LO

    @property
    def purely_linear(self):
        """GNSF: Flag indicating whether GNSF structure is purely linear."""
        return self.__purely_linear

    @property
    def ipiv_x(self):
        """GNSF: Pivot indices for x in GNSF formulation."""
        return self.__ipiv_x

    @property
    def ipiv_z(self):
        """GNSF: Pivot indices for z in GNSF formulation."""
        return self.__ipiv_z

    @property
    def c_LO(self):
        """GNSF: Vector c_LO (linear output) in GNSF formulation."""
        return self.__c_LO


    def _detect_dims(self,):
        """
        Detect dimensions of GNSF model.
        """

        nx1 = self.A.shape[1]
        nz1 = self.A.shape[0] - nx1
        nuhat = self.L_u.shape[0]
        ny = self.L_x.shape[0]
        nout = self.C.shape[1]

        self.__dims = AcadosGnsfDims(nx1, nz1, nuhat, ny, nout)


    def _make_consistent(self,):
        """
        Make sure the GNSF model is consistent.
        """

        # sanity checks

        pass


    @classmethod
    def detect(cls, model: AcadosModel, dims: Union[AcadosOcpDims, AcadosSimDims]):
        """
        Detect if given AcadosModel can be represented in GNSF form.

        Args:
            model: AcadosModel to be checked
        """

        gnsf_model = cls()


        return gnsf_model


    def to_dict(self):
        """
        Convert AcadosGnsfModel to dictionary.

        Returns:
            dict: dictionary representation of AcadosGnsfModel
        """

        model_dict = {}
        return model_dict


    @classmethod
    def from_dict(cls, model_dict: dict):
        """
        Create AcadosGnsfModel from dictionary.

        Args:
            model_dict: dictionary representation of AcadosGnsfModel
        """

        gnsf_model = cls()
        return gnsf_model