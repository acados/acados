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
#

from typing import Optional
import numpy as np
from copy import deepcopy
from scipy.linalg import block_diag

import casadi as ca
import os

from .acados_dims import MultiphaseOcpDims
from .acados_model import AcadosModel
from .acados_ocp_cost import AcadosOcpCost
from .acados_ocp_constraints import AcadosOcpConstraints
from .acados_dims import AcadosOcpDims
from .acados_ocp_options import AcadosOcpOptions

from .acados_ocp import AcadosOcp

from .utils import (get_acados_path, format_class_dict,
                    get_shared_lib_ext, is_column, is_empty, casadi_length,)
from .penalty_utils import symmetric_huber_penalty



class AcadosMultiphaseOcp:
    """
    Class containing the description of an optimal control problem with multiple phases.

    This object can be used to create an :py:class:`acados_template.acados_ocp_solver.AcadosOcpSolver`.
    """
    def __init__(self, n_phases: int, N_list: list):

        if n_phases != len(N_list):
            raise Exception('Number of phases does not match the length of N_list.')
        self.n_phases = n_phases
        self.N_list = N_list

        self.name = 'multiphase_ocp'
        self.dims = MultiphaseOcpDims()
        self.model = n_phases * [AcadosModel()]
        """Model definitions, type :py:class:`acados_template.acados_model.AcadosModel`"""
        self.cost = n_phases * [AcadosOcpCost()]
        """Cost definitions, type :py:class:`acados_template.acados_ocp.AcadosOcpCost`"""
        self.constraints = n_phases * [AcadosOcpConstraints()]
        """Constraints definitions, type :py:class:`acados_template.acados_ocp.AcadosOcpConstraints`"""

        self.phases_dims = n_phases * [AcadosOcpDims()]

        self.dummy_ocp_list = []

        # NOTE: this is the same for AcadosOcp
        self.solver_options = AcadosOcpOptions()
        """Solver Options, type :py:class:`acados_template.acados_ocp.AcadosOcpOptions`"""

        acados_path = get_acados_path()

        self.acados_include_path = os.path.join(acados_path, 'include').replace(os.sep, '/') # the replace part is important on Windows for CMake
        """Path to acados include directory (set automatically), type: `string`"""
        self.acados_lib_path = os.path.join(acados_path, 'lib').replace(os.sep, '/') # the replace part is important on Windows for CMake
        """Path to where acados library is located, type: `string`"""
        self.shared_lib_ext = get_shared_lib_ext()

        # get cython paths
        from sysconfig import get_paths
        self.cython_include_dirs = [np.get_include(), get_paths()['include']]

        self.__parameter_values = np.array([])

        self.code_export_directory = 'c_generated_code'
        """Path to where code will be exported. Default: `c_generated_code`."""

    @property
    def parameter_values(self):
        """:math:`p` - initial values for parameter vector - can be updated stagewise"""
        return self.__parameter_values

    @parameter_values.setter
    def parameter_values(self, parameter_values):
        if isinstance(parameter_values, np.ndarray):
            self.__parameter_values = parameter_values
        else:
            raise Exception('Invalid parameter_values value. ' +
                            f'Expected numpy array, got {type(parameter_values)}.')


    def set_phase(self, ocp: AcadosOcp, phase_idx: int) -> None:
        self.model[phase_idx] = ocp.model
        self.cost[phase_idx] = ocp.cost
        self.constraints[phase_idx] = ocp.constraints
        return

    def make_consistent(self) -> None:

        dims = self.dims
        opts = self.solver_options

        N_horizon = sum(self.N_list)
        if dims.N is None:
            dims.N = sum(self.N_list)
            print("AcadosMultiphaseOcp: make_consistent: N =", dims.N)
        elif dims.N != N_horizon:
            raise Exception(f"AcadosMultiphaseOcp: make_consistent: N = {dims.N} != {N_horizon} = sum(N_list).")

        phase_idx = np.cumsum([0] + self.N_list).tolist()

        self.start_idx = phase_idx[:-1]
        self.end_idx = phase_idx[1:]

        self.cost_start_idx = phase_idx.copy()
        self.cost_start_idx[0] += 1

        for i in range(self.n_phases):

            # create dummy ocp
            ocp = AcadosOcp()
            ocp.dims = self.phases_dims[i]
            ocp.dims.N = dims.N # NOTE: to not change options when making ocp consistent
            ocp.model = self.model[i]
            ocp.constraints = self.constraints[i]
            ocp.cost = self.cost[i]
            ocp.solver_options = opts

            print(f"Calling make_consistent for phase {i}.")
            ocp.make_consistent()

            self.dummy_ocp_list.append(ocp)

        return


    def to_dict(self) -> dict:
        # Copy ocp object dictionary
        ocp_dict = dict(deepcopy(self).__dict__)
        del ocp_dict['dummy_ocp_list']

        # convert acados classes to dicts
        for key, v in ocp_dict.items():
            if isinstance(v, (AcadosOcpOptions, MultiphaseOcpDims)):
                ocp_dict[key]=dict(getattr(self, key).__dict__)
            if isinstance(v, list):
                for i, item in enumerate(v):
                    if isinstance(item, (AcadosModel, AcadosOcpDims, AcadosOcpConstraints, AcadosOcpCost)):
                        ocp_dict[key][i] = format_class_dict(dict(item.__dict__))

        ocp_dict = format_class_dict(ocp_dict)
        return ocp_dict


    def remove_x0_elimination(self) -> None:
        self.constraints[0].idxbxe_0 = np.zeros((0,))
        self.dims[0].nbxe_0 = 0
        self.constraints[0].__has_x0 = False
        return


    def translate_nls_cost_to_conl(self):
        """
        Translates a NONLINEAR_LS cost to a CONVEX_OVER_NONLINEAR cost.
        """
        # TODO: make this functionality more generic, removing it from AcadosOcp
        raise NotImplementedError("translate_nls_cost_to_conl not implemented for multiphase OCPs.")
        return


    def formulate_constraint_as_L2_penalty(
        self,
        constr_expr: ca.SX,
        weight: float,
        upper_bound: Optional[float],
        lower_bound: Optional[float],
        residual_name: str = "new_residual",
    ) -> None:
        """
        Formulate a constraint as an L2 penalty and add it to the current cost.
        """
        # TODO: make this functionality more generic, removing it from AcadosOcp
        raise NotImplementedError("L2 penalty reformulation not implemented for multiphase OCPs.")
        return


    def formulate_constraint_as_Huber_penalty(
        self,
        constr_expr: ca.SX,
        weight: float,
        upper_bound: Optional[float],
        lower_bound: Optional[float],
        residual_name: str = "new_residual",
        huber_delta: float = 1.0,
        use_xgn = True,
        min_hess = 0,
    ) -> None:
        # TODO: make this functionality more generic, removing it from AcadosOcp
        raise NotImplementedError("Huber penalty reformulation not implemented for multiphase OCPs.")
        return