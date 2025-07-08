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

from typing import Optional, Union, Tuple
import numpy as np

from scipy.linalg import block_diag
from copy import deepcopy

import casadi as ca
import os, shutil
import json

from .acados_model import AcadosModel
from .acados_ocp_cost import AcadosOcpCost
from .acados_ocp_constraints import AcadosOcpConstraints
from .acados_dims import AcadosOcpDims
from .acados_ocp_options import AcadosOcpOptions
from .acados_ocp_iterate import AcadosOcpIterate

from .utils import (get_acados_path, format_class_dict, make_object_json_dumpable, render_template,
                    get_shared_lib_ext, is_column, is_empty, casadi_length, check_if_square,
                    check_casadi_version, ACADOS_INFTY)
from .penalty_utils import symmetric_huber_penalty, one_sided_huber_penalty

from .zoro_description import ZoroDescription, process_zoro_description
from .casadi_function_generation import (
    GenerateContext, AcadosCodegenOptions,
    generate_c_code_conl_cost, generate_c_code_nls_cost, generate_c_code_external_cost,
    generate_c_code_explicit_ode, generate_c_code_implicit_ode, generate_c_code_discrete_dynamics, generate_c_code_gnsf,
    generate_c_code_constraint
)

class AcadosOcp:
    """
    Class containing the full description of the optimal control problem.
    This object can be used to create an :py:class:`acados_template.acados_ocp_solver.AcadosOcpSolver`.

    The class has the following properties that can be modified to formulate a specific OCP, see below:

        - :py:attr:`dims` of type :py:class:`acados_template.acados_dims.AcadosOcpDims`
        - :py:attr:`model` of type :py:class:`acados_template.acados_model.AcadosModel`
        - :py:attr:`cost` of type :py:class:`acados_template.acados_ocp_cost.AcadosOcpCost`
        - :py:attr:`constraints` of type :py:class:`acados_template.acados_ocp_constraints.AcadosOcpConstraints`
        - :py:attr:`solver_options` of type :py:class:`acados_template.acados_ocp_options.AcadosOcpOptions`

        - :py:attr:`acados_include_path` (set automatically)
        - :py:attr:`shared_lib_ext` (set automatically)
        - :py:attr:`acados_lib_path` (set automatically)
        - :py:attr:`parameter_values` - used to initialize the parameters (can be changed)
        - :py:attr:`p_global_values` - used to initialize the global parameters (can be changed)
    """
    def __init__(self, acados_path=''):
        """
        Keyword arguments:
        acados_path -- path of your acados installation
        """
        if acados_path == '':
            acados_path = get_acados_path()

        self.dims = AcadosOcpDims()
        """Dimension definitions, type :py:class:`acados_template.acados_dims.AcadosOcpDims`"""
        self.model = AcadosModel()
        """Model definitions, type :py:class:`acados_template.acados_model.AcadosModel`"""
        self.cost = AcadosOcpCost()
        """Cost definitions, type :py:class:`acados_template.acados_ocp_cost.AcadosOcpCost`"""
        self.constraints = AcadosOcpConstraints()
        """Constraints definitions, type :py:class:`acados_template.acados_ocp_constraints.AcadosOcpConstraints`"""
        self.solver_options = AcadosOcpOptions()
        """Solver Options, type :py:class:`acados_template.acados_ocp_options.AcadosOcpOptions`"""

        self.zoro_description = None
        """zoRO - zero order robust optimization - description: for advanced users."""

        self.acados_include_path = os.path.join(acados_path, 'include').replace(os.sep, '/') # the replace part is important on Windows for CMake
        """Path to acados include directory (set automatically), type: `string`"""
        self.acados_lib_path = os.path.join(acados_path, 'lib').replace(os.sep, '/') # the replace part is important on Windows for CMake
        """Path to where acados library is located, type: `string`"""
        self.shared_lib_ext = get_shared_lib_ext()

        # get cython paths
        from sysconfig import get_paths
        self.cython_include_dirs = [np.get_include(), get_paths()['include']]

        self.__parameter_values = np.array([])
        self.__p_global_values = np.array([])
        self.__problem_class = 'OCP'
        self.__json_file = "acados_ocp.json"

        self.code_export_directory = 'c_generated_code'
        """Path to where code will be exported. Default: `c_generated_code`."""

        self.simulink_opts = None
        """Options to configure Simulink S-function blocks, mainly to activate possible Inputs and Outputs."""


    @property
    def parameter_values(self):
        """:math:`p` - initial values for parameter vector - can be updated stagewise"""
        return self.__parameter_values

    @parameter_values.setter
    def parameter_values(self, parameter_values):
        if isinstance(parameter_values, np.ndarray):
            if not is_column(parameter_values):
                raise ValueError("parameter_values should be column vector.")
            self.__parameter_values = parameter_values
        else:
            raise ValueError('Invalid parameter_values value. ' +
                            f'Expected numpy array, got {type(parameter_values)}.')

    @property
    def p_global_values(self):
        r"""initial values for :math:`p_\text{global}` vector, see `AcadosModel.p_global` - can be updated.
        Type: `numpy.ndarray` of shape `(np_global, )`.
        """
        return self.__p_global_values

    @p_global_values.setter
    def p_global_values(self, p_global_values):
        if isinstance(p_global_values, np.ndarray):
            if not is_column(p_global_values):
                raise ValueError("p_global_values should be column vector.")

            self.__p_global_values = p_global_values
        else:
            raise ValueError('Invalid p_global_values value. ' +
                            f'Expected numpy array, got {type(p_global_values)}.')

    @property
    def json_file(self):
        """Name of the json file where the problem description is stored."""
        return self.__json_file

    @json_file.setter
    def json_file(self, json_file):
        self.__json_file = json_file


    def _make_consistent_cost_initial(self):
        dims = self.dims
        cost = self.cost
        model = self.model
        opts = self.solver_options
        if opts.N_horizon == 0:
            return

        # if not set, copy fields from path constraints
        if cost.cost_type_0 is None:
            self.copy_path_cost_to_stage_0()

        if cost.cost_type_0 == 'AUTO':
            self.detect_cost_type(model, cost, dims, "initial")

        if cost.cost_type_0 in ['LINEAR_LS', 'NONLINEAR_LS']:
            if isinstance(cost.yref_0, (ca.SX, ca.MX, ca.DM)):
                raise Exception("yref_0 should be numpy array, symbolics are only supported before solver creation, to allow reformulating costs, e.g. using translate_cost_to_external_cost().")

            if isinstance(cost.W_0, (ca.SX, ca.MX, ca.DM)):
                raise Exception("W_0 should be numpy array, symbolics are only supported before solver creation, to allow reformulating costs, e.g. using translate_cost_to_external_cost().")

        if cost.cost_type_0 == 'LINEAR_LS':
            check_if_square(cost.W_0, 'W_0')
            ny_0 = cost.W_0.shape[0]
            if cost.Vx_0.shape[0] != ny_0 or cost.Vu_0.shape[0] != ny_0:
                raise ValueError('inconsistent dimension ny_0, regarding W_0, Vx_0, Vu_0.' + \
                                f'\nGot W_0[{cost.W_0.shape}], Vx_0[{cost.Vx_0.shape}], Vu_0[{cost.Vu_0.shape}]\n')
            if dims.nz != 0 and cost.Vz_0.shape[0] != ny_0:
                raise ValueError('inconsistent dimension ny_0, regarding W_0, Vx_0, Vu_0, Vz_0.' +
                                f'\nGot W_0[{cost.W_0.shape}], Vx_0[{cost.Vx_0.shape}], Vu_0[{cost.Vu_0.shape}], Vz_0[{cost.Vz_0.shape}]\n')
            if cost.Vx_0.shape[1] != dims.nx and ny_0 != 0:
                raise ValueError('inconsistent dimension: Vx_0 should have nx columns.')
            if cost.Vu_0.shape[1] != dims.nu and ny_0 != 0:
                raise ValueError('inconsistent dimension: Vu_0 should have nu columns.')
            if cost.yref_0.shape[0] != ny_0:
                raise ValueError('inconsistent dimension: regarding W_0, yref_0.' + \
                                f'\nGot W_0[{cost.W_0.shape}], yref_0[{cost.yref_0.shape}]\n')
            dims.ny_0 = ny_0

        elif cost.cost_type_0 == 'NONLINEAR_LS':
            ny_0 = cost.W_0.shape[0]
            check_if_square(cost.W_0, 'W_0')
            if (is_empty(model.cost_y_expr_0) and ny_0 != 0) or casadi_length(model.cost_y_expr_0) != ny_0 or cost.yref_0.shape[0] != ny_0:
                raise ValueError('inconsistent dimension ny_0: regarding W_0, cost_y_expr.' +
                                f'\nGot W_0[{cost.W_0.shape}], yref_0[{cost.yref_0.shape}], ',
                                f'cost_y_expr_0 [{casadi_length(model.cost_y_expr_0)}]\n')
            dims.ny_0 = ny_0

        elif cost.cost_type_0 == 'CONVEX_OVER_NONLINEAR':
            if is_empty(model.cost_y_expr_0):
                raise ValueError('cost_y_expr_0 and/or cost_y_expr not provided.')
            ny_0 = casadi_length(model.cost_y_expr_0)
            if is_empty(model.cost_r_in_psi_expr_0) or casadi_length(model.cost_r_in_psi_expr_0) != ny_0:
                raise ValueError('inconsistent dimension ny_0: regarding cost_y_expr_0 and cost_r_in_psi_0.')
            if is_empty(model.cost_psi_expr_0) or casadi_length(model.cost_psi_expr_0) != 1:
                raise ValueError('cost_psi_expr_0 not provided or not scalar-valued.')
            if cost.yref_0.shape[0] != ny_0:
                raise ValueError('inconsistent dimension: regarding yref_0 and cost_y_expr_0, cost_r_in_psi_0.')
            dims.ny_0 = ny_0

            if not (opts.hessian_approx=='EXACT' and opts.exact_hess_cost==False) and opts.hessian_approx != 'GAUSS_NEWTON':
                raise ValueError("\nWith CONVEX_OVER_NONLINEAR cost type, possible Hessian approximations are:\n"
                "GAUSS_NEWTON or EXACT with 'exact_hess_cost' == False.\n")

        elif cost.cost_type_0 == 'EXTERNAL':
            if isinstance(model.cost_expr_ext_cost_0, (float, int)):
                model.cost_expr_ext_cost_0 = ca.DM(model.cost_expr_ext_cost_0)
            if not isinstance(model.cost_expr_ext_cost_0, (ca.MX, ca.SX, ca.DM)):
                raise TypeError('cost_expr_ext_cost_0 should be casadi expression.')
            if not casadi_length(model.cost_expr_ext_cost_0) == 1:
                raise ValueError('cost_expr_ext_cost_0 should be scalar-valued.')
            if not is_empty(model.cost_expr_ext_cost_custom_hess_0):
                if model.cost_expr_ext_cost_custom_hess_0.shape != (dims.nx+dims.nu, dims.nx+dims.nu):
                    raise ValueError('cost_expr_ext_cost_custom_hess_0 should have shape (nx+nu, nx+nu).')


    def _make_consistent_cost_path(self):
        dims = self.dims
        cost = self.cost
        model = self.model
        opts = self.solver_options
        if opts.N_horizon == 0:
            return

        if cost.cost_type == 'AUTO':
            self.detect_cost_type(model, cost, dims, "path")

        if cost.cost_type in ['LINEAR_LS', 'NONLINEAR_LS']:
            if isinstance(cost.yref, (ca.SX, ca.MX, ca.DM)):
                raise Exception("yref should be numpy array, symbolics are only supported before solver creation, to allow reformulating costs, e.g. using translate_cost_to_external_cost().")

            if isinstance(cost.W, (ca.SX, ca.MX, ca.DM)):
                raise Exception("W should be numpy array, symbolics are only supported before solver creation, to allow reformulating costs, e.g. using translate_cost_to_external_cost().")

        if cost.cost_type == 'LINEAR_LS':
            ny = cost.W.shape[0]
            check_if_square(cost.W, 'W')
            if cost.Vx.shape[0] != ny or cost.Vu.shape[0] != ny:
                raise ValueError('inconsistent dimension ny, regarding W, Vx, Vu.' + \
                                f'\nGot W[{cost.W.shape}], Vx[{cost.Vx.shape}], Vu[{cost.Vu.shape}]\n')
            if dims.nz != 0 and cost.Vz.shape[0] != ny:
                raise ValueError('inconsistent dimension ny, regarding W, Vx, Vu, Vz.' + \
                                f'\nGot W[{cost.W.shape}], Vx[{cost.Vx.shape}], Vu[{cost.Vu.shape}], Vz[{cost.Vz.shape}]\n')
            if cost.Vx.shape[1] != dims.nx and ny != 0:
                raise ValueError('inconsistent dimension: Vx should have nx columns.')
            if cost.Vu.shape[1] != dims.nu and ny != 0:
                raise ValueError('inconsistent dimension: Vu should have nu columns.')
            if cost.yref.shape[0] != ny:
                raise ValueError('inconsistent dimension: regarding W, yref.' + \
                                f'\nGot W[{cost.W.shape}], yref[{cost.yref.shape}]\n')
            dims.ny = ny

        elif cost.cost_type == 'NONLINEAR_LS':
            ny = cost.W.shape[0]
            check_if_square(cost.W, 'W')
            if (is_empty(model.cost_y_expr) and ny != 0) or casadi_length(model.cost_y_expr) != ny or cost.yref.shape[0] != ny:
                raise ValueError('inconsistent dimension: regarding W, yref.' + \
                                f'\nGot W[{cost.W.shape}], yref[{cost.yref.shape}],',
                                f'cost_y_expr[{casadi_length(model.cost_y_expr)}]\n')
            dims.ny = ny

        elif cost.cost_type == 'CONVEX_OVER_NONLINEAR':
            if is_empty(model.cost_y_expr):
                raise ValueError('cost_y_expr and/or cost_y_expr not provided.')
            ny = casadi_length(model.cost_y_expr)
            if is_empty(model.cost_r_in_psi_expr) or casadi_length(model.cost_r_in_psi_expr) != ny:
                raise ValueError('inconsistent dimension ny: regarding cost_y_expr and cost_r_in_psi.')
            if is_empty(model.cost_psi_expr) or casadi_length(model.cost_psi_expr) != 1:
                raise ValueError('cost_psi_expr not provided or not scalar-valued.')
            if cost.yref.shape[0] != ny:
                raise ValueError('inconsistent dimension: regarding yref and cost_y_expr, cost_r_in_psi.')
            dims.ny = ny

            if not (opts.hessian_approx=='EXACT' and opts.exact_hess_cost==False) and opts.hessian_approx != 'GAUSS_NEWTON':
                raise ValueError("\nWith CONVEX_OVER_NONLINEAR cost type, possible Hessian approximations are:\n"
                "GAUSS_NEWTON or EXACT with 'exact_hess_cost' == False.\n")

        elif cost.cost_type == 'EXTERNAL':
            if isinstance(model.cost_expr_ext_cost, (float, int)):
                model.cost_expr_ext_cost = ca.DM(model.cost_expr_ext_cost)
            if not isinstance(model.cost_expr_ext_cost, (ca.MX, ca.SX, ca.DM)):
                raise TypeError('cost_expr_ext_cost should be casadi expression.')
            if not casadi_length(model.cost_expr_ext_cost) == 1:
                raise ValueError('cost_expr_ext_cost should be scalar-valued.')
            if not is_empty(model.cost_expr_ext_cost_custom_hess):
                if model.cost_expr_ext_cost_custom_hess.shape != (dims.nx+dims.nu, dims.nx+dims.nu):
                    raise ValueError('cost_expr_ext_cost_custom_hess should have shape (nx+nu, nx+nu).')


    def _make_consistent_cost_terminal(self):
        dims = self.dims
        cost = self.cost
        model = self.model
        opts = self.solver_options

        if cost.cost_type_e == 'AUTO':
            self.detect_cost_type(model, cost, dims, "terminal")

        if cost.cost_type_e in ['LINEAR_LS', 'NONLINEAR_LS']:
            if isinstance(cost.yref_e, (ca.SX, ca.MX, ca.DM)):
                raise Exception("yref_e should be numpy array, symbolics are only supported before solver creation, to allow reformulating costs, e.g. using translate_cost_to_external_cost().")

            if isinstance(cost.W_e, (ca.SX, ca.MX, ca.DM)):
                raise Exception("W_e should be numpy array, symbolics are only supported before solver creation, to allow reformulating costs, e.g. using translate_cost_to_external_cost().")

            ny_e = cost.W_e.shape[0]
            check_if_square(cost.W_e, 'W_e')
            dims.ny_e = ny_e

            if cost.cost_type_e == 'LINEAR_LS':
                if cost.Vx_e.shape[0] != ny_e:
                    raise ValueError('inconsistent dimension ny_e: regarding W_e, cost_y_expr_e.' + \
                        f'\nGot W_e[{cost.W_e.shape}], Vx_e[{cost.Vx_e.shape}]')
                if cost.Vx_e.shape[1] != dims.nx and ny_e != 0:
                    raise ValueError('inconsistent dimension: Vx_e should have nx columns.')
                if cost.yref_e.shape[0] != ny_e:
                    raise ValueError('inconsistent dimension: regarding W_e, yref_e.')

            elif cost.cost_type_e == 'NONLINEAR_LS':
                if (is_empty(model.cost_y_expr_e) and ny_e != 0) or casadi_length(model.cost_y_expr_e) != ny_e or cost.yref_e.shape[0] != ny_e:
                    raise ValueError('inconsistent dimension ny_e: regarding W_e, cost_y_expr.' +
                                    f'\nGot W_e[{cost.W_e.shape}], yref_e[{cost.yref_e.shape}], ',
                                    f'cost_y_expr_e [{casadi_length(model.cost_y_expr_e)}]\n')

        elif cost.cost_type_e == 'CONVEX_OVER_NONLINEAR':
            if is_empty(model.cost_y_expr_e):
                raise ValueError('cost_y_expr_e not provided.')
            ny_e = casadi_length(model.cost_y_expr_e)
            if is_empty(model.cost_r_in_psi_expr_e) or casadi_length(model.cost_r_in_psi_expr_e) != ny_e:
                raise ValueError('inconsistent dimension ny_e: regarding cost_y_expr_e and cost_r_in_psi_e.')
            if is_empty(model.cost_psi_expr_e) or casadi_length(model.cost_psi_expr_e) != 1:
                raise ValueError('cost_psi_expr_e not provided or not scalar-valued.')
            if cost.yref_e.shape[0] != ny_e:
                raise ValueError('inconsistent dimension: regarding yref_e and cost_y_expr_e, cost_r_in_psi_e.')
            dims.ny_e = ny_e

            if not (opts.hessian_approx=='EXACT' and opts.exact_hess_cost==False) and opts.hessian_approx != 'GAUSS_NEWTON':
                raise ValueError("\nWith CONVEX_OVER_NONLINEAR cost type, possible Hessian approximations are:\n"
                "GAUSS_NEWTON or EXACT with 'exact_hess_cost' == False.\n")

        elif cost.cost_type_e == 'EXTERNAL':
            if isinstance(model.cost_expr_ext_cost_e, (float, int)):
                model.cost_expr_ext_cost_e = ca.DM(model.cost_expr_ext_cost_e)
            if not isinstance(model.cost_expr_ext_cost_e, (ca.MX, ca.SX, ca.DM)):
                raise TypeError(f'cost_expr_ext_cost_e should be casadi expression, got {model.cost_expr_ext_cost_e}.')
            if not casadi_length(model.cost_expr_ext_cost_e) == 1:
                raise ValueError('cost_expr_ext_cost_e should be scalar-valued.')
            if not is_empty(model.cost_expr_ext_cost_custom_hess_e):
                if model.cost_expr_ext_cost_custom_hess_e.shape != (dims.nx, dims.nx):
                    raise ValueError('cost_expr_ext_cost_custom_hess_e should have shape (nx, nx).')


    def _make_consistent_constraints_initial(self):
        constraints = self.constraints
        dims = self.dims
        model = self.model
        opts = self.solver_options
        if opts.N_horizon == 0:
            dims.nbxe_0 = 0
            return

        nbx_0 = constraints.idxbx_0.shape[0]
        dims.nbx_0 = nbx_0
        if constraints.ubx_0.shape[0] != nbx_0 or constraints.lbx_0.shape[0] != nbx_0:
            raise ValueError('inconsistent dimension nbx_0, regarding idxbx_0, ubx_0, lbx_0.')
        if any(constraints.idxbx_0 >= dims.nx):
            raise ValueError(f'idxbx_0 = {constraints.idxbx_0} contains value >= nx = {dims.nx}.')

        if constraints.has_x0 and dims.nbx_0 != dims.nx:
            raise ValueError(f"x0 should have shape nx = {dims.nx}.")

        if constraints.has_x0 and not np.all(constraints.idxbxe_0 == np.arange(dims.nx)):
            raise ValueError(f"idxbxe_0 should be 0:{dims.nx} if x0 is set.")

        dims.nbxe_0 = constraints.idxbxe_0.shape[0]
        if any(constraints.idxbxe_0 >= dims.nbx_0):
            raise ValueError(f'idxbxe_0 = {constraints.idxbxe_0} contains value >= nbx_0 = {dims.nbx_0}.')

        nh_0 = 0 if is_empty(model.con_h_expr_0) else casadi_length(model.con_h_expr_0)

        if constraints.uh_0.shape[0] != nh_0 or constraints.lh_0.shape[0] != nh_0:
            raise ValueError('inconsistent dimension nh_0, regarding lh_0, uh_0, con_h_expr_0.')
        else:
            dims.nh_0 = nh_0

        if is_empty(model.con_phi_expr_0):
            dims.nphi_0 = 0
            dims.nr_0 = 0
        else:
            dims.nphi_0 = casadi_length(model.con_phi_expr_0)
            constraints.constr_type_0 = "BGP"
            if is_empty(model.con_r_expr_0):
                raise ValueError('convex over nonlinear constraints: con_r_expr_0 but con_phi_expr_0 is nonempty')
            else:
                dims.nr_0 = casadi_length(model.con_r_expr_0)


    def _make_consistent_constraints_path(self):
        constraints = self.constraints
        dims = self.dims
        model = self.model
        opts = self.solver_options
        if opts.N_horizon == 0:
            return

        nbx = constraints.idxbx.shape[0]
        if constraints.ubx.shape[0] != nbx or constraints.lbx.shape[0] != nbx:
            raise ValueError('inconsistent dimension nbx, regarding idxbx, ubx, lbx.')
        else:
            dims.nbx = nbx
        if any(constraints.idxbx >= dims.nx):
            raise ValueError(f'idxbx = {constraints.idxbx} contains value >= nx = {dims.nx}.')

        nbu = constraints.idxbu.shape[0]
        if constraints.ubu.shape[0] != nbu or constraints.lbu.shape[0] != nbu:
            raise ValueError('inconsistent dimension nbu, regarding idxbu, ubu, lbu.')
        else:
            dims.nbu = nbu
        if any(constraints.idxbu >= dims.nu):
            raise ValueError(f'idxbu = {constraints.idxbu} contains value >= nu = {dims.nu}.')

        # lg <= C * x + D * u <= ug
        ng = constraints.lg.shape[0]
        if constraints.ug.shape[0] != ng or constraints.C.shape[0] != ng \
        or constraints.D.shape[0] != ng:
            raise ValueError('inconsistent dimension ng, regarding lg, ug, C, D.')
        else:
            dims.ng = ng

        if ng > 0:
            if constraints.C.shape[1] != dims.nx:
                raise ValueError(f'inconsistent dimension nx, regarding C, got C.shape[1] = {constraints.C.shape[1]}.')
            if constraints.D.shape[1] != dims.nu:
                raise ValueError(f'inconsistent dimension nu, regarding D, got D.shape[1] = {constraints.D.shape[1]}.')

        if not is_empty(model.con_h_expr):
            nh = casadi_length(model.con_h_expr)
        else:
            nh = 0

        if constraints.uh.shape[0] != nh or constraints.lh.shape[0] != nh:
            raise ValueError('inconsistent dimension nh, regarding lh, uh, con_h_expr.')
        else:
            dims.nh = nh

        if is_empty(model.con_phi_expr):
            dims.nphi = 0
            dims.nr = 0
        else:
            dims.nphi = casadi_length(model.con_phi_expr)
            constraints.constr_type = "BGP"
            if is_empty(model.con_r_expr):
                raise ValueError('convex over nonlinear constraints: con_r_expr but con_phi_expr is nonempty')
            else:
                dims.nr = casadi_length(model.con_r_expr)


    def _make_consistent_constraints_terminal(self):
        dims = self.dims
        constraints = self.constraints
        model = self.model

        nbx_e = constraints.idxbx_e.shape[0]
        if constraints.ubx_e.shape[0] != nbx_e or constraints.lbx_e.shape[0] != nbx_e:
            raise ValueError('inconsistent dimension nbx_e, regarding idxbx_e, ubx_e, lbx_e.')
        else:
            dims.nbx_e = nbx_e
        if any(constraints.idxbx_e >= dims.nx):
            raise ValueError(f'idxbx_e = {constraints.idxbx_e} contains value >= nx = {dims.nx}.')

        ng_e = constraints.lg_e.shape[0]
        if constraints.ug_e.shape[0] != ng_e or constraints.C_e.shape[0] != ng_e:
            raise ValueError('inconsistent dimension ng_e, regarding_e lg_e, ug_e, C_e.')
        else:
            dims.ng_e = ng_e

        nh_e = 0 if is_empty(model.con_h_expr_e) else casadi_length(model.con_h_expr_e)

        if constraints.uh_e.shape[0] != nh_e or constraints.lh_e.shape[0] != nh_e:
            raise ValueError('inconsistent dimension nh_e, regarding lh_e, uh_e, con_h_expr_e.')
        else:
            dims.nh_e = nh_e

        if is_empty(model.con_phi_expr_e):
            dims.nphi_e = 0
            dims.nr_e = 0
        else:
            dims.nphi_e = casadi_length(model.con_phi_expr_e)
            constraints.constr_type_e = "BGP"
            if is_empty(model.con_r_expr_e):
                raise ValueError('convex over nonlinear constraints: con_r_expr_e but con_phi_expr_e is nonempty')
            else:
                dims.nr_e = casadi_length(model.con_r_expr_e)


    def _make_consistent_slacks_initial(self):
        constraints = self.constraints
        dims = self.dims
        opts = self.solver_options
        cost = self.cost
        if opts.N_horizon == 0:
            return

        nh_0 = dims.nh_0
        nsbu = dims.nsbu
        nsg = dims.nsg
        ns = dims.ns
        nsh_0 = constraints.idxsh_0.shape[0]
        if nsh_0 > nh_0:
            raise ValueError(f'inconsistent dimension nsh_0 = {nsh_0}. Is greater than nh_0 = {nh_0}.')
        if any(constraints.idxsh_0 >= nh_0):
            raise ValueError(f'idxsh_0 = {constraints.idxsh_0} contains value >= nh_0 = {nh_0}.')
        if is_empty(constraints.lsh_0):
            constraints.lsh_0 = np.zeros((nsh_0,))
        elif constraints.lsh_0.shape[0] != nsh_0:
            raise ValueError('inconsistent dimension nsh_0, regarding idxsh_0, lsh_0.')
        if is_empty(constraints.ush_0):
            constraints.ush_0 = np.zeros((nsh_0,))
        elif constraints.ush_0.shape[0] != nsh_0:
            raise ValueError('inconsistent dimension nsh_0, regarding idxsh_0, ush_0.')
        dims.nsh_0 = nsh_0

        nsphi_0 = constraints.idxsphi_0.shape[0]
        if nsphi_0 > dims.nphi_0:
            raise ValueError(f'inconsistent dimension nsphi_0 = {nsphi_0}. Is greater than nphi_0 = {dims.nphi_0}.')
        if any(constraints.idxsphi_0 >= dims.nphi_0):
            raise ValueError(f'idxsphi_0 = {constraints.idxsphi_0} contains value >= nphi_0 = {dims.nphi_0}.')
        if is_empty(constraints.lsphi_0):
            constraints.lsphi_0 = np.zeros((nsphi_0,))
        elif constraints.lsphi_0.shape[0] != nsphi_0:
            raise ValueError('inconsistent dimension nsphi_0, regarding idxsphi_0, lsphi_0.')
        if is_empty(constraints.usphi_0):
            constraints.usphi_0 = np.zeros((nsphi_0,))
        elif constraints.usphi_0.shape[0] != nsphi_0:
            raise ValueError('inconsistent dimension nsphi_0, regarding idxsphi_0, usphi_0.')
        dims.nsphi_0 = nsphi_0

        # Note: at stage 0 bounds on x are not slacked!
        ns_0 = nsbu + nsg + nsphi_0 + nsh_0  # NOTE: nsbx not supported at stage 0

        if cost.zl_0 is None and cost.zu_0 is None and cost.Zl_0 is None and cost.Zu_0 is None:
            if ns_0 == 0:
                cost.zl_0 = np.array([])
                cost.zu_0 = np.array([])
                cost.Zl_0 = np.array([])
                cost.Zu_0 = np.array([])
            elif ns_0 == ns:
                cost.zl_0 = cost.zl
                cost.zu_0 = cost.zu
                cost.Zl_0 = cost.Zl
                cost.Zu_0 = cost.Zu
                print("Fields cost.[zl_0, zu_0, Zl_0, Zu_0] are not provided.")
                print("Using entries [zl, zu, Zl, Zu] at initial node for slack penalties.\n")
            else:
                raise ValueError("Fields cost.[zl_0, zu_0, Zl_0, Zu_0] are not provided and cannot be inferred from other fields.\n")

        for field in ("Zl_0", "Zu_0", "zl_0", "zu_0"):
            dim = getattr(cost, field).shape[0]
            if dim != ns_0:
                raise Exception(f'Inconsistent size for fields {field}, with dimension {dim}, \n\t'\
                + f'Detected ns_0 = {ns_0} = nsbu + nsg + nsh_0 + nsphi_0.\n\t'\
                + f'With nsbu = {nsbu}, nsg = {nsg}, nsh_0 = {nsh_0}, nsphi_0 = {nsphi_0}.')
        dims.ns_0 = ns_0


    def _make_consistent_slacks_path(self):
        constraints = self.constraints
        dims = self.dims
        opts = self.solver_options
        cost = self.cost
        if opts.N_horizon == 0:
            return

        nbx = dims.nbx
        nbu = dims.nbu
        nh = dims.nh
        ng = dims.ng
        nsbx = constraints.idxsbx.shape[0]
        if nsbx > nbx:
            raise ValueError(f'inconsistent dimension nsbx = {nsbx}. Is greater than nbx = {nbx}.')
        if any(constraints.idxsbx >= nbx):
            raise ValueError(f'idxsbx = {constraints.idxsbx} contains value >= nbx = {nbx}.')
        if is_empty(constraints.lsbx):
            constraints.lsbx = np.zeros((nsbx,))
        elif constraints.lsbx.shape[0] != nsbx:
            raise ValueError('inconsistent dimension nsbx, regarding idxsbx, lsbx.')
        if is_empty(constraints.usbx):
            constraints.usbx = np.zeros((nsbx,))
        elif constraints.usbx.shape[0] != nsbx:
            raise ValueError('inconsistent dimension nsbx, regarding idxsbx, usbx.')
        dims.nsbx = nsbx

        nsbu = constraints.idxsbu.shape[0]
        if nsbu > nbu:
            raise ValueError(f'inconsistent dimension nsbu = {nsbu}. Is greater than nbu = {nbu}.')
        if any(constraints.idxsbu >= nbu):
            raise ValueError(f'idxsbu = {constraints.idxsbu} contains value >= nbu = {nbu}.')
        if is_empty(constraints.lsbu):
            constraints.lsbu = np.zeros((nsbu,))
        elif constraints.lsbu.shape[0] != nsbu:
            raise ValueError('inconsistent dimension nsbu, regarding idxsbu, lsbu.')
        if is_empty(constraints.usbu):
            constraints.usbu = np.zeros((nsbu,))
        elif constraints.usbu.shape[0] != nsbu:
            raise ValueError('inconsistent dimension nsbu, regarding idxsbu, usbu.')
        dims.nsbu = nsbu

        nsh = constraints.idxsh.shape[0]
        if nsh > nh:
            raise ValueError(f'inconsistent dimension nsh = {nsh}. Is greater than nh = {nh}.')
        if any(constraints.idxsh >= nh):
            raise ValueError(f'idxsh = {constraints.idxsh} contains value >= nh = {nh}.')
        if is_empty(constraints.lsh):
            constraints.lsh = np.zeros((nsh,))
        elif constraints.lsh.shape[0] != nsh:
            raise ValueError('inconsistent dimension nsh, regarding idxsh, lsh.')
        if is_empty(constraints.ush):
            constraints.ush = np.zeros((nsh,))
        elif constraints.ush.shape[0] != nsh:
            raise ValueError('inconsistent dimension nsh, regarding idxsh, ush.')
        dims.nsh = nsh

        nsphi = constraints.idxsphi.shape[0]
        if nsphi > dims.nphi:
            raise ValueError(f'inconsistent dimension nsphi = {nsphi}. Is greater than nphi = {dims.nphi}.')
        if any(constraints.idxsphi >= dims.nphi):
            raise ValueError(f'idxsphi = {constraints.idxsphi} contains value >= nphi = {dims.nphi}.')
        if is_empty(constraints.lsphi):
            constraints.lsphi = np.zeros((nsphi,))
        elif constraints.lsphi.shape[0] != nsphi:
            raise ValueError('inconsistent dimension nsphi, regarding idxsphi, lsphi.')
        if is_empty(constraints.usphi):
            constraints.usphi = np.zeros((nsphi,))
        elif constraints.usphi.shape[0] != nsphi:
            raise ValueError('inconsistent dimension nsphi, regarding idxsphi, usphi.')
        dims.nsphi = nsphi

        nsg = constraints.idxsg.shape[0]
        if nsg > ng:
            raise ValueError(f'inconsistent dimension nsg = {nsg}. Is greater than ng = {ng}.')
        if any(constraints.idxsg >= ng):
            raise ValueError(f'idxsg = {constraints.idxsg} contains value >= ng = {ng}.')
        if is_empty(constraints.lsg):
            constraints.lsg = np.zeros((nsg,))
        elif constraints.lsg.shape[0] != nsg:
            raise ValueError('inconsistent dimension nsg, regarding idxsg, lsg.')
        if is_empty(constraints.usg):
            constraints.usg = np.zeros((nsg,))
        elif constraints.usg.shape[0] != nsg:
            raise ValueError('inconsistent dimension nsg, regarding idxsg, usg.')
        dims.nsg = nsg

        ns = nsbx + nsbu + nsh + nsg + nsphi
        for field in ("Zl", "Zu", "zl", "zu"):
            dim = getattr(cost, field).shape[0]
            if dim != ns:
                raise Exception(f'Inconsistent size for fields {field}, with dimension {dim}, \n\t'\
                    + f'Detected ns = {ns} = nsbx + nsbu + nsg + nsh + nsphi.\n\t'\
                    + f'With nsbx = {nsbx}, nsbu = {nsbu}, nsg = {nsg}, nsh = {nsh}, nsphi = {nsphi}.')
        dims.ns = ns


    def _make_consistent_slacks_terminal(self):
        constraints = self.constraints
        dims = self.dims
        cost = self.cost

        nbx_e = dims.nbx_e
        nh_e = dims.nh_e
        ng_e = dims.ng_e
        nsbx_e = constraints.idxsbx_e.shape[0]
        if nsbx_e > nbx_e:
            raise ValueError(f'inconsistent dimension nsbx_e = {nsbx_e}. Is greater than nbx_e = {nbx_e}.')
        if any(constraints.idxsbx_e >= nbx_e):
            raise ValueError(f'idxsbx_e = {constraints.idxsbx_e} contains value >= nbx_e = {nbx_e}.')
        if is_empty(constraints.lsbx_e):
            constraints.lsbx_e = np.zeros((nsbx_e,))
        elif constraints.lsbx_e.shape[0] != nsbx_e:
            raise ValueError('inconsistent dimension nsbx_e, regarding idxsbx_e, lsbx_e.')
        if is_empty(constraints.usbx_e):
            constraints.usbx_e = np.zeros((nsbx_e,))
        elif constraints.usbx_e.shape[0] != nsbx_e:
            raise ValueError('inconsistent dimension nsbx_e, regarding idxsbx_e, usbx_e.')
        dims.nsbx_e = nsbx_e

        nsh_e = constraints.idxsh_e.shape[0]
        if nsh_e > nh_e:
            raise ValueError(f'inconsistent dimension nsh_e = {nsh_e}. Is greater than nh_e = {nh_e}.')
        if any(constraints.idxsh_e >= nh_e):
            raise ValueError(f'idxsh_e = {constraints.idxsh_e} contains value >= nh_e = {nh_e}.')
        if is_empty(constraints.lsh_e):
            constraints.lsh_e = np.zeros((nsh_e,))
        elif constraints.lsh_e.shape[0] != nsh_e:
            raise ValueError('inconsistent dimension nsh_e, regarding idxsh_e, lsh_e.')
        if is_empty(constraints.ush_e):
            constraints.ush_e = np.zeros((nsh_e,))
        elif constraints.ush_e.shape[0] != nsh_e:
            raise ValueError('inconsistent dimension nsh_e, regarding idxsh_e, ush_e.')
        dims.nsh_e = nsh_e

        nsphi_e = constraints.idxsphi_e.shape[0]
        if nsphi_e > dims.nphi_e:
            raise ValueError(f'inconsistent dimension nsphi_e = {nsphi_e}. Is greater than nphi_e = {dims.nphi_e}.')
        if any(constraints.idxsphi_e >= dims.nphi_e):
            raise ValueError(f'idxsphi_e = {constraints.idxsphi_e} contains value >= nphi_e = {dims.nphi_e}.')
        if is_empty(constraints.lsphi_e):
            constraints.lsphi_e = np.zeros((nsphi_e,))
        elif constraints.lsphi_e.shape[0] != nsphi_e:
            raise ValueError('inconsistent dimension nsphi_e, regarding idxsphi_e, lsphi_e.')
        if is_empty(constraints.usphi_e):
            constraints.usphi_e = np.zeros((nsphi_e,))
        elif constraints.usphi_e.shape[0] != nsphi_e:
            raise ValueError('inconsistent dimension nsphi_e, regarding idxsphi_e, usphi_e.')
        dims.nsphi_e = nsphi_e

        nsg_e = constraints.idxsg_e.shape[0]
        if nsg_e > ng_e:
            raise ValueError(f'inconsistent dimension nsg_e = {nsg_e}. Is greater than ng_e = {ng_e}.')
        if any(constraints.idxsg_e >= ng_e):
            raise ValueError(f'idxsg_e = {constraints.idxsg_e} contains value >= ng_e = {ng_e}.')
        if is_empty(constraints.lsg_e):
            constraints.lsg_e = np.zeros((nsg_e,))
        elif constraints.lsg_e.shape[0] != nsg_e:
            raise ValueError('inconsistent dimension nsg_e, regarding idxsg_e, lsg_e.')
        if is_empty(constraints.usg_e):
            constraints.usg_e = np.zeros((nsg_e,))
        elif constraints.usg_e.shape[0] != nsg_e:
            raise ValueError('inconsistent dimension nsg_e, regarding idxsg_e, usg_e.')
        dims.nsg_e = nsg_e

        # terminal
        ns_e = nsbx_e + nsh_e + nsg_e + nsphi_e
        for field in ("Zl_e", "Zu_e", "zl_e", "zu_e"):
            dim = getattr(cost, field).shape[0]
            if dim != ns_e:
                raise Exception(f'Inconsistent size for fields {field}, with dimension {dim}, \n\t'\
                + f'Detected ns_e = {ns_e} = nsbx_e + nsg_e + nsh_e + nsphi_e.\n\t'\
                + f'With nsbx_e = {nsbx_e}, nsg_e = {nsg_e}, nsh_e = {nsh_e}, nsphi_e = {nsphi_e}.')

        dims.ns_e = ns_e


    def _make_consistent_discretization(self):
        opts = self.solver_options
        if opts.N_horizon == 0:
            opts.shooting_nodes = np.array([0.])
            opts.time_steps = np.array([])
            return

        if not isinstance(opts.tf, (float, int)):
            raise TypeError(f'Time horizon tf should be float provided, got tf = {opts.tf}.')

        if is_empty(opts.time_steps) and is_empty(opts.shooting_nodes):
            # uniform discretization
            opts.time_steps = opts.tf / opts.N_horizon * np.ones((opts.N_horizon,))
            opts.shooting_nodes = np.concatenate((np.array([0.]), np.cumsum(opts.time_steps)))

        elif not is_empty(opts.shooting_nodes):
            if np.shape(opts.shooting_nodes)[0] != opts.N_horizon+1:
                raise ValueError('inconsistent dimension N, regarding shooting_nodes.')

            time_steps = opts.shooting_nodes[1:] - opts.shooting_nodes[0:-1]
            # identify constant time_steps: due to numerical reasons the content of time_steps might vary a bit
            avg_time_steps = np.average(time_steps)
            # criterion for constant time step detection: the min/max difference in values normalized by the average
            check_const_time_step = (np.max(time_steps)-np.min(time_steps)) / avg_time_steps
            # if the criterion is small, we have a constant time_step
            if check_const_time_step < 1e-9:
                time_steps[:] = avg_time_steps  # if we have a constant time_step: apply the average time_step

            opts.time_steps = time_steps

        elif not is_empty(opts.time_steps) and is_empty(opts.shooting_nodes):
            # compute shooting nodes from time_steps for convenience
            opts.shooting_nodes = np.concatenate((np.array([0.]), np.cumsum(opts.time_steps)))

        elif (not is_empty(opts.time_steps)) and (not is_empty(opts.shooting_nodes)):
            ValueError('Please provide either time_steps or shooting_nodes for nonuniform discretization')

        tf = np.sum(opts.time_steps)
        if (tf - opts.tf) / tf > 1e-13:
            raise ValueError(f'Inconsistent discretization: {opts.tf}'
                f' = tf != sum(opts.time_steps) = {tf}.')


    def _make_consistent_simulation(self):
        opts = self.solver_options
        if opts.N_horizon == 0:
            return

        # set integrator time automatically
        opts.Tsim = opts.time_steps[0]

        # num_steps
        if isinstance(opts.sim_method_num_steps, np.ndarray) and opts.sim_method_num_steps.size == 1:
            opts.sim_method_num_steps = opts.sim_method_num_steps.item()

        if isinstance(opts.sim_method_num_steps, (int, float)) and opts.sim_method_num_steps % 1 == 0:
            opts.sim_method_num_steps = opts.sim_method_num_steps * np.ones((opts.N_horizon,), dtype=np.int64)
        elif isinstance(opts.sim_method_num_steps, np.ndarray) and opts.sim_method_num_steps.size == opts.N_horizon \
            and np.all(np.equal(np.mod(opts.sim_method_num_steps, 1), 0)):
            opts.sim_method_num_steps = np.reshape(opts.sim_method_num_steps, (opts.N_horizon,)).astype(np.int64)
        else:
            raise TypeError("Wrong value for sim_method_num_steps. Should be either int or array of ints of shape (N,).")

        # num_stages
        if isinstance(opts.sim_method_num_stages, np.ndarray) and opts.sim_method_num_stages.size == 1:
            opts.sim_method_num_stages = opts.sim_method_num_stages.item()

        if isinstance(opts.sim_method_num_stages, (int, float)) and opts.sim_method_num_stages % 1 == 0:
            opts.sim_method_num_stages = opts.sim_method_num_stages * np.ones((opts.N_horizon,), dtype=np.int64)
        elif isinstance(opts.sim_method_num_stages, np.ndarray) and opts.sim_method_num_stages.size == opts.N_horizon \
            and np.all(np.equal(np.mod(opts.sim_method_num_stages, 1), 0)):
            opts.sim_method_num_stages = np.reshape(opts.sim_method_num_stages, (opts.N_horizon,)).astype(np.int64)
        else:
            raise ValueError("Wrong value for sim_method_num_stages. Should be either int or array of ints of shape (N,).")

        # jac_reuse
        if isinstance(opts.sim_method_jac_reuse, np.ndarray) and opts.sim_method_jac_reuse.size == 1:
            opts.sim_method_jac_reuse = opts.sim_method_jac_reuse.item()

        if isinstance(opts.sim_method_jac_reuse, (int, float)) and opts.sim_method_jac_reuse % 1 == 0:
            opts.sim_method_jac_reuse = opts.sim_method_jac_reuse * np.ones((opts.N_horizon,), dtype=np.int64)
        elif isinstance(opts.sim_method_jac_reuse, np.ndarray) and opts.sim_method_jac_reuse.size == opts.N_horizon \
            and np.all(np.equal(np.mod(opts.sim_method_jac_reuse, 1), 0)):
            opts.sim_method_jac_reuse = np.reshape(opts.sim_method_jac_reuse, (opts.N_horizon,)).astype(np.int64)
        else:
            raise ValueError("Wrong value for sim_method_jac_reuse. Should be either int or array of ints of shape (N,).")


    def make_consistent(self, is_mocp_phase: bool=False, verbose: bool=True) -> None:
        """
        Detect dimensions, perform sanity checks
        """
        dims = self.dims
        cost = self.cost
        constraints = self.constraints
        model = self.model
        opts = self.solver_options

        model.make_consistent(dims)
        self.name = model.name

        if opts.N_horizon is None and dims.N is None:
            raise ValueError('N_horizon not provided.')
        elif opts.N_horizon is None and dims.N is not None:
            opts.N_horizon = dims.N
            print("field AcadosOcpDims.N has been migrated to AcadosOcpOptions.N_horizon. setting AcadosOcpOptions.N_horizon = N. For future comppatibility, please use AcadosOcpOptions.N_horizon directly.")
        elif opts.N_horizon is not None and dims.N is not None and opts.N_horizon != dims.N:
            raise ValueError(f'Inconsistent dimension N, regarding N = {dims.N}, N_horizon = {opts.N_horizon}.')
        else:
            dims.N = opts.N_horizon

        # check if nx != nx_next
        if not is_mocp_phase and dims.nx != dims.nx_next and opts.N_horizon > 1:
            raise ValueError('nx_next should be equal to nx if more than one shooting interval is used.')

        # parameters
        if self.parameter_values.shape[0] != dims.np:
            raise ValueError('inconsistent dimension np, regarding model.p and parameter_values.' + \
                f'\nGot np = {dims.np}, self.parameter_values.shape = {self.parameter_values.shape[0]}\n')

        # p_global_values
        if self.p_global_values.shape[0] != dims.np_global:
            raise ValueError('inconsistent dimension np_global, regarding model.p_global and p_global_values.' + \
                f'\nGot np_global = {dims.np_global}, self.p_global_values.shape = {self.p_global_values.shape[0]}\n')

        ## cost
        self._make_consistent_cost_initial()
        self._make_consistent_cost_path()
        self._make_consistent_cost_terminal()

        # GN check
        if verbose:
            gn_warning_0 = (opts.N_horizon > 0 and cost.cost_type_0 == 'EXTERNAL' and opts.hessian_approx == 'GAUSS_NEWTON' and opts.ext_cost_num_hess == 0 and is_empty(model.cost_expr_ext_cost_custom_hess_0))
            gn_warning_path = (opts.N_horizon > 0 and cost.cost_type == 'EXTERNAL' and opts.hessian_approx == 'GAUSS_NEWTON' and opts.ext_cost_num_hess == 0 and is_empty(model.cost_expr_ext_cost_custom_hess))
            gn_warning_terminal = (cost.cost_type_e == 'EXTERNAL' and opts.hessian_approx == 'GAUSS_NEWTON' and opts.ext_cost_num_hess == 0 and is_empty(model.cost_expr_ext_cost_custom_hess_e))
            if any([gn_warning_0, gn_warning_path, gn_warning_terminal]):
                external_cost_types = []
                if gn_warning_0:
                    external_cost_types.append('cost_type_0')
                if gn_warning_path:
                    external_cost_types.append('cost_type')
                if gn_warning_terminal:
                    external_cost_types.append('cost_type_e')
                print("\nWARNING: Gauss-Newton Hessian approximation with EXTERNAL cost type not well defined!\n"
                f"got cost_type EXTERNAL for {', '.join(external_cost_types)}, hessian_approx: 'GAUSS_NEWTON'.\n"
                "With this setting, acados will proceed computing the exact Hessian for the cost term and no Hessian contribution from constraints and dynamics.\n"
                "If the external cost is a linear least squares cost, this coincides with the Gauss-Newton Hessian.\n"
                "Note: There is also the option to use the external cost module with a numerical Hessian approximation (see `ext_cost_num_hess`).\n"
                "OR the option to provide a symbolic custom Hessian approximation (see `cost_expr_ext_cost_custom_hess`).\n")

        # cost integration
        if opts.N_horizon > 0:
            supports_cost_integration = lambda type : type in ['NONLINEAR_LS', 'CONVEX_OVER_NONLINEAR']
            if opts.cost_discretization == 'INTEGRATOR':
                if any([not supports_cost_integration(cost) for cost in [cost.cost_type_0, cost.cost_type]]):
                    raise ValueError('cost_discretization == INTEGRATOR only works with cost in ["NONLINEAR_LS", "CONVEX_OVER_NONLINEAR"] costs.')
                if opts.nlp_solver_type == "SQP_WITH_FEASIBLE_QP":
                    raise ValueError('cost_discretization == INTEGRATOR is not compatible with SQP_WITH_FEASIBLE_QP yet.')

        ## constraints
        if opts.qp_solver == 'PARTIAL_CONDENSING_QPDUNES':
            self.remove_x0_elimination()
        self._make_consistent_constraints_initial()
        self._make_consistent_constraints_path()
        self._make_consistent_constraints_terminal()

        self._make_consistent_slacks_path()
        self._make_consistent_slacks_initial()
        self._make_consistent_slacks_terminal()

        # check for ACADOS_INFTY
        if opts.qp_solver not in ["PARTIAL_CONDENSING_HPIPM", "FULL_CONDENSING_HPIPM", "FULL_CONDENSING_DAQP"]:
            # loop over all bound vectors
            if opts.N_horizon > 0:
                fields_to_check = ['lbx_0', 'ubx_0', 'lbx', 'ubx', 'lbx_e', 'ubx_e', 'lg', 'ug', 'lg_e', 'ug_e', 'lh', 'uh', 'lh_e', 'uh_e', 'lbu', 'ubu', 'lphi', 'uphi', 'lphi_e', 'uphi_e']
            else:
                fields_to_check = ['lbx_0', 'ubx_0', 'lbx_e', 'ubx_e', 'lg_e', 'ug_e', 'lh_e', 'uh_e''lphi_e', 'uphi_e']
            for field in fields_to_check:
                bound = getattr(constraints, field)
                if any(bound >= ACADOS_INFTY) or any(bound <= -ACADOS_INFTY):
                    raise ValueError(f"Field {field} contains values outside the interval (-ACADOS_INFTY, ACADOS_INFTY) with ACADOS_INFTY = {ACADOS_INFTY:.2e}. One-sided constraints are not supported by the chosen QP solver {opts.qp_solver}.")

        self._make_consistent_discretization()

        # cost scaling
        if opts.cost_scaling is None:
            opts.cost_scaling = np.append(opts.time_steps, 1.0)
        if opts.cost_scaling.shape[0] != opts.N_horizon + 1:
            raise ValueError(f'cost_scaling should be of length N+1 = {opts.N_horizon+1}, got {opts.cost_scaling.shape[0]}.')

        self._make_consistent_simulation()

        # fixed hessian
        if opts.fixed_hess:
            if opts.hessian_approx == 'EXACT':
                raise ValueError('fixed_hess is not compatible with hessian_approx == EXACT.')
            if cost.cost_type != "LINEAR_LS" and opts.N_horizon > 0:
                raise ValueError('fixed_hess is only compatible LINEAR_LS cost_type.')
            if cost.cost_type_0 != "LINEAR_LS" and opts.N_horizon > 0:
                raise ValueError('fixed_hess is only compatible LINEAR_LS cost_type_0.')
            if cost.cost_type_e != "LINEAR_LS":
                raise ValueError('fixed_hess is only compatible LINEAR_LS cost_type_e.')

        # solution sensitivities
        if opts.N_horizon > 0:
            bgp_type_constraint_pairs = [
                ("path", model.con_phi_expr), ("initial", model.con_phi_expr_0), ("terminal", model.con_phi_expr_e),
                ("path", model.con_r_expr), ("initial", model.con_r_expr_0), ("terminal", model.con_r_expr_e)
            ]
            cost_types_to_check = [cost.cost_type, cost.cost_type_0, cost.cost_type_e]
            suffix = f", got cost_types {cost.cost_type_0, cost.cost_type, cost.cost_type_e}."
        else:
            bgp_type_constraint_pairs = [
                ("terminal", model.con_phi_expr_e), ("terminal", model.con_r_expr_e)
            ]
            cost_types_to_check = [cost.cost_type_e]
            suffix = f", got cost_type_e {cost.cost_type_e}."

        if opts.with_solution_sens_wrt_params:
            if dims.np_global == 0:
                raise ValueError('with_solution_sens_wrt_params is only compatible if global parameters `p_global` are provided. Sensitivities wrt parameters have been refactored to use p_global instead of p in https://github.com/acados/acados/pull/1316. Got emty p_global.')
            if any([cost_type not in ["EXTERNAL", "LINEAR_LS"] for cost_type in cost_types_to_check]):
                raise ValueError('with_solution_sens_wrt_params is only compatible with EXTERNAL and LINEAR_LS cost_type' + suffix)
            if opts.N_horizon > 0 and opts.integrator_type != "DISCRETE":
                raise NotImplementedError('with_solution_sens_wrt_params is only compatible with DISCRETE dynamics.')
            for horizon_type, constraint in bgp_type_constraint_pairs:
                if constraint is not None and any(ca.which_depends(constraint, model.p_global)):
                    raise NotImplementedError(f"with_solution_sens_wrt_params is not supported for BGP constraints that depend on p_global. Got dependency on p_global for {horizon_type} constraint.")
            if opts.qp_solver_cond_N != opts.N_horizon or opts.qp_solver.startswith("FULL_CONDENSING"):
                if opts.qp_solver_cond_ric_alg != 0:
                    print("Warning: Parametric sensitivities with condensing should be used with qp_solver_cond_ric_alg=0, as otherwise the full space Hessian needs to be factorized and the algorithm cannot handle indefinite ones.")

        if opts.with_value_sens_wrt_params:
            if dims.np_global == 0:
                raise ValueError('with_value_sens_wrt_params is only compatible if global parameters `p_global` are provided. Sensitivities wrt parameters have been refactored to use p_global instead of p in https://github.com/acados/acados/pull/1316. Got emty p_global.')
            if any([cost_type not in ["EXTERNAL", "LINEAR_LS"] for cost_type in cost_types_to_check]):
                raise ValueError('with_value_sens_wrt_params is only compatible with EXTERNAL cost_type' + suffix)
            if opts.N_horizon > 0 and opts.integrator_type != "DISCRETE":
                raise NotImplementedError('with_value_sens_wrt_params is only compatible with DISCRETE dynamics.')
            for horizon_type, constraint in bgp_type_constraint_pairs:
                if constraint is not None and any(ca.which_depends(constraint, model.p_global)):
                    raise NotImplementedError(f"with_value_sens_wrt_params is not supported for BGP constraints that depend on p_global. Got dependency on p_global for {horizon_type} constraint.")

        if opts.tau_min > 0 and "HPIPM" not in opts.qp_solver:
            raise ValueError('tau_min > 0 is only compatible with HPIPM.')

        if opts.qp_solver_cond_N is None:
            opts.qp_solver_cond_N = opts.N_horizon
        if opts.qp_solver_cond_N > opts.N_horizon:
            raise ValueError("qp_solver_cond_N > N_horizon is not supported.")

        if opts.qp_solver_cond_block_size is not None:
            if sum(opts.qp_solver_cond_block_size) != opts.N_horizon:
                raise ValueError(f'sum(qp_solver_cond_block_size) = {sum(opts.qp_solver_cond_block_size)} != N = {opts.N_horizon}.')
            if len(opts.qp_solver_cond_block_size) != opts.qp_solver_cond_N+1:
                raise ValueError(f'qp_solver_cond_block_size = {opts.qp_solver_cond_block_size} should have length qp_solver_cond_N+1 = {opts.qp_solver_cond_N+1}.')

        if opts.nlp_solver_type == "DDP":
            if opts.N_horizon == 0:
                raise ValueError("DDP solver only supported for N_horizon > 0.")
            if opts.qp_solver != "PARTIAL_CONDENSING_HPIPM" or opts.qp_solver_cond_N != opts.N_horizon:
                raise ValueError(f'DDP solver only supported for PARTIAL_CONDENSING_HPIPM with qp_solver_cond_N == N, got qp solver {opts.qp_solver} and qp_solver_cond_N {opts.qp_solver_cond_N}, N {opts.N_horizon}.')
            if any([dims.nbu, dims.nbx, dims.ng, dims.nh, dims.nphi]):
                raise ValueError(f'DDP only supports initial state constraints, got path constraints. Dimensions: dims.nbu = {dims.nbu}, dims.nbx = {dims.nbx}, dims.ng = {dims.ng}, dims.nh = {dims.nh}, dims.nphi = {dims.nphi}')
            if any([dims.ng_e, dims.nphi_e, dims.nh_e]):
                raise ValueError('DDP only supports initial state constraints, got terminal constraints.')

        if opts.qpscaling_scale_constraints != "NO_CONSTRAINT_SCALING" or opts.qpscaling_scale_objective != "NO_OBJECTIVE_SCALING":
            if opts.nlp_solver_type == "SQP_RTI":
                raise NotImplementedError('qpscaling_scale_constraints and qpscaling_scale_objective not supported for SQP_RTI solver.')

        # Set default parameters for globalization
        ddp_with_merit_or_funnel = opts.globalization == 'FUNNEL_L1PEN_LINESEARCH' or (opts.nlp_solver_type == "DDP" and opts.globalization == 'MERIT_BACKTRACKING')
        if opts.globalization_alpha_min is None:
            if ddp_with_merit_or_funnel:
                opts.globalization_alpha_min = 1e-17
            else:
                opts.globalization_alpha_min = 0.05

        if opts.globalization_alpha_reduction is None:
            if ddp_with_merit_or_funnel:
                opts.globalization_alpha_reduction = 0.5
            else:
                opts.globalization_alpha_reduction = 0.7

        if opts.globalization_eps_sufficient_descent is None:
            if ddp_with_merit_or_funnel:
                opts.globalization_eps_sufficient_descent = 1e-6
            else:
                opts.globalization_eps_sufficient_descent = 1e-4

        if opts.eval_residual_at_max_iter is None:
            if ddp_with_merit_or_funnel:
                opts.eval_residual_at_max_iter = True
            else:
                opts.eval_residual_at_max_iter = False

        if opts.globalization_full_step_dual is None:
            if ddp_with_merit_or_funnel:
                opts.globalization_full_step_dual = 1
            else:
                opts.globalization_full_step_dual = 0

        # AS-RTI
        if opts.as_rti_level in [1, 2] and any([cost_type.endswith("LINEAR_LS") for cost_type in cost_types_to_check]):
            raise NotImplementedError('as_rti_level in [1, 2] not supported for LINEAR_LS and NONLINEAR_LS cost type.')

        # sanity check for Funnel globalization and SQP
        if opts.globalization == 'FUNNEL_L1PEN_LINESEARCH' and opts.nlp_solver_type not in ['SQP', 'SQP_WITH_FEASIBLE_QP']:
            raise NotImplementedError('FUNNEL_L1PEN_LINESEARCH only supports SQP.')

        # termination
        if opts.nlp_solver_tol_min_step_norm is None:
            if ddp_with_merit_or_funnel:
                opts.nlp_solver_tol_min_step_norm = 1e-12
            else:
                opts.nlp_solver_tol_min_step_norm = 0.0

        # zoRO
        if self.zoro_description is not None:
            if opts.N_horizon == 0:
                raise ValueError('zoRO only supported for N_horizon > 0.')
            if not isinstance(self.zoro_description, ZoroDescription):
                raise TypeError('zoro_description should be of type ZoroDescription or None')
            else:
                self.zoro_description = process_zoro_description(self.zoro_description)

        # nlp_solver_warm_start_first_qp_from_nlp
        if opts.nlp_solver_warm_start_first_qp_from_nlp and (opts.qp_solver != "PARTIAL_CONDENSING_HPIPM" or opts.qp_solver_cond_N != opts.N_horizon):
            raise NotImplementedError('nlp_solver_warm_start_first_qp_from_nlp only supported for PARTIAL_CONDENSING_HPIPM with qp_solver_cond_N == N.')

        # Anderson acceleration
        if opts.with_anderson_acceleration:
            if opts.nlp_solver_type == "DDP":
                raise NotImplementedError('Anderson acceleration not supported for DDP solver.')
            if opts.globalization != "FIXED_STEP":
                raise NotImplementedError('Anderson acceleration only supported for FIXED_STEP globalization for now.')

        # check terminal stage
        for field in ('cost_expr_ext_cost_e', 'cost_expr_ext_cost_custom_hess_e',
                      'cost_y_expr_e', 'cost_psi_expr_e', 'cost_conl_custom_outer_hess_e',
                      'con_h_expr_e', 'con_phi_expr_e', 'con_r_expr_e',):
            val = getattr(model, field)
            if not is_empty(val) and (ca.depends_on(val, model.u) or ca.depends_on(val, model.z)):
                raise ValueError(f'{field} can not depend on u or z.')

        return


    def _get_external_function_header_templates(self, ) -> list:
        dims = self.dims
        name = self.model.name
        opts = self.solver_options
        template_list = []

        # dynamics
        if opts.N_horizon > 0:
            model_dir = os.path.join(self.code_export_directory, f'{name}_model')
            template_list.append(('model.in.h', f'{name}_model.h', model_dir))
        # constraints
        if any(np.array([dims.nh, dims.nh_e, dims.nh_0, dims.nphi, dims.nphi_e, dims.nphi_0]) > 0):
            constraints_dir = os.path.join(self.code_export_directory, f'{name}_constraints')
            template_list.append(('constraints.in.h', f'{name}_constraints.h', constraints_dir))
        # cost
        if any([self.cost.cost_type != 'LINEAR_LS', self.cost.cost_type_0 != 'LINEAR_LS', self.cost.cost_type_e != 'LINEAR_LS']):
            cost_dir = os.path.join(self.code_export_directory, f'{name}_cost')
            template_list.append(('cost.in.h', f'{name}_cost.h', cost_dir))

        return template_list


    def __get_template_list(self, cmake_builder=None) -> list:
        """
        returns a list of tuples in the form:
        (input_filename, output_filname)
        or
        (input_filename, output_filname, output_directory)
        """
        name = self.model.name
        opts = self.solver_options
        template_list = []

        template_list.append(('main.in.c', f'main_{name}.c'))
        template_list.append(('acados_solver.in.c', f'acados_solver_{name}.c'))
        template_list.append(('acados_solver.in.h', f'acados_solver_{name}.h'))
        template_list.append(('acados_solver.in.pxd', f'acados_solver.pxd'))
        if cmake_builder is not None:
            template_list.append(('CMakeLists.in.txt', 'CMakeLists.txt'))
        else:
            template_list.append(('Makefile.in', 'Makefile'))

        # sim
        if opts.N_horizon > 0 and self.solver_options.integrator_type != 'DISCRETE':
            template_list.append(('acados_sim_solver.in.c', f'acados_sim_solver_{name}.c'))
            template_list.append(('acados_sim_solver.in.h', f'acados_sim_solver_{name}.h'))
            template_list.append(('main_sim.in.c', f'main_sim_{name}.c'))

        # model
        template_list += self._get_external_function_header_templates()

        if self.dims.n_global_data > 0:
            template_list.append(('p_global_precompute_fun.in.h', f'{self.name}_p_global_precompute_fun.h'))

        # Simulink
        if self.simulink_opts is not None:
            template_list += self._get_matlab_simulink_template_list(name)
            template_list += self._get_integrator_simulink_template_list(name)

        return template_list


    @classmethod
    def _get_matlab_simulink_template_list(cls, name: str) -> list:
        template_list = []
        template_file = os.path.join('matlab_templates', 'acados_solver_sfun.in.c')
        template_list.append((template_file, f'acados_solver_sfunction_{name}.c'))
        template_file = os.path.join('matlab_templates', 'make_sfun.in.m')
        template_list.append((template_file, f'make_sfun_{name}.m'))
        # MEX wrapper files
        template_file = os.path.join('matlab_templates', 'mex_solver.in.m')
        template_list.append((template_file, f'{name}_mex_solver.m'))
        template_file = os.path.join('matlab_templates', 'make_mex.in.m')
        template_list.append((template_file, f'make_mex_{name}.m'))
        template_file = os.path.join('matlab_templates', 'acados_mex_create.in.c')
        template_list.append((template_file, f'acados_mex_create_{name}.c'))
        template_file = os.path.join('matlab_templates', 'acados_mex_free.in.c')
        template_list.append((template_file, f'acados_mex_free_{name}.c'))
        template_file = os.path.join('matlab_templates', 'acados_mex_solve.in.c')
        template_list.append((template_file, f'acados_mex_solve_{name}.c'))
        template_file = os.path.join('matlab_templates', 'acados_mex_set.in.c')
        template_list.append((template_file, f'acados_mex_set_{name}.c'))
        return template_list

    # dont render sim sfunctions for MOCP
    @classmethod
    def _get_integrator_simulink_template_list(cls, name: str) -> list:
        template_list = []
        template_file = os.path.join('matlab_templates', 'acados_sim_solver_sfun.in.c')
        template_list.append((template_file, f'acados_sim_solver_sfunction_{name}.c'))
        template_file = os.path.join('matlab_templates', 'make_sfun_sim.in.m')
        template_list.append((template_file, f'make_sfun_sim_{name}.m'))
        return template_list

    def render_templates(self, cmake_builder=None):

        # check json file
        json_path = os.path.abspath(self.json_file)
        if not os.path.exists(json_path):
            raise FileNotFoundError(f'Path "{json_path}" not found!')

        template_list = self.__get_template_list(cmake_builder=cmake_builder)

        # Render templates
        for tup in template_list:
            output_dir = self.code_export_directory if len(tup) <= 2 else tup[2]
            render_template(tup[0], tup[1], output_dir, json_path)

        # Custom templates
        acados_template_path = os.path.dirname(os.path.abspath(__file__))
        custom_template_glob = os.path.join(acados_template_path, 'custom_update_templates', '*')
        for tup in self.solver_options.custom_templates:
            render_template(tup[0], tup[1], self.code_export_directory, json_path, template_glob=custom_template_glob)
        return


    def dump_to_json(self) -> None:
        with open(self.json_file, 'w') as f:
            json.dump(self.to_dict(), f, default=make_object_json_dumpable, indent=4, sort_keys=True)
        return

    def generate_external_functions(self, context: Optional[GenerateContext] = None) -> GenerateContext:

        if context is None:
            # options for code generation
            code_gen_opts = AcadosCodegenOptions(
                ext_fun_expand_constr = self.solver_options.ext_fun_expand_constr,
                ext_fun_expand_cost = self.solver_options.ext_fun_expand_cost,
                ext_fun_expand_precompute = self.solver_options.ext_fun_expand_precompute,
                ext_fun_expand_dyn = self.solver_options.ext_fun_expand_dyn,
                code_export_directory = self.code_export_directory,
                with_solution_sens_wrt_params = self.solver_options.with_solution_sens_wrt_params,
                with_value_sens_wrt_params = self.solver_options.with_value_sens_wrt_params,
                generate_hess = self.solver_options.hessian_approx == 'EXACT',
            )

            context = GenerateContext(self.model.p_global, self.name, code_gen_opts)

        context = self._setup_code_generation_context(context)
        context.finalize()
        self.__external_function_files_model = context.get_external_function_file_list(ocp_specific=False)
        self.__external_function_files_ocp = context.get_external_function_file_list(ocp_specific=True)
        self.dims.n_global_data = context.get_n_global_data()

        return context


    def _setup_code_generation_context(self, context: GenerateContext, ignore_initial: bool = False, ignore_terminal: bool = False) -> GenerateContext:

        model = self.model
        constraints = self.constraints
        opts = self.solver_options

        check_casadi_version()
        self._setup_code_generation_context_dynamics(context)

        if opts.N_horizon > 0:
            if ignore_initial and ignore_terminal:
                stage_type_indices = [1]
            elif ignore_initial:
                stage_type_indices = [1, 2]
            elif ignore_terminal:
                stage_type_indices = [0, 1]
            else:
                stage_type_indices = [0, 1, 2]
        else:
            stage_type_indices = [2]

        stage_types = [val for i, val in enumerate(['initial', 'path', 'terminal']) if i in stage_type_indices]
        nhs = [val for i, val in enumerate(['nh_0', 'nh', 'nh_e']) if i in stage_type_indices]
        nphis = [val for i, val in enumerate(['nphi_0', 'nphi', 'nphi_e']) if i in stage_type_indices]
        cost_types = [val for i, val in enumerate(['cost_type_0', 'cost_type', 'cost_type_e']) if i in stage_type_indices]

        for attr_nh, attr_nphi, stage_type in zip(nhs, nphis, stage_types):
            if getattr(self.dims, attr_nh) > 0 or getattr(self.dims, attr_nphi) > 0:
                generate_c_code_constraint(context, model, constraints, stage_type)

        for attr, stage_type in zip(cost_types, stage_types):
            if getattr(self.cost, attr) == 'NONLINEAR_LS':
                generate_c_code_nls_cost(context, model, stage_type)
            elif getattr(self.cost, attr) == 'CONVEX_OVER_NONLINEAR':
                generate_c_code_conl_cost(context, model, stage_type)
            elif getattr(self.cost, attr) == 'EXTERNAL':
                generate_c_code_external_cost(context, model, stage_type)
            # TODO: generic

        return context


    def _setup_code_generation_context_dynamics(self, context: GenerateContext):
        opts = self.solver_options
        model = self.model

        if opts.N_horizon == 0:
            return

        code_gen_opts = context.opts

        # create code_export_dir, model_dir
        model_dir = os.path.join(code_gen_opts.code_export_directory, model.name + '_model')
        if not os.path.exists(model_dir):
            os.makedirs(model_dir)

        if self.model.dyn_ext_fun_type == 'casadi':
            if self.solver_options.integrator_type == 'ERK':
                generate_c_code_explicit_ode(context, model, model_dir)
            elif self.solver_options.integrator_type == 'IRK':
                generate_c_code_implicit_ode(context, model, model_dir)
            elif self.solver_options.integrator_type == 'LIFTED_IRK':
                if model.t != []:
                    raise NotImplementedError("LIFTED_IRK with time-varying dynamics not implemented yet.")
                generate_c_code_implicit_ode(context, model, model_dir)
            elif self.solver_options.integrator_type == 'GNSF':
                generate_c_code_gnsf(context, model, model_dir)
            elif self.solver_options.integrator_type == 'DISCRETE':
                generate_c_code_discrete_dynamics(context, model, model_dir)
            else:
                raise ValueError("ocp_generate_external_functions: unknown integrator type.")
        else:
            target_dir = os.path.join(code_gen_opts.code_export_directory, model_dir)
            target_location = os.path.join(target_dir, model.dyn_generic_source)
            shutil.copyfile(model.dyn_generic_source, target_location)
            context.add_external_function_file(model.dyn_generic_source, target_dir)


    def remove_x0_elimination(self) -> None:
        """Remove the elimination of x0 from the constraints, bounds on x0 are handled as general bounds on x."""
        self.constraints.remove_x0_elimination()

    def to_dict(self) -> dict:
        # Copy ocp object dictionary
        ocp_dict = dict(deepcopy(self).__dict__)

        # convert acados classes to dicts
        for key, v in ocp_dict.items():
            if isinstance(v, (AcadosModel, AcadosOcpDims, AcadosOcpConstraints, AcadosOcpCost, AcadosOcpOptions, ZoroDescription)):
                ocp_dict[key]=dict(getattr(self, key).__dict__)

        ocp_dict = format_class_dict(ocp_dict)
        return ocp_dict


    def copy_path_cost_to_stage_0(self):
        """Set all cost definitions at stage 0 to the corresponding path cost definitions."""
        cost = self.cost
        model = self.model

        cost.cost_type_0 = cost.cost_type
        cost.W_0 = cost.W
        cost.Vx_0 = cost.Vx
        cost.Vu_0 = cost.Vu
        cost.Vz_0 = cost.Vz
        cost.yref_0 = cost.yref
        cost.cost_ext_fun_type_0 = cost.cost_ext_fun_type

        model.cost_y_expr_0 = model.cost_y_expr
        model.cost_expr_ext_cost_0 = model.cost_expr_ext_cost
        model.cost_expr_ext_cost_custom_hess_0 = model.cost_expr_ext_cost_custom_hess
        model.cost_psi_expr_0 = model.cost_psi_expr
        model.cost_r_in_psi_expr_0 = model.cost_r_in_psi_expr
        return

    def translate_nls_cost_to_conl(self):
        """
        Translates a NONLINEAR_LS cost to a CONVEX_OVER_NONLINEAR cost.
        """
        casadi_symbol = self.model.get_casadi_symbol()
        # initial cost
        if self.cost.cost_type_0 is None:
            print("Initial cost is None, skipping.")
        elif self.cost.cost_type_0 == "CONVEX_OVER_NONLINEAR":
            print("Initial cost is already CONVEX_OVER_NONLINEAR, skipping.")
        elif self.cost.cost_type_0 == "NONLINEAR_LS":
            print("Translating initial NONLINEAR_LS cost to CONVEX_OVER_NONLINEAR.")
            self.cost.cost_type_0 = "CONVEX_OVER_NONLINEAR"
            ny_0 = self.model.cost_y_expr_0.shape[0]
            conl_res_0 = casadi_symbol('residual_conl', ny_0)
            self.model.cost_r_in_psi_expr_0 = conl_res_0
            self.model.cost_psi_expr_0 = .5 * conl_res_0.T @ ca.sparsify(ca.DM(self.cost.W_0)) @ conl_res_0
        else:
            raise TypeError(f"Terminal cost type must be NONLINEAR_LS, got cost_type_0 {self.cost.cost_type_0}.")

        # path cost
        if self.cost.cost_type == "CONVEX_OVER_NONLINEAR":
            print("Path cost is already CONVEX_OVER_NONLINEAR, skipping.")
        elif self.cost.cost_type == "NONLINEAR_LS":
            print("Translating path NONLINEAR_LS cost to CONVEX_OVER_NONLINEAR.")
            self.cost.cost_type = "CONVEX_OVER_NONLINEAR"
            ny = self.model.cost_y_expr.shape[0]
            conl_res = casadi_symbol('residual_conl', ny)
            self.model.cost_r_in_psi_expr = conl_res
            self.model.cost_psi_expr = .5 * conl_res.T @ ca.sparsify(ca.DM(self.cost.W)) @ conl_res
        else:
            raise TypeError(f"Path cost type must be NONLINEAR_LS, got cost_type {self.cost.cost_type}.")

        # terminal cost
        if self.cost.cost_type_e == "CONVEX_OVER_NONLINEAR":
            print("Terminal cost is already CONVEX_OVER_NONLINEAR, skipping.")
        elif self.cost.cost_type_e == "NONLINEAR_LS":
            print("Translating terminal NONLINEAR_LS cost to CONVEX_OVER_NONLINEAR.")
            self.cost.cost_type_e = "CONVEX_OVER_NONLINEAR"
            ny_e = self.model.cost_y_expr_e.shape[0]
            conl_res_e = casadi_symbol('residual_conl', ny_e)
            self.model.cost_r_in_psi_expr_e = conl_res_e
            self.model.cost_psi_expr_e = .5 * conl_res_e.T @ ca.sparsify(ca.DM(self.cost.W_e)) @ conl_res_e
        else:
            raise ValueError(f"Initial cost type must be NONLINEAR_LS, got cost_type_e {self.cost.cost_type_e}.")
        return


    def translate_cost_to_external_cost(self,
                                        p: Optional[Union[ca.SX, ca.MX]] = None,
                                        p_values: Optional[np.ndarray] = None,
                                        p_global: Optional[Union[ca.SX, ca.MX]] = None,
                                        p_global_values: Optional[np.ndarray] = None,
                                        yref_0: Optional[Union[ca.SX, ca.MX]] = None,
                                        yref: Optional[Union[ca.SX, ca.MX]] = None,
                                        yref_e: Optional[Union[ca.SX, ca.MX]] = None,
                                        W_0: Optional[Union[ca.SX, ca.MX]] = None,
                                        W: Optional[Union[ca.SX, ca.MX]] = None,
                                        W_e: Optional[Union[ca.SX, ca.MX]] = None,
                                        cost_hessian: str = 'EXACT',
                                        ):
        """
        Translates cost to EXTERNAL cost and optionally provide parametrization of references and weighting matrices.

        :param p: Optional CasADi symbolics with additional stagewise parameters which are used to define yref_0, yref, yref_e, W_0, W, W_e. Will be appended to model.p.
        :param p_values: numpy array with the same shape as p providing initial parameter values.
        :param p_global: Optional CasADi symbolics with additional global parameters which are used to define yref_0, yref, yref_e, W_0, W, W_e. Will be appended to model.p_global.
        :param p_global_values: numpy array with the same shape as p_global providing initial global parameter values.
        :param W_0, W, W_e: Optional CasADi expressions which should be used instead of the numerical values provided by the cost module, shapes should be (ny_0, ny_0), (ny, ny), (ny_e, ny_e).
        :param yref_0, yref, yref_e: Optional CasADi expressions which should be used instead of the numerical values provided by the cost module, shapes should be (ny_0, 1), (ny, 1), (ny_e, 1).
        :param cost_hessian: 'EXACT' or 'GAUSS_NEWTON', determines how the cost hessian is computed.
        """
        if cost_hessian not in ['EXACT', 'GAUSS_NEWTON']:
            raise ValueError(f"Invalid cost_hessian {cost_hessian}, should be 'EXACT' or 'GAUSS_NEWTON'.")

        casadi_symbolics_type = type(self.model.x)

        # check p, p_values and append
        if p is not None:
            if p_values is None:
                raise ValueError("If p is not None, also p_values need to be provided.")
            if not (is_column(p) and is_column(p_values)):
                raise ValueError("p, p_values need to be column vectors.")
            if p.shape[0] != p_values.shape[0]:
                raise ValueError(f"Mismatching shapes regarding p, p_values: p has shape {p.shape}, p_values has shape {p_values.shape}.")
            if not isinstance(p, casadi_symbolics_type):
                raise TypeError(f"p has wrong type, got {type(p)}, expected {casadi_symbolics_type}.")

            self.model.p = ca.vertcat(self.model.p, p)
            self.parameter_values = np.concatenate((self.parameter_values, p_values))

        if p_global is not None:
            if p_global_values is None:
                raise ValueError("If p_global is not None, also p_global_values need to be provided.")
            if not (is_column(p_global) and is_column(p_global_values)):
                raise ValueError("p_global, p_global_values need to be column vectors.")
            if p_global.shape[0] != p_global_values.shape[0]:
                raise ValueError(f"Mismatching shapes regarding p_global, p_global_values: p_global has shape {p_global.shape}, p_global_values has shape {p_global_values.shape}.")
            if not isinstance(p_global, casadi_symbolics_type):
                raise TypeError(f"p_global has wrong type, got {type(p_global)}, expected {casadi_symbolics_type}.")

            self.model.p_global = ca.vertcat(self.model.p_global, p_global)
            self.p_global_values = np.concatenate((self.p_global_values, p_global_values))

        self.translate_initial_cost_term_to_external(yref_0, W_0, cost_hessian)
        self.translate_intermediate_cost_term_to_external(yref, W, cost_hessian)
        self.translate_terminal_cost_term_to_external(yref_e, W_e, cost_hessian)


    def translate_initial_cost_term_to_external(self, yref_0: Optional[Union[ca.SX, ca.MX]] = None, W_0: Optional[Union[ca.SX, ca.MX]] = None, cost_hessian: str = 'EXACT'):

        if cost_hessian not in ['EXACT', 'GAUSS_NEWTON']:
            raise Exception(f"Invalid cost_hessian {cost_hessian}, should be 'EXACT' or 'GAUSS_NEWTON'.")

        if cost_hessian == 'GAUSS_NEWTON':
            if self.cost.cost_type_0 not in ['LINEAR_LS', 'NONLINEAR_LS', None]:
                raise Exception(f"cost_hessian 'GAUSS_NEWTON' is only supported for LINEAR_LS, NONLINEAR_LS cost types, got cost_type_0 = {self.cost.cost_type_0}.")

        casadi_symbolics_type = type(self.model.x)

        if yref_0 is None:
            yref_0 = self.cost.yref_0
        else:
            if yref_0.shape[0] != self.cost.yref_0.shape[0]:
                raise ValueError(f"yref_0 has wrong shape, got {yref_0.shape}, expected {self.cost.yref_0.shape}.")

            if not isinstance(yref_0, casadi_symbolics_type):
                raise TypeError(f"yref_0 has wrong type, got {type(yref_0)}, expected {casadi_symbolics_type}.")

        if W_0 is None:
            W_0 = self.cost.W_0
        else:
            if W_0.shape != self.cost.W_0.shape:
                raise ValueError(f"W_0 has wrong shape, got {W_0.shape}, expected {self.cost.W_0.shape}.")

            if not isinstance(W_0, casadi_symbolics_type):
                raise TypeError(f"W_0 has wrong type, got {type(W_0)}, expected {casadi_symbolics_type}.")

        if self.cost.cost_type_0 == "LINEAR_LS":
            self.model.cost_expr_ext_cost_0 = \
                self.__translate_ls_cost_to_external_cost(self.model.x, self.model.u, self.model.z,
                                                          self.cost.Vx_0, self.cost.Vu_0, self.cost.Vz_0,
                                                          yref_0, W_0)
        elif self.cost.cost_type_0 == "NONLINEAR_LS":
            self.model.cost_expr_ext_cost_0 = \
                self.__translate_nls_cost_to_external_cost(self.model.cost_y_expr_0, yref_0, W_0)

            if cost_hessian == 'GAUSS_NEWTON':
                self.model.cost_expr_ext_cost_custom_hess_0 = self.__get_gn_hessian_expression_from_nls_cost(self.model.cost_y_expr_0, yref_0, W_0, self.model.x, self.model.u, self.model.z)

        elif self.cost.cost_type_0 == "CONVEX_OVER_NONLINEAR":
            self.model.cost_expr_ext_cost_0 = \
                self.__translate_conl_cost_to_external_cost(self.model.cost_r_in_psi_expr_0, self.model.cost_psi_expr_0,
                                                            self.model.cost_y_expr_0, yref_0)

        if self.cost.cost_type_0 is not None:
            self.cost.cost_type_0 = 'EXTERNAL'


    def translate_intermediate_cost_term_to_external(self, yref: Optional[Union[ca.SX, ca.MX]] = None, W: Optional[Union[ca.SX, ca.MX]] = None, cost_hessian: str = 'EXACT'):

        if cost_hessian not in ['EXACT', 'GAUSS_NEWTON']:
            raise ValueError(f"Invalid cost_hessian {cost_hessian}, should be 'EXACT' or 'GAUSS_NEWTON'.")

        if cost_hessian == 'GAUSS_NEWTON':
            if self.cost.cost_type not in ['LINEAR_LS', 'NONLINEAR_LS']:
                raise ValueError(f"cost_hessian 'GAUSS_NEWTON' is only supported for LINEAR_LS, NONLINEAR_LS cost types, got cost_type = {self.cost.cost_type}.")

        casadi_symbolics_type = type(self.model.x)

        if yref is None:
            yref = self.cost.yref
        else:
            if yref.shape[0] != self.cost.yref.shape[0]:
                raise ValueError(f"yref has wrong shape, got {yref.shape}, expected {self.cost.yref.shape}.")

            if not isinstance(yref, casadi_symbolics_type):
                raise TypeError(f"yref has wrong type, got {type(yref)}, expected {casadi_symbolics_type}.")

        if W is None:
            W = self.cost.W
        else:
            if W.shape != self.cost.W.shape:
                raise ValueError(f"W has wrong shape, got {W.shape}, expected {self.cost.W.shape}.")

            if not isinstance(W, casadi_symbolics_type):
                raise TypeError(f"W has wrong type, got {type(W)}, expected {casadi_symbolics_type}.")

        if self.cost.cost_type == "LINEAR_LS":
            self.model.cost_expr_ext_cost = \
                self.__translate_ls_cost_to_external_cost(self.model.x, self.model.u, self.model.z,
                                                          self.cost.Vx, self.cost.Vu, self.cost.Vz,
                                                          yref, W)
        elif self.cost.cost_type == "NONLINEAR_LS":
            self.model.cost_expr_ext_cost = \
                self.__translate_nls_cost_to_external_cost(self.model.cost_y_expr, yref, W)
            if cost_hessian == 'GAUSS_NEWTON':
                self.model.cost_expr_ext_cost_custom_hess = self.__get_gn_hessian_expression_from_nls_cost(self.model.cost_y_expr, yref, W, self.model.x, self.model.u, self.model.z)

        elif self.cost.cost_type == "CONVEX_OVER_NONLINEAR":
            self.model.cost_expr_ext_cost = \
                self.__translate_conl_cost_to_external_cost(self.model.cost_r_in_psi_expr, self.model.cost_psi_expr,
                                                            self.model.cost_y_expr, yref)

        self.cost.cost_type = 'EXTERNAL'


    def translate_terminal_cost_term_to_external(self, yref_e: Optional[Union[ca.SX, ca.MX]] = None, W_e: Optional[Union[ca.SX, ca.MX]] = None, cost_hessian: str = 'EXACT'):
        if cost_hessian not in ['EXACT', 'GAUSS_NEWTON']:
            raise ValueError(f"Invalid cost_hessian {cost_hessian}, should be 'EXACT' or 'GAUSS_NEWTON'.")

        if cost_hessian == 'GAUSS_NEWTON':
            if self.cost.cost_type_e not in ['LINEAR_LS', 'NONLINEAR_LS']:
                raise ValueError(f"cost_hessian 'GAUSS_NEWTON' is only supported for LINEAR_LS, NONLINEAR_LS cost types, got cost_type_e = {self.cost.cost_type_e}.")

        casadi_symbolics_type = type(self.model.x)

        if yref_e is None:
            yref_e = self.cost.yref_e
        else:
            if yref_e.shape[0] != self.cost.yref_e.shape[0]:
                raise ValueError(f"yref_e has wrong shape, got {yref_e.shape}, expected {self.cost.yref_e.shape}.")

            if not isinstance(yref_e, casadi_symbolics_type):
                raise TypeError(f"yref_e has wrong type, got {type(yref_e)}, expected {casadi_symbolics_type}.")

        if W_e is None:
            W_e = self.cost.W_e
        else:
            if W_e.shape != self.cost.W_e.shape:
                raise ValueError(f"W_e has wrong shape, got {W_e.shape}, expected {self.cost.W_e.shape}.")

            if not isinstance(W_e, casadi_symbolics_type):
                raise TypeError(f"W_e has wrong type, got {type(W_e)}, expected {casadi_symbolics_type}.")

        if self.cost.cost_type_e == "LINEAR_LS":
            self.model.cost_expr_ext_cost_e = \
                self.__translate_ls_cost_to_external_cost(self.model.x, self.model.u, self.model.z,
                                                          self.cost.Vx_e, None, None,
                                                          yref_e, W_e)
        elif self.cost.cost_type_e == "NONLINEAR_LS":
            self.model.cost_expr_ext_cost_e = \
                self.__translate_nls_cost_to_external_cost(self.model.cost_y_expr_e, yref_e, W_e)
            if cost_hessian == 'GAUSS_NEWTON':
                self.model.cost_expr_ext_cost_custom_hess_e = self.__get_gn_hessian_expression_from_nls_cost(self.model.cost_y_expr_e, yref_e, W_e, self.model.x, [], self.model.z)

        elif self.cost.cost_type_e == "CONVEX_OVER_NONLINEAR":
            self.model.cost_expr_ext_cost_e = \
                self.__translate_conl_cost_to_external_cost(self.model.cost_r_in_psi_expr_e, self.model.cost_psi_expr_e,
                                                            self.model.cost_y_expr_e, yref_e)
        self.cost.cost_type_e = 'EXTERNAL'


    @staticmethod
    def __translate_ls_cost_to_external_cost(x, u, z, Vx, Vu, Vz, yref, W):
        res = 0
        if not is_empty(Vx):
            res += Vx @ x
        if not is_empty(Vu):
            res += Vu @ u
        if not is_empty(Vz):
            res += Vz @ z
        res -= yref

        return 0.5 * (res.T @ W @ res)

    @staticmethod
    def __translate_nls_cost_to_external_cost(y_expr, yref, W):
        res = y_expr - yref
        return 0.5 * (res.T @ W @ res)

    @staticmethod
    def __get_gn_hessian_expression_from_nls_cost(y_expr, yref, W, x, u, z):
        res = y_expr - yref
        ux = ca.vertcat(u, x)
        inner_jac = ca.jacobian(res, ux)
        gn_hess = inner_jac.T @ W @ inner_jac
        return gn_hess

    @staticmethod
    def __translate_conl_cost_to_external_cost(r, psi, y_expr, yref):
        return ca.substitute(psi, r, y_expr - yref)

    def formulate_constraint_as_L2_penalty(
        self,
        constr_expr: ca.SX,
        weight: float,
        upper_bound: Optional[float],
        lower_bound: Optional[float],
        residual_name: str = "new_residual",
        constraint_type: str = "path",
    ) -> None:
        """
        Formulate a constraint as an L2 penalty and add it to the current cost.
        """

        casadi_symbol = self.model.get_casadi_symbol()

        if upper_bound is None and lower_bound is None:
            raise ValueError("Either upper or lower bound must be provided.")

        # compute violation expression
        violation_expr = 0.0
        y_ref_new = np.zeros(1)
        if upper_bound is not None:
            violation_expr = ca.fmax(violation_expr, (constr_expr - upper_bound))
        if lower_bound is not None:
            violation_expr = ca.fmax(violation_expr, (lower_bound - constr_expr))

        # add penalty as cost
        if constraint_type == "path":
            self.cost.yref = np.concatenate((self.cost.yref, y_ref_new))
            self.model.cost_y_expr = ca.vertcat(self.model.cost_y_expr, violation_expr)
            if self.cost.cost_type == "NONLINEAR_LS":
                self.cost.W = block_diag(self.cost.W, weight)
            elif self.cost.cost_type == "CONVEX_OVER_NONLINEAR":
                new_residual = casadi_symbol(residual_name, constr_expr.shape)
                self.model.cost_r_in_psi_expr = ca.vertcat(self.model.cost_r_in_psi_expr, new_residual)
                self.model.cost_psi_expr += .5 * weight * new_residual**2
            elif self.cost.cost_type == "EXTERNAL":
                self.model.cost_expr_ext_cost += .5 * weight * violation_expr**2
            else:
                raise NotImplementedError(f"formulate_constraint_as_L2_penalty not implemented for path cost with cost_type {self.cost.cost_type}.")
        elif constraint_type == "initial":
            self.cost.yref_0 = np.concatenate((self.cost.yref_0, y_ref_new))
            self.model.cost_y_expr_0 = ca.vertcat(self.model.cost_y_expr_0, violation_expr)
            if self.cost.cost_type_0 == "NONLINEAR_LS":
                self.cost.W_0 = block_diag(self.cost.W_0, weight)
            elif self.cost.cost_type_0 == "CONVEX_OVER_NONLINEAR":
                new_residual = casadi_symbol(residual_name, constr_expr.shape)
                self.model.cost_r_in_psi_expr_0 = ca.vertcat(self.model.cost_r_in_psi_expr_0, new_residual)
                self.model.cost_psi_expr_0 += .5 * weight * new_residual**2
            elif self.cost.cost_type_0 == "EXTERNAL":
                self.model.cost_expr_ext_cost_0 += .5 * weight * violation_expr**2
            else:
                raise NotImplementedError(f"formulate_constraint_as_L2_penalty not implemented for initial cost with cost_type_0 {self.cost.cost_type_0}.")
        elif constraint_type == "terminal":
            self.cost.yref_e = np.concatenate((self.cost.yref_e, y_ref_new))
            self.model.cost_y_expr_e = ca.vertcat(self.model.cost_y_expr_e, violation_expr)
            if self.cost.cost_type_e == "NONLINEAR_LS":
                self.cost.W_e = block_diag(self.cost.W_e, weight)
            elif self.cost.cost_type_e == "CONVEX_OVER_NONLINEAR":
                new_residual = casadi_symbol(residual_name, constr_expr.shape)
                self.model.cost_r_in_psi_expr_e = ca.vertcat(self.model.cost_r_in_psi_expr_e, new_residual)
                self.model.cost_psi_expr_e += .5 * weight * new_residual**2
            elif self.cost.cost_type_e == "EXTERNAL":
                self.model.cost_expr_ext_cost_e += .5 * weight * violation_expr**2
            else:
                raise NotImplementedError(f"formulate_constraint_as_L2_penalty not implemented for terminal cost with cost_type_e {self.cost.cost_type_e}.")
        return


    def formulate_constraint_as_Huber_penalty(
        self,
        constr_expr: Union[ca.SX, ca.MX],
        weight: float,
        upper_bound: Optional[float]=None,
        lower_bound: Optional[float]=None,
        residual_name: str = "new_residual",
        huber_delta: float = 1.0,
        use_xgn = True,
        min_hess = 0,
    ) -> None:
        """
        Formulate a constraint as Huber penalty and add it to the current cost.

        use_xgn: if true an XGN Hessian is used, if false a GGN Hessian (= exact Hessian, in this case) is used.
        min_hess: provide a minimum value for the hessian
        weight: weight of the penalty corresponding to Hessian in quadratic region
        """
        if isinstance(constr_expr, ca.MX):
            casadi_symbol = ca.MX.sym
            casadi_zeros = ca.MX.zeros
        elif isinstance(constr_expr, ca.SX):
            casadi_symbol = ca.SX.sym
            casadi_zeros = ca.SX.zeros

        # if (upper_bound is None or lower_bound is None):
        #     raise NotImplementedError("only symmetric Huber for now")
        if upper_bound is None and lower_bound is None:
            raise ValueError("Either upper or lower bound must be provided.")

        if self.cost.cost_type != "CONVEX_OVER_NONLINEAR":
            raise ValueError("Huber penalty is only supported for CONVEX_OVER_NONLINEAR cost type.")

        if use_xgn and is_empty(self.model.cost_conl_custom_outer_hess):
            # switch to XGN Hessian start with exact Hessian of previously defined cost
            exact_cost_hess = ca.hessian(self.model.cost_psi_expr, self.model.cost_r_in_psi_expr)[0]
            self.model.cost_conl_custom_outer_hess = exact_cost_hess

        # define residual
        new_residual = casadi_symbol(residual_name, constr_expr.shape)

        if upper_bound is not None and lower_bound is not None:
            if upper_bound < lower_bound:
                raise ValueError("Upper bound must be greater than lower bound.")
            # normalize constraint to [-1, 1]
            width = upper_bound - lower_bound
            center = lower_bound + 0.5 * width
            constr_expr = 2 * (constr_expr - center) / width

            # define penalty
            penalty, penalty_grad, penalty_hess, penalty_hess_xgn = \
                symmetric_huber_penalty(new_residual, delta=huber_delta, w=weight*width**2, min_hess=min_hess)
        elif upper_bound is not None:
            # define penalty
            constr_expr -= upper_bound
            penalty, penalty_grad, penalty_hess, penalty_hess_xgn = \
                one_sided_huber_penalty(new_residual, delta=huber_delta, w=weight, min_hess=min_hess)
        elif lower_bound is not None:
            raise NotImplementedError("lower bound only not implemented, please change sign on constraint and use upper bound.")

        # add penalty to cost
        self.model.cost_r_in_psi_expr = ca.vertcat(self.model.cost_r_in_psi_expr, new_residual)
        self.model.cost_psi_expr += penalty
        self.model.cost_y_expr = ca.vertcat(self.model.cost_y_expr, constr_expr)
        self.cost.yref = np.concatenate((self.cost.yref, np.zeros(1)))

        # add Hessian term
        if use_xgn:
            zero_offdiag = casadi_zeros(self.model.cost_conl_custom_outer_hess.shape[0], penalty_hess_xgn.shape[1])
            self.model.cost_conl_custom_outer_hess = ca.blockcat(self.model.cost_conl_custom_outer_hess,
                                                                zero_offdiag, zero_offdiag.T, penalty_hess_xgn)
        elif not is_empty(self.model.cost_conl_custom_outer_hess):
            zero_offdiag = casadi_zeros(self.model.cost_conl_custom_outer_hess.shape[0], penalty_hess_xgn.shape[1])
            # add penalty Hessian to existing Hessian
            self.model.cost_conl_custom_outer_hess = ca.blockcat(self.model.cost_conl_custom_outer_hess,
                                                                zero_offdiag, zero_offdiag.T, penalty_hess)

        return


    def add_linear_constraint(self, C: np.ndarray, D: np.ndarray, lg: np.ndarray, ug: np.ndarray) -> None:
        """
        Add a linear constraint of the form lg <= C * x + D * u <= ug to the OCP.
        """

        if C.shape[0] != lg.shape[0] or C.shape[0] != ug.shape[0]:
            raise ValueError("C, lg, ug must have the same number of rows.")

        if D.shape[0] != C.shape[0]:
            raise ValueError("C and D must have the same number of rows.")

        if self.constraints.C.shape[0] == 0:
            # no linear constraints have been added yet
            self.constraints.C = C
            self.constraints.D = D
            self.constraints.lg = lg
            self.constraints.ug = ug
        else:
            self.constraints.C = ca.vertcat(self.constraints.C, C)
            self.constraints.D = ca.vertcat(self.constraints.D, D)
            self.constraints.lg = ca.vertcat(self.constraints.lg, lg)
            self.constraints.ug = ca.vertcat(self.constraints.ug, ug)

        return

    def translate_to_feasibility_problem(self,
                                        keep_x0: bool=False,
                                        keep_cost: bool=False,
                                        parametric_x0: bool=False) -> None:
        """
        Translate an OCP to a feasibility problem by removing all cost term and then formulating all constraints as L2 penalties.

        Note: all weights are set to 1.0 for now.
        Options to specify weights should be implemented later for advanced use cases.

        :param keep_x0: if True, x0 constraint is kept as a constraint
        :param keep_cost: if True, cost is not removed before formulating constraints as penalties
        :param parametric_x0: if True, replace the value of the initial state constraint with a parameter that is appended to the model parameters.
        """

        self.model.make_consistent(self.dims) # sets the correct MX/SX defaults
        model = self.model
        cost = self.cost
        constraints = self.constraints
        new_constraints = AcadosOcpConstraints()

        if keep_cost:
            # initial stage - if not set, copy fields from path constraints
            if cost.cost_type_0 is None:
                self.copy_path_cost_to_stage_0()
        else:
            # set cost to zero
            cost.cost_type = "NONLINEAR_LS"
            cost.cost_type_e = "NONLINEAR_LS"
            cost.cost_type_0 = "NONLINEAR_LS"

            cost.yref = np.array([])
            cost.yref_0 = np.array([])
            cost.yref_e = np.array([])

            zeros = model.get_casadi_zeros()
            model.cost_y_expr = zeros((0, 0))
            model.cost_y_expr_e = zeros((0, 0))
            model.cost_y_expr_0 = zeros((0, 0))

            cost.W = np.zeros((0, 0))
            cost.W_e = np.zeros((0, 0))
            cost.W_0 = np.zeros((0, 0))

        expr_bound_list = [
            (model.x[constraints.idxbx], constraints.lbx, constraints.ubx),
            (model.u[constraints.idxbu], constraints.lbu, constraints.ubu),
            (model.con_h_expr, constraints.lh, constraints.uh),
        ]

        if casadi_length(model.con_phi_expr) > 0:
            phi_o_r_expr = ca.substitute(model.con_phi_expr, model.con_r_in_phi, model.con_r_expr)
            expr_bound_list.append((phi_o_r_expr, constraints.lphi, constraints.uphi))
            # NOTE: for now, we don't exploit convex over nonlinear structure of phi

        for constr_expr, lower_bound, upper_bound in expr_bound_list:
            for i in range(casadi_length(constr_expr)):
                self.formulate_constraint_as_L2_penalty(constr_expr[i], weight=1.0, upper_bound=upper_bound[i], lower_bound=lower_bound[i])

        model.con_h_expr = None
        model.con_phi_expr = None
        model.con_r_expr = None
        model.con_r_in_phi = None

        # formulate **terminal** constraints as L2 penalties
        expr_bound_list_e = [
            (model.x[constraints.idxbx_e], constraints.lbx_e, constraints.ubx_e),
            (model.con_h_expr_e, constraints.lh_e, constraints.uh_e),
        ]

        if casadi_length(model.con_phi_expr_e) > 0:
            phi_o_r_expr_e = ca.substitute(model.con_phi_expr_e, model.con_r_in_phi_e, model.con_r_expr_e)
            expr_bound_list_e.append((phi_o_r_expr_e, constraints.lphi_e, constraints.uphi_e))
            # NOTE: for now, we don't exploit convex over nonlinear structure of phi

        for constr_expr, lower_bound, upper_bound in expr_bound_list_e:
            for i in range(casadi_length(constr_expr)):
                self.formulate_constraint_as_L2_penalty(constr_expr[i], weight=1.0, upper_bound=upper_bound[i], lower_bound=lower_bound[i], constraint_type="terminal")

        model.con_h_expr_e = None
        model.con_phi_expr_e = None
        model.con_r_expr_e = None
        model.con_r_in_phi_e = None

        # Convert initial conditions to l2 penalty
        # Expressions for control constraints on u
        expr_bound_list_0 = [
            (model.u[constraints.idxbu], constraints.lbu, constraints.ubu),
            (model.con_h_expr_0, constraints.lh_0, constraints.uh_0),
        ]

        # initial state constraint
        if (keep_x0 or parametric_x0) and not constraints.has_x0:
            raise NotImplementedError("translate_to_feasibility_problem: options keep_x0, parametric_x0 not defined for problems without x0 constraints.")
        if parametric_x0 and keep_x0:
            raise NotImplementedError("translate_to_feasibility_problem: parametric_x0 and keep cannot both be True.")
        if keep_x0:
            new_constraints.x0 = constraints.lbx_0
        elif parametric_x0:
            symbol = model.get_casadi_symbol()
            param_x0 = symbol('param_x0', len(constraints.idxbx_0))
            new_params = constraints.lbx_0
            model.p = ca.vertcat(model.p, param_x0)
            self.parameter_values = np.concatenate((self.parameter_values, new_params))
            expr_bound_list_0.append((model.x[constraints.idxbx_0], param_x0, param_x0))
        else:
            expr_bound_list_0.append((model.x[constraints.idxbx_0], constraints.lbx_0, constraints.ubx_0))

        if casadi_length(model.con_phi_expr_0) > 0:
            phi_o_r_expr_0 = ca.substitute(model.con_phi_expr_0, model.con_r_in_phi_0, model.con_r_expr_0)
            expr_bound_list_0.append((phi_o_r_expr_0, constraints.lphi_0, constraints.uphi_0))
            # NOTE: for now, we don't exploit convex over nonlinear structure of phi

        for constr_expr, lower_bound, upper_bound in expr_bound_list_0:
            for i in range(casadi_length(constr_expr)):
                self.formulate_constraint_as_L2_penalty(constr_expr[i], weight=1.0, upper_bound=upper_bound[i], lower_bound=lower_bound[i], constraint_type="initial")

        model.con_h_expr_0 = None
        model.con_phi_expr_0 = None
        model.con_r_expr_0 = None
        model.con_r_in_phi_0 = None

        # delete constraint fromulation from constraints object
        self.constraints = new_constraints

    def augment_with_t0_param(self) -> None:
        """Add a parameter t0 to the model and set it to 0.0."""
        # TODO only needed in benchmark for problems with time-varying references.
        # maybe remove this function and model.t0 from acados (and move to benchmark)
        if self.model.t0 is not None:
            raise ValueError("Parameter t0 is already present in the model.")
        self.model.t0 = ca.SX.sym("t0")
        self.model.p = ca.vertcat(self.model.p, self.model.t0)
        self.parameter_values = np.append(self.parameter_values, [0.0])
        self.p_global_values = np.append(self.p_global_values, [0.0])
        return


    def detect_cost_type(self, model: AcadosModel, cost: AcadosOcpCost, dims: AcadosOcpDims, stage_type: str) -> None:
        """
        If the cost type of a stage (initial, path or terminal) is set to AUTO, try to reformulate it as a LINEAR_LS cost.
        If that is not possible (cost is not quadratic or includes parameters), use the EXTERNAL cost type.
        """
        # Extract model variables
        x = model.x
        u = model.u
        z = model.z
        p = model.p

        # Check type
        if not isinstance(x, ca.SX):
            raise ValueError("Cost type detection only works for casadi.SX!")

        nx = casadi_length(x)
        nu = casadi_length(u)
        nz = casadi_length(z)

        print('--------------------------------------------------------------')
        if stage_type == 'terminal':
            expr_cost = model.cost_expr_ext_cost_e
            print('Structure detection for terminal cost term')
        elif stage_type == 'path':
            expr_cost = model.cost_expr_ext_cost
            print('Structure detection for path cost')
        elif stage_type == 'initial':
            expr_cost = model.cost_expr_ext_cost_0
            print('Structure detection for initial cost term')

        if not (isinstance(expr_cost, ca.SX) or isinstance(expr_cost, ca.MX)):
            print('expr_cost =', expr_cost)
            raise ValueError("Cost type detection requires definition of cost term as CasADi SX or MX.")

        if ca.is_quadratic(expr_cost, x) and ca.is_quadratic(expr_cost, u) and ca.is_quadratic(expr_cost, z) \
                and not any(ca.which_depends(expr_cost, p)) and not any(ca.which_depends(expr_cost, model.p_global)) \
                and not any(ca.which_depends(expr_cost, model.t)):

            if expr_cost.is_zero():
                print('Cost function is zero -> Reformulating as LINEAR_LS cost.')
                ny = 0
                Vx, Vu, Vz, W, y_ref, y = [], [], [], [], [], []
            else:
                cost_fun = ca.Function('cost_fun', [x, u, z], [expr_cost])
                dummy = ca.SX.sym('dummy', 1, 1)

                print('Cost function is quadratic -> Reformulating as LINEAR_LS cost.')

                Hxuz_fun = ca.Function('Hxuz_fun', [dummy], [ca.hessian(expr_cost, ca.vertcat(x, u, z))[0]])
                H_xuz = np.array(Hxuz_fun(0))

                xuz_idx = []
                for i in range(nx + nu + nz):
                    if np.any(H_xuz[i, :]):
                        xuz_idx.append(i)

                x_idx = list(set(range(nx)) & set(xuz_idx))
                u_idx = list(set(range(nx, nx + nu)) & set(xuz_idx))
                z_idx = list(set(range(nx + nu, nx + nu + nz)) & set(xuz_idx))

                ny = len(xuz_idx)

                Vx = np.zeros((ny, nx))
                Vu = np.zeros((ny, nu))
                Vz = np.zeros((ny, nz))
                W = np.zeros((ny, ny))

                i = 0
                for idx in x_idx:
                    Vx[i, idx] = 1
                    W[i, :] = H_xuz[idx, xuz_idx] / 2
                    i += 1

                for idx in u_idx:
                    iu = idx - nx
                    Vu[i, iu] = 1
                    W[i, :] = H_xuz[idx, xuz_idx] / 2
                    i += 1

                for idx in z_idx:
                    iz = idx - nx - nu
                    Vz[i, iz] = 1
                    W[i, :] = H_xuz[idx, xuz_idx] / 2
                    i += 1

                xuz = ca.vertcat(x, u, z)
                y = xuz[xuz_idx]
                jac_fun = ca.Function('jac_fun', [y], [ca.jacobian(expr_cost, y).T])
                y_ref = np.linalg.solve(W, -0.5 * np.array(jac_fun(np.zeros((ny, 1)))))

                y = -y_ref + Vx @ x + Vu @ u
                if nz > 0:
                    y += Vz @ z

                lls_cost_fun = ca.Function('lls_cost_fun', [x, u, z], [ca.mtimes(y.T, ca.mtimes(W, y))])

                rel_err_tol = 1e-13
                for jj in range(5):
                    x0 = np.random.rand(nx)
                    u0 = np.random.rand(nu)
                    z0 = np.random.rand(nz)

                    val1 = np.array(lls_cost_fun(x0, u0, z0))
                    val2 = np.array(cost_fun(x0, u0, z0))
                    diff_eval = np.abs(val1 - val2)
                    rel_error = diff_eval / np.maximum(np.abs(val1), np.abs(val2))
                    if rel_error > rel_err_tol:
                        print(f'Something went wrong when reformulating with linear least square cost, '
                            f'got relative error {rel_error:.2e}, should be < {rel_err_tol:.2e}')
                        raise RuntimeError('Reformulation error')

                W = 2 * W

            # Extract output
            if stage_type == 'terminal':
                if np.any(Vu):
                    raise ValueError('Terminal cost term cannot depend on the control input (u)!')
                if np.any(Vz):
                    raise ValueError('Terminal cost term cannot depend on the algebraic variables (z)!')
                cost.cost_type_e = 'LINEAR_LS'
                dims.ny_e = ny
                cost.Vx_e = Vx
                cost.W_e = W
                cost.yref_e = y_ref
            elif stage_type == 'path':
                cost.cost_type = 'LINEAR_LS'
                dims.ny = ny
                cost.Vx = Vx
                cost.Vu = Vu
                cost.Vz = Vz
                cost.W = W
                cost.yref = y_ref
            elif stage_type == 'initial':
                cost.cost_type_0 = 'LINEAR_LS'
                dims.ny_0 = ny
                cost.Vx_0 = Vx
                cost.Vu_0 = Vu
                cost.Vz_0 = Vz
                cost.W_0 = W
                cost.yref_0 = y_ref

            print('\n\nReformulated cost term in linear least squares form with:')
            print('cost = 0.5 * || Vx * x + Vu * u + Vz * z - y_ref ||_W\n')
            print('Vx\n', Vx)
            print('Vu\n', Vu)
            print('Vz\n', Vz)
            print('W\n', W)
            print('y_ref\n', y_ref)
            print('y (symbolic)\n', y)
            print('NOTE: These numerical values can be updated online using the appropriate setters.')

        else:
            print('\n\nCost function is not quadratic or includes parameters -> Using external cost\n\n')
            if stage_type == 'terminal':
                cost.cost_type_e = 'EXTERNAL'
            elif stage_type == 'path':
                cost.cost_type = 'EXTERNAL'
            elif stage_type == 'initial':
                cost.cost_type_0 = 'EXTERNAL'

        print('--------------------------------------------------------------')

    def ensure_solution_sensitivities_available(self, parametric=True) -> None:
        """
        Check if the options are set correctly for calculating sensitivities.

        :param parametric: if True, check also if parametric sensitivities are available.

        :raises NotImplementedError: if the QP solver is not HPIPM.
        :raises ValueError: if the Hessian approximation or regularization method is not set correctly for parametric sensitivities.
        """
        has_custom_hess = self.model._has_custom_hess()

        self.solver_options._ensure_solution_sensitivities_available(
            parametric,
            has_custom_hess
        )

    def get_initial_cost_expression(self):
        model = self.model
        if self.cost.cost_type == "LINEAR_LS":
            if is_empty(self.cost.Vx_0):
                return 0

            y = self.cost.Vx_0 @ model.x + self.cost.Vu_0 @ model.u

            if not is_empty(self.cost.Vz_0):
                y += self.cost.Vz @ model.z
            residual = y - self.cost.yref_0
            cost_dot = 0.5 * (residual.T @ self.cost.W_0 @ residual)

        elif self.cost.cost_type == "NONLINEAR_LS":
            residual = model.cost_y_expr_0 - self.cost.yref_0
            cost_dot = 0.5 * (residual.T @ self.cost.W_0 @ residual)

        elif self.cost.cost_type == "EXTERNAL":
            cost_dot = model.cost_expr_ext_cost_0

        elif self.cost.cost_type == "CONVEX_OVER_NONLINEAR":
            cost_dot = ca.substitute(
            model.cost_psi_expr_0, model.cost_r_in_psi_expr_0, model.cost_y_expr_0)
        else:
            raise ValueError("create_model_with_cost_state: Unknown cost type.")

        return cost_dot


    def get_path_cost_expression(self):
        model = self.model
        if self.cost.cost_type == "LINEAR_LS":
            if is_empty(self.cost.Vx):
                return 0

            y = self.cost.Vx @ model.x + self.cost.Vu @ model.u

            if not is_empty(self.cost.Vz):
                y += self.cost.Vz @ model.z
            residual = y - self.cost.yref
            cost_dot = 0.5 * (residual.T @ self.cost.W @ residual)

        elif self.cost.cost_type == "NONLINEAR_LS":
            residual = model.cost_y_expr - self.cost.yref
            cost_dot = 0.5 * (residual.T @ self.cost.W @ residual)

        elif self.cost.cost_type == "EXTERNAL":
            cost_dot = model.cost_expr_ext_cost

        elif self.cost.cost_type == "CONVEX_OVER_NONLINEAR":
            cost_dot = ca.substitute(
            model.cost_psi_expr, model.cost_r_in_psi_expr, model.cost_y_expr)
        else:
            raise ValueError("create_model_with_cost_state: Unknown cost type.")

        return cost_dot


    def get_terminal_cost_expression(self):
        model = self.model
        if self.cost.cost_type_e == "LINEAR_LS":
            if is_empty(self.cost.Vx_e):
                return 0.0
            y = self.cost.Vx_e @ model.x
            residual = y - self.cost.yref_e
            cost_dot = 0.5 * (residual.T @ self.cost.W_e @ residual)

        elif self.cost.cost_type_e == "NONLINEAR_LS":
            residual = model.cost_y_expr_e - self.cost.yref_e
            cost_dot = 0.5 * (residual.T @ self.cost.W_e @ residual)

        elif self.cost.cost_type_e == "EXTERNAL":
            cost_dot = model.cost_expr_ext_cost_e

        elif self.cost.cost_type_e == "CONVEX_OVER_NONLINEAR":
            cost_dot = ca.substitute(
            model.cost_psi_expr_e, model.cost_r_in_psi_expr_e, model.cost_y_expr_e)
        else:
            raise ValueError(f"create_model_with_cost_state: Unknown terminal cost type {self.cost.cost_type_e}.")

        return cost_dot


    def create_default_initial_iterate(self) -> AcadosOcpIterate:
        """
        Create a default initial iterate for the OCP.
        """
        self.make_consistent()
        dims = self.dims

        if self.constraints.has_x0:
            x_traj = (self.solver_options.N_horizon+1) * [self.constraints.x0]
        else:
            x_traj = (self.solver_options.N_horizon+1) * [np.zeros(dims.nx)]
        u_traj = self.solver_options.N_horizon * [np.zeros(self.dims.nu)]
        z_traj = self.solver_options.N_horizon * [np.zeros(self.dims.nz)]
        sl_traj = [np.zeros(self.dims.ns_0)] + (self.solver_options.N_horizon-1) * [np.zeros(self.dims.ns)] + [np.zeros(self.dims.ns_e)]
        su_traj = [np.zeros(self.dims.ns_0)] + (self.solver_options.N_horizon-1) * [np.zeros(self.dims.ns)] + [np.zeros(self.dims.ns_e)]

        pi_traj = self.solver_options.N_horizon * [np.zeros(self.dims.nx)]

        ni_0 = dims.nbu + dims.nbx_0 + dims.nh_0 + dims.nphi_0 + dims.ng + dims.ns_0
        ni = dims.nbu + dims.nbx + dims.nh + dims.nphi + dims.ng + dims.ns
        ni_e = dims.nbx_e + dims.nh_e + dims.nphi_e + dims.ng_e + dims.ns_e
        lam_traj = [np.zeros(2*ni_0)] + (self.solver_options.N_horizon-1) * [np.zeros(2*ni)] + [np.zeros(2*ni_e)]

        iterate = AcadosOcpIterate(
            x_traj=x_traj,
            u_traj=u_traj,
            z_traj=z_traj,
            sl_traj=sl_traj,
            su_traj=su_traj,
            pi_traj=pi_traj,
            lam_traj=lam_traj,
        )
        return iterate
