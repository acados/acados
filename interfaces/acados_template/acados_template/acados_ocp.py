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

from typing import Optional, Union
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

from .utils import (get_acados_path, format_class_dict, make_object_json_dumpable, render_template,
                    get_shared_lib_ext, is_column, is_empty, casadi_length, check_if_square,
                    check_casadi_version)
from .penalty_utils import symmetric_huber_penalty, one_sided_huber_penalty

from .zoro_description import ZoroDescription, process_zoro_description
from .casadi_function_generation import (
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
        self.__problem_class = 'OCP'

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
            self.__parameter_values = parameter_values
        else:
            raise Exception('Invalid parameter_values value. ' +
                            f'Expected numpy array, got {type(parameter_values)}.')


    def make_consistent(self) -> None:
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

        # parameters
        if self.parameter_values.shape[0] != dims.np:
            raise Exception('inconsistent dimension np, regarding model.p and parameter_values.' + \
                f'\nGot np = {dims.np}, self.parameter_values.shape = {self.parameter_values.shape[0]}\n')

        ## cost
        # initial stage - if not set, copy fields from path constraints
        if cost.cost_type_0 is None:
            self.copy_path_cost_to_stage_0()

        if cost.cost_type_0 == 'LINEAR_LS':
            check_if_square(cost.W_0, 'W_0')
            ny_0 = cost.W_0.shape[0]
            if cost.Vx_0.shape[0] != ny_0 or cost.Vu_0.shape[0] != ny_0:
                raise Exception('inconsistent dimension ny_0, regarding W_0, Vx_0, Vu_0.' + \
                                f'\nGot W_0[{cost.W_0.shape}], Vx_0[{cost.Vx_0.shape}], Vu_0[{cost.Vu_0.shape}]\n')
            if dims.nz != 0 and cost.Vz_0.shape[0] != ny_0:
                raise Exception('inconsistent dimension ny_0, regarding W_0, Vx_0, Vu_0, Vz_0.' + \
                                f'\nGot W_0[{cost.W_0.shape}], Vx_0[{cost.Vx_0.shape}], Vu_0[{cost.Vu_0.shape}], Vz_0[{cost.Vz_0.shape}]\n')
            if cost.Vx_0.shape[1] != dims.nx and ny_0 != 0:
                raise Exception('inconsistent dimension: Vx_0 should have nx columns.')
            if cost.Vu_0.shape[1] != dims.nu and ny_0 != 0:
                raise Exception('inconsistent dimension: Vu_0 should have nu columns.')
            if cost.yref_0.shape[0] != ny_0:
                raise Exception('inconsistent dimension: regarding W_0, yref_0.' + \
                                f'\nGot W_0[{cost.W_0.shape}], yref_0[{cost.yref_0.shape}]\n')
            dims.ny_0 = ny_0

        elif cost.cost_type_0 == 'NONLINEAR_LS':
            ny_0 = cost.W_0.shape[0]
            check_if_square(cost.W_0, 'W_0')
            if (is_empty(model.cost_y_expr_0) and ny_0 != 0) or casadi_length(model.cost_y_expr_0) != ny_0 or cost.yref_0.shape[0] != ny_0:
                raise Exception('inconsistent dimension ny_0: regarding W_0, cost_y_expr.' +
                                f'\nGot W_0[{cost.W_0.shape}], yref_0[{cost.yref_0.shape}], ',
                                f'cost_y_expr_0 [{casadi_length(model.cost_y_expr_0)}]\n')
            dims.ny_0 = ny_0

        elif cost.cost_type_0 == 'CONVEX_OVER_NONLINEAR':
            if is_empty(model.cost_y_expr_0):
                raise Exception('cost_y_expr_0 and/or cost_y_expr not provided.')
            ny_0 = casadi_length(model.cost_y_expr_0)
            if is_empty(model.cost_r_in_psi_expr_0) or casadi_length(model.cost_r_in_psi_expr_0) != ny_0:
                raise Exception('inconsistent dimension ny_0: regarding cost_y_expr_0 and cost_r_in_psi_0.')
            if is_empty(model.cost_psi_expr_0) or casadi_length(model.cost_psi_expr_0) != 1:
                raise Exception('cost_psi_expr_0 not provided or not scalar-valued.')
            if cost.yref_0.shape[0] != ny_0:
                raise Exception('inconsistent dimension: regarding yref_0 and cost_y_expr_0, cost_r_in_psi_0.')
            dims.ny_0 = ny_0

            if not (opts.hessian_approx=='EXACT' and opts.exact_hess_cost==False) and opts.hessian_approx != 'GAUSS_NEWTON':
                raise Exception("\nWith CONVEX_OVER_NONLINEAR cost type, possible Hessian approximations are:\n"
                "GAUSS_NEWTON or EXACT with 'exact_hess_cost' == False.\n")

        elif cost.cost_type_0 == 'EXTERNAL':
            if opts.hessian_approx == 'GAUSS_NEWTON' and opts.ext_cost_num_hess == 0 and model.cost_expr_ext_cost_custom_hess_0 is None:
                print("\nWARNING: Gauss-Newton Hessian approximation with EXTERNAL cost type not possible!\n"
                "got cost_type_0: EXTERNAL, hessian_approx: 'GAUSS_NEWTON.'\n"
                "GAUSS_NEWTON hessian is not defined for EXTERNAL cost formulation.\n"
                "If you continue, acados will proceed computing the exact hessian for the cost term.\n"
                "Note: There is also the option to use the external cost module with a numerical hessian approximation (see `ext_cost_num_hess`).\n"
                "OR the option to provide a symbolic custom hessian approximation (see `cost_expr_ext_cost_custom_hess`).\n")

        # path
        if cost.cost_type == 'LINEAR_LS':
            ny = cost.W.shape[0]
            check_if_square(cost.W, 'W')
            if cost.Vx.shape[0] != ny or cost.Vu.shape[0] != ny:
                raise Exception('inconsistent dimension ny, regarding W, Vx, Vu.' + \
                                f'\nGot W[{cost.W.shape}], Vx[{cost.Vx.shape}], Vu[{cost.Vu.shape}]\n')
            if dims.nz != 0 and cost.Vz.shape[0] != ny:
                raise Exception('inconsistent dimension ny, regarding W, Vx, Vu, Vz.' + \
                                f'\nGot W[{cost.W.shape}], Vx[{cost.Vx.shape}], Vu[{cost.Vu.shape}], Vz[{cost.Vz.shape}]\n')
            if cost.Vx.shape[1] != dims.nx and ny != 0:
                raise Exception('inconsistent dimension: Vx should have nx columns.')
            if cost.Vu.shape[1] != dims.nu and ny != 0:
                raise Exception('inconsistent dimension: Vu should have nu columns.')
            if cost.yref.shape[0] != ny:
                raise Exception('inconsistent dimension: regarding W, yref.' + \
                                f'\nGot W[{cost.W.shape}], yref[{cost.yref.shape}]\n')
            dims.ny = ny

        elif cost.cost_type == 'NONLINEAR_LS':
            ny = cost.W.shape[0]
            check_if_square(cost.W, 'W')
            if (is_empty(model.cost_y_expr) and ny != 0) or casadi_length(model.cost_y_expr) != ny or cost.yref.shape[0] != ny:
                raise Exception('inconsistent dimension: regarding W, yref.' + \
                                f'\nGot W[{cost.W.shape}], yref[{cost.yref.shape}],',
                                f'cost_y_expr[{casadi_length(model.cost_y_expr)}]\n')
            dims.ny = ny

        elif cost.cost_type == 'CONVEX_OVER_NONLINEAR':
            if is_empty(model.cost_y_expr):
                raise Exception('cost_y_expr and/or cost_y_expr not provided.')
            ny = casadi_length(model.cost_y_expr)
            if is_empty(model.cost_r_in_psi_expr) or casadi_length(model.cost_r_in_psi_expr) != ny:
                raise Exception('inconsistent dimension ny: regarding cost_y_expr and cost_r_in_psi.')
            if is_empty(model.cost_psi_expr) or casadi_length(model.cost_psi_expr) != 1:
                raise Exception('cost_psi_expr not provided or not scalar-valued.')
            if cost.yref.shape[0] != ny:
                raise Exception('inconsistent dimension: regarding yref and cost_y_expr, cost_r_in_psi.')
            dims.ny = ny

            if not (opts.hessian_approx=='EXACT' and opts.exact_hess_cost==False) and opts.hessian_approx != 'GAUSS_NEWTON':
                raise Exception("\nWith CONVEX_OVER_NONLINEAR cost type, possible Hessian approximations are:\n"
                "GAUSS_NEWTON or EXACT with 'exact_hess_cost' == False.\n")

        elif cost.cost_type == 'EXTERNAL':
            if opts.hessian_approx == 'GAUSS_NEWTON' and opts.ext_cost_num_hess == 0 and model.cost_expr_ext_cost_custom_hess is None:
                print("\nWARNING: Gauss-Newton Hessian approximation with EXTERNAL cost type not possible!\n"
                "got cost_type: EXTERNAL, hessian_approx: 'GAUSS_NEWTON.'\n"
                "GAUSS_NEWTON hessian is only supported for cost_types [NON]LINEAR_LS.\n"
                "If you continue, acados will proceed computing the exact hessian for the cost term.\n"
                "Note: There is also the option to use the external cost module with a numerical hessian approximation (see `ext_cost_num_hess`).\n"
                "OR the option to provide a symbolic custom hessian approximation (see `cost_expr_ext_cost_custom_hess`).\n")

        # terminal
        if cost.cost_type_e == 'LINEAR_LS':
            ny_e = cost.W_e.shape[0]
            check_if_square(cost.W_e, 'W_e')
            if cost.Vx_e.shape[0] != ny_e:
                raise Exception('inconsistent dimension ny_e: regarding W_e, cost_y_expr_e.' + \
                    f'\nGot W_e[{cost.W_e.shape}], Vx_e[{cost.Vx_e.shape}]')
            if cost.Vx_e.shape[1] != dims.nx and ny_e != 0:
                raise Exception('inconsistent dimension: Vx_e should have nx columns.')
            if cost.yref_e.shape[0] != ny_e:
                raise Exception('inconsistent dimension: regarding W_e, yref_e.')
            dims.ny_e = ny_e

        elif cost.cost_type_e == 'NONLINEAR_LS':
            ny_e = cost.W_e.shape[0]
            check_if_square(cost.W_e, 'W_e')
            if (is_empty(model.cost_y_expr_e) and ny_e != 0) or casadi_length(model.cost_y_expr_e) != ny_e or cost.yref_e.shape[0] != ny_e:
                raise Exception('inconsistent dimension ny_e: regarding W_e, cost_y_expr.' +
                                f'\nGot W_e[{cost.W_e.shape}], yref_e[{cost.yref_e.shape}], ',
                                f'cost_y_expr_e [{casadi_length(model.cost_y_expr_e)}]\n')
            dims.ny_e = ny_e

        elif cost.cost_type_e == 'CONVEX_OVER_NONLINEAR':
            if is_empty(model.cost_y_expr_e):
                raise Exception('cost_y_expr_e not provided.')
            ny_e = casadi_length(model.cost_y_expr_e)
            if is_empty(model.cost_r_in_psi_expr_e) or casadi_length(model.cost_r_in_psi_expr_e) != ny_e:
                raise Exception('inconsistent dimension ny_e: regarding cost_y_expr_e and cost_r_in_psi_e.')
            if is_empty(model.cost_psi_expr_e) or casadi_length(model.cost_psi_expr_e) != 1:
                raise Exception('cost_psi_expr_e not provided or not scalar-valued.')
            if cost.yref_e.shape[0] != ny_e:
                raise Exception('inconsistent dimension: regarding yref_e and cost_y_expr_e, cost_r_in_psi_e.')
            dims.ny_e = ny_e

            if not (opts.hessian_approx=='EXACT' and opts.exact_hess_cost==False) and opts.hessian_approx != 'GAUSS_NEWTON':
                raise Exception("\nWith CONVEX_OVER_NONLINEAR cost type, possible Hessian approximations are:\n"
                "GAUSS_NEWTON or EXACT with 'exact_hess_cost' == False.\n")

        elif cost.cost_type_e == 'EXTERNAL':
            if opts.hessian_approx == 'GAUSS_NEWTON' and opts.ext_cost_num_hess == 0 and model.cost_expr_ext_cost_custom_hess_e is None:
                print("\nWARNING: Gauss-Newton Hessian approximation with EXTERNAL cost type not possible!\n"
                "got cost_type_e: EXTERNAL, hessian_approx: 'GAUSS_NEWTON.'\n"
                "GAUSS_NEWTON hessian is only supported for cost_types [NON]LINEAR_LS.\n"
                "If you continue, acados will proceed computing the exact hessian for the cost term.\n"
                "Note: There is also the option to use the external cost module with a numerical hessian approximation (see `ext_cost_num_hess`).\n"
                "OR the option to provide a symbolic custom hessian approximation (see `cost_expr_ext_cost_custom_hess`).\n")

        supports_cost_integration = lambda type : type in ['NONLINEAR_LS', 'CONVEX_OVER_NONLINEAR']
        if opts.cost_discretization == 'INTEGRATOR' and \
            any([not supports_cost_integration(cost) for cost in [cost.cost_type_0, cost.cost_type, cost.cost_type_e]]):
            raise Exception('cost_discretization == INTEGRATOR only works with cost in ["NONLINEAR_LS", "CONVEX_OVER_NONLINEAR"] costs.')

        ## constraints
        # initial
        this_shape = constraints.lbx_0.shape
        other_shape = constraints.ubx_0.shape
        if not this_shape == other_shape:
            raise Exception('lbx_0, ubx_0 have different shapes!')
        if not is_column(constraints.lbx_0):
            raise Exception('lbx_0, ubx_0 must be column vectors!')
        dims.nbx_0 = constraints.lbx_0.size

        if constraints.has_x0 and dims.nbx_0 != dims.nx:
            raise Exception(f"x0 should have shape nx = {dims.nx}.")

        if constraints.has_x0:
            # case: x0 was set: nbx0 are all equalities.
            dims.nbxe_0 = dims.nbx_0
        elif constraints.idxbxe_0 is not None:
            dims.nbxe_0 = constraints.idxbxe_0.shape[0]
            if any(constraints.idxbxe_0 > dims.nbx_0):
                raise Exception(f'idxbxe_0 = {constraints.idxbxe_0} contains value > nbx_0 = {dims.nbx_0}.')
        elif dims.nbxe_0 is None:
            # case: x0 and idxbxe_0 were not set -> dont assume nbx0 to be equality constraints.
            dims.nbxe_0 = 0

        if not is_empty(model.con_h_expr_0):
            nh_0 = casadi_length(model.con_h_expr_0)
        else:
            nh_0 = 0

        if constraints.uh_0.shape[0] != nh_0 or constraints.lh_0.shape[0] != nh_0:
            raise Exception('inconsistent dimension nh_0, regarding lh_0, uh_0, con_h_expr_0.')
        else:
            dims.nh_0 = nh_0

        if is_empty(model.con_phi_expr_0):
            dims.nphi_0 = 0
            dims.nr_0 = 0
        else:
            dims.nphi_0 = casadi_length(model.con_phi_expr_0)
            constraints.constr_type_0 = "BGP"
            if is_empty(model.con_r_expr_0):
                raise Exception('convex over nonlinear constraints: con_r_expr_0 but con_phi_expr_0 is nonempty')
            else:
                dims.nr_0 = casadi_length(model.con_r_expr_0)

        # path
        nbx = constraints.idxbx.shape[0]
        if constraints.ubx.shape[0] != nbx or constraints.lbx.shape[0] != nbx:
            raise Exception('inconsistent dimension nbx, regarding idxbx, ubx, lbx.')
        else:
            dims.nbx = nbx
        if any(constraints.idxbx > dims.nx):
            raise Exception(f'idxbx = {constraints.idxbx} contains value > nx = {dims.nx}.')

        nbu = constraints.idxbu.shape[0]
        if constraints.ubu.shape[0] != nbu or constraints.lbu.shape[0] != nbu:
            raise Exception('inconsistent dimension nbu, regarding idxbu, ubu, lbu.')
        else:
            dims.nbu = nbu
        if any(constraints.idxbu > dims.nu):
            raise Exception(f'idxbu = {constraints.idxbu} contains value > nu = {dims.nu}.')

        # lg <= C * x + D * u <= ug
        ng = constraints.lg.shape[0]
        if constraints.ug.shape[0] != ng or constraints.C.shape[0] != ng \
        or constraints.D.shape[0] != ng:
            raise Exception('inconsistent dimension ng, regarding lg, ug, C, D.')
        else:
            dims.ng = ng

        if ng > 0:
            if constraints.C.shape[1] != dims.nx:
                raise Exception(f'inconsistent dimension nx, regarding C, got C.shape[1] = {constraints.C.shape[1]}.')
            if constraints.D.shape[1] != dims.nu:
                raise Exception(f'inconsistent dimension nu, regarding D, got D.shape[1] = {constraints.D.shape[1]}.')

        if not is_empty(model.con_h_expr):
            nh = casadi_length(model.con_h_expr)
        else:
            nh = 0

        if constraints.uh.shape[0] != nh or constraints.lh.shape[0] != nh:
            raise Exception('inconsistent dimension nh, regarding lh, uh, con_h_expr.')
        else:
            dims.nh = nh

        if is_empty(model.con_phi_expr):
            dims.nphi = 0
            dims.nr = 0
        else:
            dims.nphi = casadi_length(model.con_phi_expr)
            constraints.constr_type = "BGP"
            if is_empty(model.con_r_expr):
                raise Exception('convex over nonlinear constraints: con_r_expr but con_phi_expr is nonempty')
            else:
                dims.nr = casadi_length(model.con_r_expr)


        # terminal
        nbx_e = constraints.idxbx_e.shape[0]
        if constraints.ubx_e.shape[0] != nbx_e or constraints.lbx_e.shape[0] != nbx_e:
            raise Exception('inconsistent dimension nbx_e, regarding idxbx_e, ubx_e, lbx_e.')
        else:
            dims.nbx_e = nbx_e
        if any(constraints.idxbx_e > dims.nx):
            raise Exception(f'idxbx_e = {constraints.idxbx_e} contains value > nx = {dims.nx}.')

        ng_e = constraints.lg_e.shape[0]
        if constraints.ug_e.shape[0] != ng_e or constraints.C_e.shape[0] != ng_e:
            raise Exception('inconsistent dimension ng_e, regarding_e lg_e, ug_e, C_e.')
        else:
            dims.ng_e = ng_e

        if not is_empty(model.con_h_expr_e):
            nh_e = casadi_length(model.con_h_expr_e)
        else:
            nh_e = 0

        if constraints.uh_e.shape[0] != nh_e or constraints.lh_e.shape[0] != nh_e:
            raise Exception('inconsistent dimension nh_e, regarding lh_e, uh_e, con_h_expr_e.')
        else:
            dims.nh_e = nh_e

        if is_empty(model.con_phi_expr_e):
            dims.nphi_e = 0
            dims.nr_e = 0
        else:
            dims.nphi_e = casadi_length(model.con_phi_expr_e)
            constraints.constr_type_e = "BGP"
            if is_empty(model.con_r_expr_e):
                raise Exception('convex over nonlinear constraints: con_r_expr_e but con_phi_expr_e is nonempty')
            else:
                dims.nr_e = casadi_length(model.con_r_expr_e)

        # Slack dimensions
        nsbx = constraints.idxsbx.shape[0]
        if nsbx > nbx:
            raise Exception(f'inconsistent dimension nsbx = {nsbx}. Is greater than nbx = {nbx}.')
        if any(constraints.idxsbx > nbx):
            raise Exception(f'idxsbx = {constraints.idxsbx} contains value > nbx = {nbx}.')
        if is_empty(constraints.lsbx):
            constraints.lsbx = np.zeros((nsbx,))
        elif constraints.lsbx.shape[0] != nsbx:
            raise Exception('inconsistent dimension nsbx, regarding idxsbx, lsbx.')
        if is_empty(constraints.usbx):
            constraints.usbx = np.zeros((nsbx,))
        elif constraints.usbx.shape[0] != nsbx:
            raise Exception('inconsistent dimension nsbx, regarding idxsbx, usbx.')
        dims.nsbx = nsbx

        nsbu = constraints.idxsbu.shape[0]
        if nsbu > nbu:
            raise Exception(f'inconsistent dimension nsbu = {nsbu}. Is greater than nbu = {nbu}.')
        if any(constraints.idxsbu > nbu):
            raise Exception(f'idxsbu = {constraints.idxsbu} contains value > nbu = {nbu}.')
        if is_empty(constraints.lsbu):
            constraints.lsbu = np.zeros((nsbu,))
        elif constraints.lsbu.shape[0] != nsbu:
            raise Exception('inconsistent dimension nsbu, regarding idxsbu, lsbu.')
        if is_empty(constraints.usbu):
            constraints.usbu = np.zeros((nsbu,))
        elif constraints.usbu.shape[0] != nsbu:
            raise Exception('inconsistent dimension nsbu, regarding idxsbu, usbu.')
        dims.nsbu = nsbu

        nsh = constraints.idxsh.shape[0]
        if nsh > nh:
            raise Exception(f'inconsistent dimension nsh = {nsh}. Is greater than nh = {nh}.')
        if any(constraints.idxsh > nh):
            raise Exception(f'idxsh = {constraints.idxsh} contains value > nh = {nh}.')
        if is_empty(constraints.lsh):
            constraints.lsh = np.zeros((nsh,))
        elif constraints.lsh.shape[0] != nsh:
            raise Exception('inconsistent dimension nsh, regarding idxsh, lsh.')
        if is_empty(constraints.ush):
            constraints.ush = np.zeros((nsh,))
        elif constraints.ush.shape[0] != nsh:
            raise Exception('inconsistent dimension nsh, regarding idxsh, ush.')
        dims.nsh = nsh

        nsphi = constraints.idxsphi.shape[0]
        if nsphi > dims.nphi:
            raise Exception(f'inconsistent dimension nsphi = {nsphi}. Is greater than nphi = {dims.nphi}.')
        if any(constraints.idxsphi > dims.nphi):
            raise Exception(f'idxsphi = {constraints.idxsphi} contains value > nphi = {dims.nphi}.')
        if is_empty(constraints.lsphi):
            constraints.lsphi = np.zeros((nsphi,))
        elif constraints.lsphi.shape[0] != nsphi:
            raise Exception('inconsistent dimension nsphi, regarding idxsphi, lsphi.')
        if is_empty(constraints.usphi):
            constraints.usphi = np.zeros((nsphi,))
        elif constraints.usphi.shape[0] != nsphi:
            raise Exception('inconsistent dimension nsphi, regarding idxsphi, usphi.')
        dims.nsphi = nsphi

        nsg = constraints.idxsg.shape[0]
        if nsg > ng:
            raise Exception(f'inconsistent dimension nsg = {nsg}. Is greater than ng = {ng}.')
        if any(constraints.idxsg > ng):
            raise Exception(f'idxsg = {constraints.idxsg} contains value > ng = {ng}.')
        if is_empty(constraints.lsg):
            constraints.lsg = np.zeros((nsg,))
        elif constraints.lsg.shape[0] != nsg:
            raise Exception('inconsistent dimension nsg, regarding idxsg, lsg.')
        if is_empty(constraints.usg):
            constraints.usg = np.zeros((nsg,))
        elif constraints.usg.shape[0] != nsg:
            raise Exception('inconsistent dimension nsg, regarding idxsg, usg.')
        dims.nsg = nsg

        ns = nsbx + nsbu + nsh + nsg + nsphi
        wrong_fields = []
        if cost.Zl.shape[0] != ns:
            wrong_fields += ["Zl"]
            dim = cost.Zl.shape[0]
        elif cost.Zu.shape[0] != ns:
            wrong_fields += ["Zu"]
            dim = cost.Zu.shape[0]
        elif cost.zl.shape[0] != ns:
            wrong_fields += ["zl"]
            dim = cost.zl.shape[0]
        elif cost.zu.shape[0] != ns:
            wrong_fields += ["zu"]
            dim = cost.zu.shape[0]

        if wrong_fields != []:
            raise Exception(f'Inconsistent size for fields {", ".join(wrong_fields)}, with dimension {dim}, \n\t'\
                + f'Detected ns = {ns} = nsbx + nsbu + nsg + nsh + nsphi.\n\t'\
                + f'With nsbx = {nsbx}, nsbu = {nsbu}, nsg = {nsg}, nsh = {nsh}, nsphi = {nsphi}.')
        dims.ns = ns

        # slack dimensions at initial node
        nsh_0 = constraints.idxsh_0.shape[0]
        if nsh_0 > nh_0:
            raise Exception(f'inconsistent dimension nsh_0 = {nsh_0}. Is greater than nh_0 = {nh_0}.')
        if any(constraints.idxsh_0 > nh_0):
            raise Exception(f'idxsh_0 = {constraints.idxsh_0} contains value > nh_0 = {nh_0}.')
        if is_empty(constraints.lsh_0):
            constraints.lsh_0 = np.zeros((nsh_0,))
        elif constraints.lsh_0.shape[0] != nsh_0:
            raise Exception('inconsistent dimension nsh_0, regarding idxsh_0, lsh_0.')
        if is_empty(constraints.ush_0):
            constraints.ush_0 = np.zeros((nsh_0,))
        elif constraints.ush_0.shape[0] != nsh_0:
            raise Exception('inconsistent dimension nsh_0, regarding idxsh_0, ush_0.')
        dims.nsh_0 = nsh_0

        nsphi_0 = constraints.idxsphi_0.shape[0]
        if nsphi_0 > dims.nphi_0:
            raise Exception(f'inconsistent dimension nsphi_0 = {nsphi_0}. Is greater than nphi_0 = {dims.nphi_0}.')
        if any(constraints.idxsphi_0 > dims.nphi_0):
            raise Exception(f'idxsphi_0 = {constraints.idxsphi_0} contains value > nphi_0 = {dims.nphi_0}.')
        if is_empty(constraints.lsphi_0):
            constraints.lsphi_0 = np.zeros((nsphi_0,))
        elif constraints.lsphi_0.shape[0] != nsphi_0:
            raise Exception('inconsistent dimension nsphi_0, regarding idxsphi_0, lsphi_0.')
        if is_empty(constraints.usphi_0):
            constraints.usphi_0 = np.zeros((nsphi_0,))
        elif constraints.usphi_0.shape[0] != nsphi_0:
            raise Exception('inconsistent dimension nsphi_0, regarding idxsphi_0, usphi_0.')
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
                print("Using entries [zl, zu, Zl, Zu] at intial node for slack penalties.\n")
            else:
                raise ValueError("Fields cost.[zl_0, zu_0, Zl_0, Zu_0] are not provided and cannot be inferred from other fields.\n")

        wrong_fields = []
        if cost.Zl_0.shape[0] != ns_0:
            wrong_fields += ["Zl_0"]
            dim = cost.Zl_0.shape[0]
        elif cost.Zu_0.shape[0] != ns_0:
            wrong_fields += ["Zu_0"]
            dim = cost.Zu_0.shape[0]
        elif cost.zl_0.shape[0] != ns_0:
            wrong_fields += ["zl_0"]
            dim = cost.zl_0.shape[0]
        elif cost.zu_0.shape[0] != ns_0:
            wrong_fields += ["zu_0"]
            dim = cost.zu_0.shape[0]

        if wrong_fields != []:
            raise Exception(f'Inconsistent size for fields {", ".join(wrong_fields)}, with dimension {dim}, \n\t'\
                + f'Detected ns_0 = {ns_0} = nsbu + nsg + nsh_0 + nsphi_0.\n\t'\
                + f'With nsbu = {nsbu}, nsg = {nsg}, nsh_0 = {nsh_0}, nsphi_0 = {nsphi_0}.')
        dims.ns_0 = ns_0

        # slacks at terminal node
        nsbx_e = constraints.idxsbx_e.shape[0]
        if nsbx_e > nbx_e:
            raise Exception(f'inconsistent dimension nsbx_e = {nsbx_e}. Is greater than nbx_e = {nbx_e}.')
        if any(constraints.idxsbx_e > nbx_e):
            raise Exception(f'idxsbx_e = {constraints.idxsbx_e} contains value > nbx_e = {nbx_e}.')
        if is_empty(constraints.lsbx_e):
            constraints.lsbx_e = np.zeros((nsbx_e,))
        elif constraints.lsbx_e.shape[0] != nsbx_e:
            raise Exception('inconsistent dimension nsbx_e, regarding idxsbx_e, lsbx_e.')
        if is_empty(constraints.usbx_e):
            constraints.usbx_e = np.zeros((nsbx_e,))
        elif constraints.usbx_e.shape[0] != nsbx_e:
            raise Exception('inconsistent dimension nsbx_e, regarding idxsbx_e, usbx_e.')
        dims.nsbx_e = nsbx_e

        nsh_e = constraints.idxsh_e.shape[0]
        if nsh_e > nh_e:
            raise Exception(f'inconsistent dimension nsh_e = {nsh_e}. Is greater than nh_e = {nh_e}.')
        if nsh_e > nh_e:
            raise Exception(f'inconsistent dimension nsh_e = {nsh_e}. Is greater than nh_e = {nh_e}.')
        if is_empty(constraints.lsh_e):
            constraints.lsh_e = np.zeros((nsh_e,))
        elif constraints.lsh_e.shape[0] != nsh_e:
            raise Exception('inconsistent dimension nsh_e, regarding idxsh_e, lsh_e.')
        if is_empty(constraints.ush_e):
            constraints.ush_e = np.zeros((nsh_e,))
        elif constraints.ush_e.shape[0] != nsh_e:
            raise Exception('inconsistent dimension nsh_e, regarding idxsh_e, ush_e.')
        dims.nsh_e = nsh_e

        nsphi_e = constraints.idxsphi_e.shape[0]
        if nsphi_e > dims.nphi_e:
            raise Exception(f'inconsistent dimension nsphi_e = {nsphi_e}. Is greater than nphi_e = {dims.nphi_e}.')
        if nsphi_e > dims.nphi_e:
            raise Exception(f'inconsistent dimension nsphi_e = {nsphi_e}. Is greater than nphi_e = {dims.nphi_e}.')
        if is_empty(constraints.lsphi_e):
            constraints.lsphi_e = np.zeros((nsphi_e,))
        elif constraints.lsphi_e.shape[0] != nsphi_e:
            raise Exception('inconsistent dimension nsphi_e, regarding idxsphi_e, lsphi_e.')
        if is_empty(constraints.usphi_e):
            constraints.usphi_e = np.zeros((nsphi_e,))
        elif constraints.usphi_e.shape[0] != nsphi_e:
            raise Exception('inconsistent dimension nsphi_e, regarding idxsphi_e, usphi_e.')
        dims.nsphi_e = nsphi_e

        nsg_e = constraints.idxsg_e.shape[0]
        if nsg_e > ng_e:
            raise Exception(f'inconsistent dimension nsg_e = {nsg_e}. Is greater than ng_e = {ng_e}.')
        if nsg_e > ng_e:
            raise Exception(f'inconsistent dimension nsg_e = {nsg_e}. Is greater than ng_e = {ng_e}.')
        if is_empty(constraints.lsg_e):
            constraints.lsg_e = np.zeros((nsg_e,))
        elif constraints.lsg_e.shape[0] != nsg_e:
            raise Exception('inconsistent dimension nsg_e, regarding idxsg_e, lsg_e.')
        if is_empty(constraints.usg_e):
            constraints.usg_e = np.zeros((nsg_e,))
        elif constraints.usg_e.shape[0] != nsg_e:
            raise Exception('inconsistent dimension nsg_e, regarding idxsg_e, usg_e.')
        dims.nsg_e = nsg_e

        # terminal
        ns_e = nsbx_e + nsh_e + nsg_e + nsphi_e
        wrong_field = ""
        if cost.Zl_e.shape[0] != ns_e:
            wrong_field = "Zl_e"
            dim = cost.Zl_e.shape[0]
        elif cost.Zu_e.shape[0] != ns_e:
            wrong_field = "Zu_e"
            dim = cost.Zu_e.shape[0]
        elif cost.zl_e.shape[0] != ns_e:
            wrong_field = "zl_e"
            dim = cost.zl_e.shape[0]
        elif cost.zu_e.shape[0] != ns_e:
            wrong_field = "zu_e"
            dim = cost.zu_e.shape[0]

        if wrong_field != "":
            raise Exception(f'Inconsistent size for field {wrong_field}, with dimension {dim}, \n\t'\
                + f'Detected ns_e = {ns_e} = nsbx_e + nsg_e + nsh_e + nsphi_e.\n\t'\
                + f'With nsbx_e = {nsbx_e}, nsg_e = {nsg_e}, nsh_e = {nsh_e}, nsphi_e = {nsphi_e}.')

        dims.ns_e = ns_e

        # discretization
        if not isinstance(opts.tf, (float, int)):
            raise Exception(f'Time horizon tf should be float provided, got tf = {opts.tf}.')

        if is_empty(opts.time_steps) and is_empty(opts.shooting_nodes):
            # uniform discretization
            opts.time_steps = opts.tf / dims.N * np.ones((dims.N,))
            opts.shooting_nodes = np.concatenate((np.array([0.]), np.cumsum(opts.time_steps)))

        elif not is_empty(opts.shooting_nodes):
            if np.shape(opts.shooting_nodes)[0] != dims.N+1:
                raise Exception('inconsistent dimension N, regarding shooting_nodes.')

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
            Exception('Please provide either time_steps or shooting_nodes for nonuniform discretization')

        tf = np.sum(opts.time_steps)
        if (tf - opts.tf) / tf > 1e-13:
            raise Exception(f'Inconsistent discretization: {opts.tf}'\
                f' = tf != sum(opts.time_steps) = {tf}.')

        # num_steps
        if isinstance(opts.sim_method_num_steps, np.ndarray) and opts.sim_method_num_steps.size == 1:
            opts.sim_method_num_steps = opts.sim_method_num_steps.item()

        if isinstance(opts.sim_method_num_steps, (int, float)) and opts.sim_method_num_steps % 1 == 0:
            opts.sim_method_num_steps = opts.sim_method_num_steps * np.ones((dims.N,), dtype=np.int64)
        elif isinstance(opts.sim_method_num_steps, np.ndarray) and opts.sim_method_num_steps.size == dims.N \
            and np.all(np.equal(np.mod(opts.sim_method_num_steps, 1), 0)):
            opts.sim_method_num_steps = np.reshape(opts.sim_method_num_steps, (dims.N,)).astype(np.int64)
        else:
            raise Exception("Wrong value for sim_method_num_steps. Should be either int or array of ints of shape (N,).")

        # num_stages
        if isinstance(opts.sim_method_num_stages, np.ndarray) and opts.sim_method_num_stages.size == 1:
            opts.sim_method_num_stages = opts.sim_method_num_stages.item()

        if isinstance(opts.sim_method_num_stages, (int, float)) and opts.sim_method_num_stages % 1 == 0:
            opts.sim_method_num_stages = opts.sim_method_num_stages * np.ones((dims.N,), dtype=np.int64)
        elif isinstance(opts.sim_method_num_stages, np.ndarray) and opts.sim_method_num_stages.size == dims.N \
            and np.all(np.equal(np.mod(opts.sim_method_num_stages, 1), 0)):
            opts.sim_method_num_stages = np.reshape(opts.sim_method_num_stages, (dims.N,)).astype(np.int64)
        else:
            raise Exception("Wrong value for sim_method_num_stages. Should be either int or array of ints of shape (N,).")

        # jac_reuse
        if isinstance(opts.sim_method_jac_reuse, np.ndarray) and opts.sim_method_jac_reuse.size == 1:
            opts.sim_method_jac_reuse = opts.sim_method_jac_reuse.item()

        if isinstance(opts.sim_method_jac_reuse, (int, float)) and opts.sim_method_jac_reuse % 1 == 0:
            opts.sim_method_jac_reuse = opts.sim_method_jac_reuse * np.ones((dims.N,), dtype=np.int64)
        elif isinstance(opts.sim_method_jac_reuse, np.ndarray) and opts.sim_method_jac_reuse.size == dims.N \
            and np.all(np.equal(np.mod(opts.sim_method_jac_reuse, 1), 0)):
            opts.sim_method_jac_reuse = np.reshape(opts.sim_method_jac_reuse, (dims.N,)).astype(np.int64)
        else:
            raise Exception("Wrong value for sim_method_jac_reuse. Should be either int or array of ints of shape (N,).")

        # fixed hessian
        if opts.fixed_hess:
            if opts.hessian_approx == 'EXACT':
                raise Exception('fixed_hess is not compatible with hessian_approx == EXACT.')
            if cost.cost_type != "LINEAR_LS":
                raise Exception('fixed_hess is only compatible LINEAR_LS cost_type.')
            if cost.cost_type_0 != "LINEAR_LS":
                raise Exception('fixed_hess is only compatible LINEAR_LS cost_type_0.')
            if cost.cost_type_e != "LINEAR_LS":
                raise Exception('fixed_hess is only compatible LINEAR_LS cost_type_e.')

        # solution sensitivities
        type_constraint_pairs = [("path", model.con_h_expr), ("initial", model.con_h_expr_0),
                                 ("terminal", model.con_h_expr_e),
                                 ("path", model.con_phi_expr), ("initial", model.con_phi_expr_0), ("terminal", model.con_phi_expr_e),
                                 ("path", model.con_r_expr), ("initial", model.con_r_expr_0), ("terminal", model.con_r_expr_e)]

        if opts.with_solution_sens_wrt_params:
            if cost.cost_type != "EXTERNAL" or cost.cost_type_0 != "EXTERNAL" or cost.cost_type_e != "EXTERNAL":
                raise Exception('with_solution_sens_wrt_params is only compatible with EXTERNAL cost_type.')
            if opts.integrator_type != "DISCRETE":
                raise Exception('with_solution_sens_wrt_params is only compatible with DISCRETE dynamics.')
            for horizon_type, constraint in type_constraint_pairs:
                if constraint is not None and any(ca.which_depends(constraint, model.p)):
                    raise Exception(f'with_solution_sens_wrt_params is only implemented if constraints depend not on parameters. Got parameter dependency for {horizon_type} constraint.')

        if opts.with_value_sens_wrt_params:
            if cost.cost_type != "EXTERNAL" or cost.cost_type_0 != "EXTERNAL" or cost.cost_type_e != "EXTERNAL":
                raise Exception('with_value_sens_wrt_params is only compatible with EXTERNAL cost_type.')
            if opts.integrator_type != "DISCRETE":
                raise Exception('with_value_sens_wrt_params is only compatible with DISCRETE dynamics.')
            for horizon_type, constraint in type_constraint_pairs:
                if constraint is not None and any(ca.which_depends(constraint, model.p)):
                    raise Exception(f'with_value_sens_wrt_params is only implemented if constraints depend not on parameters. Got parameter dependency for {horizon_type} constraint.')

        if opts.qp_solver_cond_N is None:
            opts.qp_solver_cond_N = dims.N

        if opts.qp_solver_cond_block_size is not None:
            if sum(opts.qp_solver_cond_block_size) != dims.N:
                raise Exception(f'sum(qp_solver_cond_block_size) = {sum(opts.qp_solver_cond_block_size)} != N = {dims.N}.')
            if len(opts.qp_solver_cond_block_size) != opts.qp_solver_cond_N+1:
                raise Exception(f'qp_solver_cond_block_size = {opts.qp_solver_cond_block_size} should have length qp_solver_cond_N+1 = {opts.qp_solver_cond_N+1}.')

        if opts.nlp_solver_type == "DDP":
            if opts.qp_solver != "PARTIAL_CONDENSING_HPIPM" or opts.qp_solver_cond_N != dims.N:
                raise Exception(f'DDP solver only supported for PARTIAL_CONDENSING_HPIPM with qp_solver_cond_N == N, got qp solver {opts.qp_solver} and qp_solver_cond_N {opts.qp_solver_cond_N}, N {dims.N}.')
            if any([dims.nbu, dims.nbx, dims.ng, dims.nh, dims.nphi]):
                raise Exception('DDP only supports initial state constraints, got path constraints.')
            if  any([dims.ng_e, dims.nphi_e, dims.nh_e]):
                raise Exception('DDP only supports initial state constraints, got terminal constraints.')

        # Set default parameters for globalization
        if opts.alpha_min == None:
            if opts.globalization == 'FUNNEL_L1PEN_LINESEARCH':
                opts.alpha_min = 1e-17
            else:
                opts.alpha_min = 0.05

        if opts.alpha_reduction == None:
            if opts.globalization == 'FUNNEL_L1PEN_LINESEARCH':
                opts.alpha_reduction = 0.5
            else:
                opts.alpha_reduction = 0.7

        if opts.eps_sufficient_descent == None:
            if opts.globalization == 'FUNNEL_L1PEN_LINESEARCH':
                opts.eps_sufficient_descent = 1e-6
            else:
                opts.eps_sufficient_descent = 1e-4

        if opts.eval_residual_at_max_iter == None:
            if opts.globalization == 'FUNNEL_L1PEN_LINESEARCH':
                opts.eval_residual_at_max_iter = True
            else:
                opts.eval_residual_at_max_iter = False

        if opts.full_step_dual == None:
            if opts.globalization == 'FUNNEL_L1PEN_LINESEARCH':
                opts.full_step_dual = 1
            else:
                opts.full_step_dual = 0

        # sanity check for Funnel globalization and SQP
        if opts.globalization == 'FUNNEL_L1PEN_LINESEARCH' and opts.nlp_solver_type != 'SQP':
            raise Exception('FUNNEL_L1PEN_LINESEARCH only supports SQP.')

        # termination
        if opts.nlp_solver_tol_min_step_norm == None:
            if opts.globalization == 'FUNNEL_L1PEN_LINESEARCH':
                opts.nlp_solver_tol_min_step_norm = 1e-12
            else:
                opts.nlp_solver_tol_min_step_norm = 0.0

        # zoRO
        if self.zoro_description is not None:
            if not isinstance(self.zoro_description, ZoroDescription):
                raise Exception('zoro_description should be of type ZoroDescription or None')
            else:
                self.zoro_description = process_zoro_description(self.zoro_description)

        return


    def _get_external_function_header_templates(self, ) -> list:
        dims = self.dims
        name = self.model.name
        template_list = []

        # dynamics
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
        template_list.append(('acados_sim_solver.in.c', f'acados_sim_solver_{name}.c'))
        template_list.append(('acados_sim_solver.in.h', f'acados_sim_solver_{name}.h'))
        template_list.append(('main_sim.in.c', f'main_sim_{name}.c'))

        # model
        template_list += self._get_external_function_header_templates()

        # Simulink
        if self.simulink_opts is not None:
            template_file = os.path.join('matlab_templates', 'acados_solver_sfun.in.c')
            template_list.append((template_file, f'acados_solver_sfunction_{name}.c'))
            template_file = os.path.join('matlab_templates', 'make_sfun.in.m')
            template_list.append((template_file, f'make_sfun_{name}.m'))
            template_file = os.path.join('matlab_templates', 'acados_sim_solver_sfun.in.c')
            template_list.append((template_file, f'acados_sim_solver_sfunction_{name}.c'))
            template_file = os.path.join('matlab_templates', 'make_sfun_sim.in.m')
            template_list.append((template_file, f'make_sfun_sim_{name}.m'))

        return template_list


    def render_templates(self, json_file: str, cmake_builder=None):

        # check json file
        json_path = os.path.abspath(json_file)
        if not os.path.exists(json_path):
            raise Exception(f'Path "{json_path}" not found!')

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


    def dump_to_json(self, json_file: str) -> None:
        with open(json_file, 'w') as f:
            json.dump(self.to_dict(), f, default=make_object_json_dumpable, indent=4, sort_keys=True)
        return


    def generate_external_functions(self):
        model = self.model
        constraints = self.constraints

        # options for code generation
        code_gen_opts = dict()
        code_gen_opts['generate_hess'] = self.solver_options.hessian_approx == 'EXACT'
        code_gen_opts['with_solution_sens_wrt_params'] = self.solver_options.with_solution_sens_wrt_params
        code_gen_opts['with_value_sens_wrt_params'] = self.solver_options.with_value_sens_wrt_params
        code_gen_opts['code_export_directory'] = self.code_export_directory

        # create code_export_dir, model_dir
        model_dir = os.path.join(code_gen_opts['code_export_directory'], model.name + '_model')
        if not os.path.exists(model_dir):
            os.makedirs(model_dir)

        check_casadi_version()
        if self.model.dyn_ext_fun_type == 'casadi':
            if self.solver_options.integrator_type == 'ERK':
                generate_c_code_explicit_ode(model, code_gen_opts)
            elif self.solver_options.integrator_type == 'IRK':
                generate_c_code_implicit_ode(model, code_gen_opts)
            elif self.solver_options.integrator_type == 'LIFTED_IRK':
                if model.t != []:
                    raise NotImplementedError("LIFTED_IRK with time-varying dynamics not implemented yet.")
                generate_c_code_implicit_ode(model, code_gen_opts)
            elif self.solver_options.integrator_type == 'GNSF':
                generate_c_code_gnsf(model, code_gen_opts)
            elif self.solver_options.integrator_type == 'DISCRETE':
                generate_c_code_discrete_dynamics(model, code_gen_opts)
            else:
                raise Exception("ocp_generate_external_functions: unknown integrator type.")
        else:
            target_location = os.path.join(code_gen_opts['code_export_directory'], model_dir, model.dyn_generic_source)
            shutil.copyfile(model.dyn_generic_source, target_location)

        stage_types = ['initial', 'path', 'terminal']

        for attr_nh, attr_nphi, stage_type in zip(['nh_0', 'nh', 'nh_e'], ['nphi_0', 'nphi', 'nphi_e'], stage_types):
            if getattr(self.dims, attr_nh) > 0 or getattr(self.dims, attr_nphi) > 0:
                generate_c_code_constraint(model, constraints, stage_type, code_gen_opts)

        for attr, stage_type in zip(['cost_type_0', 'cost_type', 'cost_type_e'], stage_types):
            if getattr(self.cost, attr) == 'NONLINEAR_LS':
                generate_c_code_nls_cost(model, stage_type, code_gen_opts)
            elif getattr(self.cost, attr) == 'CONVEX_OVER_NONLINEAR':
                generate_c_code_conl_cost(model, stage_type, code_gen_opts)
            elif getattr(self.cost, attr) == 'EXTERNAL':
                generate_c_code_external_cost(model, stage_type, code_gen_opts)


    def remove_x0_elimination(self) -> None:
        self.constraints.idxbxe_0 = np.zeros((0,))
        self.dims.nbxe_0 = 0
        self.constraints.__has_x0 = False
        return


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
            raise Exception(f"Terminal cost type must be NONLINEAR_LS, got cost_type_0 {self.cost.cost_type_0}.")

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
            raise Exception(f"Path cost type must be NONLINEAR_LS, got cost_type {self.cost.cost_type}.")

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
            raise Exception(f"Initial cost type must be NONLINEAR_LS, got cost_type_e {self.cost.cost_type_e}.")
        return


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
            raise Exception("Huber penalty is only supported for CONVEX_OVER_NONLINEAR cost type.")

        if use_xgn and self.model.cost_conl_custom_outer_hess is None:
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
        elif self.model.cost_conl_custom_outer_hess is not None:
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
            raise Exception("Parameter t0 is already present in the model.")
        self.model.t0 = ca.SX.sym("t0")
        self.model.p = ca.vertcat(self.model.p, self.model.t0)
        self.parameter_values = np.append(self.parameter_values, [0.0])
        return
