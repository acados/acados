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

from typing import Union

import casadi as ca
import numpy as np

from casadi import MX, SX

from .utils import is_empty, casadi_length, is_column
from .acados_dims import AcadosOcpDims, AcadosSimDims


class AcadosModel():
    """
    Class containing all the information to code generate the external CasADi functions
    that are needed when creating an acados ocp solver or acados integrator.
    Thus, this class contains:

    a) the :py:attr:`name` of the model,
    b) all CasADi variables/expressions needed in the CasADi function generation process.
    """
    def __init__(self):
        ## common for OCP and Integrator
        self.name = None
        """The model name is used for code generation. Type: string. Default: :code:`None`"""
        self.x = None
        """CasADi variable describing the state of the system; Default: :code:`None`"""
        self.xdot = None
        """CasADi variable describing the derivative of the state wrt time; Default: :code:`None`"""
        self.u = None
        """CasADi variable describing the input of the system; Default: :code:`None`"""
        self.z = []
        """CasADi variable describing the algebraic variables of the DAE; Default: :code:`[]`"""
        self.p = []
        """CasADi variable describing parameters of the DAE; Default: :code:`[]`"""
        self.t = []
        """
        CasADi variable representing time t in functions; Default: :code:`[]`
        NOTE:
        - For integrators, the start time has to be explicitly set via :py:attr:`acados_template.AcadosSimSolver.set`('t0').
        - For OCPs, the start time is set to 0. on each stage.
        The time dependency can be used within cost formulations and is relevant when cost integration is used.
        Start times of shooting intervals can be added using parameters.
        """

        ## dynamics
        self.f_impl_expr = None
        """
        CasADi expression for the implicit dynamics :math:`f_\\text{impl}(\dot{x}, x, u, z, p) = 0`.
        Used if :py:attr:`acados_template.acados_ocp_options.AcadosOcpOptions.integrator_type` == 'IRK'.
        Default: :code:`None`
        """
        self.f_expl_expr = None
        """
        CasADi expression for the explicit dynamics :math:`\dot{x} = f_\\text{expl}(x, u, p)`.
        Used if :py:attr:`acados_template.acados_ocp_options.AcadosOcpOptions.integrator_type` == 'ERK'.
        Default: :code:`None`
        """
        self.disc_dyn_expr = None
        """
        CasADi expression for the discrete dynamics :math:`x_{+} = f_\\text{disc}(x, u, p)`.
        Used if :py:attr:`acados_template.acados_ocp_options.AcadosOcpOptions.integrator_type` == 'DISCRETE'.
        Default: :code:`None`
        """

        self.dyn_ext_fun_type = 'casadi'
        """type of external functions for dynamics module; 'casadi' or 'generic'; Default: 'casadi'"""
        self.dyn_generic_source = None
        """name of source file for discrete dynamics, only relevant if :code:`dyn_ext_fun_type` is :code:`'generic'`; Default: :code:`None`"""
        self.dyn_disc_fun_jac_hess = None
        """name of function discrete dynamics + jacobian and hessian, only relevant if :code:`dyn_ext_fun_type` is :code:`'generic'`; Default: :code:`None`"""
        self.dyn_disc_fun_jac = None
        """name of function discrete dynamics + jacobian, only relevant if :code:`dyn_ext_fun_type` is :code:`'generic'`; Default: :code:`None`"""
        self.dyn_disc_fun = None
        """name of function discrete dynamics, only relevant if :code:`dyn_ext_fun_type` is :code:`'generic'`; Default: :code:`None`"""

        self.dyn_impl_dae_fun_jac = None
        """name of source files for implicit DAE function value and jacobian, only relevant if :code:`dyn_ext_fun_type` is :code:`'generic'`; Default: :code:`None`"""
        self.dyn_impl_dae_jac = None
        """name of source files for implicit DAE jacobian, only relevant if :code:`dyn_ext_fun_type` is :code:`'generic'`; Default: :code:`None`"""
        self.dyn_impl_dae_fun = None
        """name of source files for implicit DAE function value, only relevant if :code:`dyn_ext_fun_type` is :code:`'generic'`; Default: :code:`None`"""

        # for GNSF models
        self.gnsf = {'nontrivial_f_LO': 1, 'purely_linear': 0}
        """
        dictionary containing information on GNSF structure needed when rendering templates.
        Contains integers `nontrivial_f_LO`, `purely_linear`.
        """

        ### for OCP only.
        # NOTE: These could be moved to cost / constraints

        # constraints at initial stage
        self.con_h_expr_0 = None
        """CasADi expression for the initial constraint :math:`h^0`; Default: :code:`None`"""
        self.con_phi_expr_0 = None
        """CasADi expression for the terminal constraint :math:`\phi_0`; Default: :code:`None`"""
        self.con_r_expr_0 = None
        """CasADi expression for the terminal constraint :math:`\phi_0(r)`,
        dummy input for outer function; Default: :code:`None`"""
        self.con_r_in_phi_0 = None
        """CasADi expression for the terminal constraint :math:`\phi_0(r)`, input for outer function; Default: :code:`None`"""


        # path constraints
        # BGH(default): lh <= h(x, u) <= uh
        self.con_h_expr = None
        """CasADi expression for the constraint :math:`h`; Default: :code:`None`"""
        # BGP(convex over nonlinear): lphi <= phi(r(x, u)) <= uphi
        self.con_phi_expr = None
        """CasADi expression for the constraint phi; Default: :code:`None`"""
        self.con_r_expr = None
        """CasADi expression for the constraint phi(r),
        dummy input for outer function; Default: :code:`None`"""
        self.con_r_in_phi = None
        """CasADi expression for the terminal constraint :math:`\phi(r)`,
        input for outer function; Default: :code:`None`"""

        # terminal
        self.con_h_expr_e = None
        """CasADi expression for the terminal constraint :math:`h^e`; Default: :code:`None`"""
        self.con_phi_expr_e = None
        """CasADi expression for the terminal constraint :math:`\phi_e`; Default: :code:`None`"""
        self.con_r_expr_e = None
        """CasADi expression for the terminal constraint :math:`\phi_e(r)`,
        dummy input for outer function; Default: :code:`None`"""
        self.con_r_in_phi_e = None
        """CasADi expression for the terminal constraint :math:`\phi_e(r)`, input for outer function; Default: :code:`None`"""

        # cost
        self.cost_y_expr = None
        """CasADi expression for nonlinear least squares; Default: :code:`None`"""
        self.cost_y_expr_e = None
        """CasADi expression for nonlinear least squares, terminal; Default: :code:`None`"""
        self.cost_y_expr_0 = None
        """CasADi expression for nonlinear least squares, initial; Default: :code:`None`"""
        self.cost_expr_ext_cost = None
        """CasADi expression for external cost; Default: :code:`None`"""
        self.cost_expr_ext_cost_e = None
        """CasADi expression for external cost, terminal; Default: :code:`None`"""
        self.cost_expr_ext_cost_0 = None
        """CasADi expression for external cost, initial; Default: :code:`None`"""
        self.cost_expr_ext_cost_custom_hess = None
        """CasADi expression for custom hessian (only for external cost); Default: :code:`None`"""
        self.cost_expr_ext_cost_custom_hess_e = None
        """CasADi expression for custom hessian (only for external cost), terminal; Default: :code:`None`"""
        self.cost_expr_ext_cost_custom_hess_0 = None
        """CasADi expression for custom hessian (only for external cost), initial; Default: :code:`None`"""

        ## CONVEX_OVER_NONLINEAR convex-over-nonlinear cost: psi(y(x, u, p) - y_ref; p)
        self.cost_psi_expr_0 = None
        """
        CasADi expression for the outer loss function :math:`\psi(r, p)`, initial; Default: :code:`None`
        Used if :py:attr:`acados_template.acados_ocp_options.AcadosOcpOptions.cost_type_0` is 'CONVEX_OVER_NONLINEAR'.
        """
        self.cost_psi_expr = None
        """
        CasADi expression for the outer loss function :math:`\psi(r, p)`; Default: :code:`None`
        Used if :py:attr:`acados_template.acados_ocp_options.AcadosOcpOptions.cost_type` is 'CONVEX_OVER_NONLINEAR'.
        """
        self.cost_psi_expr_e = None
        """
        CasADi expression for the outer loss function :math:`\psi(r, p)`, terminal; Default: :code:`None`
        Used if :py:attr:`acados_template.acados_ocp_options.AcadosOcpOptions.cost_type_e` is 'CONVEX_OVER_NONLINEAR'.
        """
        self.cost_r_in_psi_expr_0 = None
        """
        CasADi symbolic input variable for the argument :math:`r` to the outer loss function :math:`\psi(r, p)`, initial; Default: :code:`None`
        Used if :py:attr:`acados_template.acados_ocp_options.AcadosOcpOptions.cost_type_0` is 'CONVEX_OVER_NONLINEAR'.
        """
        self.cost_r_in_psi_expr = None
        """
        CasADi symbolic input variable for the argument :math:`r` to the outer loss function :math:`\psi(r, p)`; Default: :code:`None`
        Used if :py:attr:`acados_template.acados_ocp_options.AcadosOcpOptions.cost_type` is 'CONVEX_OVER_NONLINEAR'.
        """
        self.cost_r_in_psi_expr_e = None
        """
        CasADi symbolic input variable for the argument :math:`r` to the outer loss function :math:`\psi(r, p)`, terminal; Default: :code:`None`
        Used if :py:attr:`acados_template.acados_ocp_options.AcadosOcpOptions.cost_type_e` is 'CONVEX_OVER_NONLINEAR'.
        """
        self.cost_conl_custom_outer_hess_0 = None
        """
        CasADi expression for the custom hessian of the outer loss function (only for convex-over-nonlinear cost), initial; Default: :code:`None`
        Used if :py:attr:`acados_template.acados_ocp_options.AcadosOcpOptions.cost_type_0` is 'CONVEX_OVER_NONLINEAR'.
        """
        self.cost_conl_custom_outer_hess = None
        """
        CasADi expression for the custom hessian of the outer loss function (only for convex-over-nonlinear cost); Default: :code:`None`
        Used if :py:attr:`acados_template.acados_ocp_options.AcadosOcpOptions.cost_type` is 'CONVEX_OVER_NONLINEAR'.
        """
        self.cost_conl_custom_outer_hess_e = None
        """
        CasADi expression for the custom hessian of the outer loss function (only for convex-over-nonlinear cost), terminal; Default: :code:`None`
        Used if :py:attr:`acados_template.acados_ocp_options.AcadosOcpOptions.cost_type_e` is 'CONVEX_OVER_NONLINEAR'.
        """
        self.nu_original = None # TODO: remove? only used by benchmark
        """
        Number of original control inputs (before polynomial control augmentation); Default: :code:`None`
        """
        self.t0 = None # TODO: remove? only used by benchmark
        """CasADi variable representing the start time of an interval; Default: :code:`None`"""
        self.__x_labels = None
        self.__u_labels = None
        self.__t_label = "t"

    @property
    def x_labels(self):
        """Contains list of labels for the states. Default: :code:`None`"""
        if self.__x_labels is None:
            return [f"x{i}" for i in range(self.x.size()[0])]
        else:
            return self.__x_labels

    @x_labels.setter
    def x_labels(self, x_labels):
        self.__x_labels = x_labels


    @property
    def u_labels(self):
        """Contains list of labels for the controls. Default: :code:`None`"""
        if self.__u_labels is None:
            return [f"x{i}" for i in range(self.x.size()[0])]
        else:
            return self.__u_labels

    @u_labels.setter
    def u_labels(self, u_labels):
        self.__u_labels = u_labels

    @property
    def t_label(self):
        """Label for the time variable. Default: :code:'t'"""
        return self.__t_label

    @t_label.setter
    def t_label(self, t_label):
        self.__t_label = t_label


    def get_casadi_symbol(self):
        if isinstance(self.x, MX):
            return MX.sym
        elif isinstance(self.x, SX):
            return SX.sym
        else:
            raise Exception(f"model.x must be casadi.SX or casadi.MX, got {type(self.x)}")

    def get_casadi_zeros(self):
        if isinstance(self.x, MX):
            return MX.zeros
        elif isinstance(self.x, SX):
            return SX.zeros
        else:
            raise Exception(f"model.x must be casadi.SX or casadi.MX, got {type(self.x)}")

    def make_consistent(self, dims: Union[AcadosOcpDims, AcadosSimDims]) -> None:

        casadi_symbol = self.get_casadi_symbol()
        if is_empty(self.p):
            self.p = casadi_symbol('p', 0, 0)

        if is_empty(self.z):
            self.z = casadi_symbol('z', 0, 0)

        # nx
        if is_column(self.x):
            dims.nx = casadi_length(self.x)
        else:
            raise Exception('model.x should be column vector!')

        # nu
        if is_empty(self.u):
            dims.nu = 0
            self.u = casadi_symbol('u', 0, 1)
        else:
            dims.nu = casadi_length(self.u)

        # nz
        if is_empty(self.z):
            dims.nz = 0
        else:
            dims.nz = casadi_length(self.z)

        # np
        if is_empty(self.p):
            dims.np = 0
        else:
            dims.np = casadi_length(self.p)

        return


    def substitute(self, var: Union[ca.SX, ca.MX], expr_new: Union[ca.SX, ca.MX]) -> None:
        """
        Substitutes the variables var with expr_new in all symbolic CasADi expressions within AcadosModel
        """
        for attr, value in self.__dict__.items():
            if isinstance(value, (ca.SX, ca.MX)):
                new = ca.substitute(value, var, expr_new)
                setattr(self, attr, new)
        return


    def augment_model_with_polynomial_control(self, degree: int) -> None:
        print("Deprecation warning: augment_model_with_polynomial_control() is deprecated and has been renamed to reformulate_with_polynomial_control().")
        self.reformulate_with_polynomial_control(degree=degree)

    def reformulate_with_polynomial_control(self, degree: int) -> None:
        """
        Augment the model with polynomial control.

        Replace the original control input :math:`u` with a polynomial control input :math:`v_{\\text{poly}} = \\sum_{i=0}^d u_i t^i`
        New controls are :math:`u_0, \\dots, u_d`.

        NOTE: bounds on controls are not changed in this function.

        :param degree: degree of the polynomial control
        :type degree: int
        """
        if self.u is None:
            raise Exception('model.u must be defined')
        if self.nu_original is not None:
            raise Exception('model.u has already been augmented')

        casadi_symbol = self.get_casadi_symbol()

        # add time to model
        if self.t == []:
            self.t = casadi_symbol('t')

        t = self.t

        u_old = self.u
        nu_original = casadi_length(self.u)

        u_coeff = casadi_symbol('u_coeff', (degree+1) * nu_original)
        u_new = np.zeros((nu_original, 1))
        for i in range(degree+1):
            u_new += t ** i * u_coeff[i*nu_original:(i+1)*nu_original]

        evaluate_polynomial_u_fun = ca.Function("evaluate_polynomial_u", [u_coeff, t], [u_new])

        self.substitute(u_old, u_new)

        self.u = u_coeff
        self.nu_original = nu_original
        self.p = ca.vertcat(self.p)

        # update name
        self.name = self.name + f"_d{degree}"

        return evaluate_polynomial_u_fun
