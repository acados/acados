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

from .utils import is_empty, casadi_length
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
        self.__name = None
        self.__x = []
        self.__xdot = []
        self.__u = []
        self.__z = []
        self.__p = []
        self.__t = []
        self.__p_global = []

        ## dynamics
        self.__f_impl_expr = []
        self.__f_expl_expr = []
        self.__disc_dyn_expr = []

        self.__dyn_ext_fun_type = 'casadi'
        self.__dyn_generic_source = None
        self.__dyn_disc_fun_jac_hess = None
        self.__dyn_disc_fun_jac = None
        self.__dyn_disc_fun = None

        self.__dyn_impl_dae_fun_jac = None
        self.__dyn_impl_dae_jac = None
        self.__dyn_impl_dae_fun = None

        # for GNSF models
        self.__gnsf_nontrivial_f_LO = 1
        self.__gnsf_purely_linear = 0

        ### for OCP only.
        # NOTE: These could be moved to cost / constraints

        # constraints at initial stage
        self.__con_h_expr_0 = []
        self.__con_phi_expr_0 = []
        self.__con_r_expr_0 = []
        self.__con_r_in_phi_0 = []


        # path constraints
        # BGH(default): lh <= h(x, u) <= uh
        self.__con_h_expr = []
        # BGP(convex over nonlinear): lphi <= phi(r(x, u)) <= uphi
        self.__con_phi_expr = []
        self.__con_r_expr = []
        self.__con_r_in_phi = []

        # terminal
        self.__con_h_expr_e = []
        self.__con_phi_expr_e = []
        self.__con_r_expr_e = []
        self.__con_r_in_phi_e = []

        # cost
        self.__cost_y_expr = []
        self.__cost_y_expr_e = []
        self.__cost_y_expr_0 = []
        self.__cost_expr_ext_cost = []
        self.__cost_expr_ext_cost_e = []
        self.__cost_expr_ext_cost_0 = []
        self.__cost_expr_ext_cost_custom_hess = []
        self.__cost_expr_ext_cost_custom_hess_e = []
        self.__cost_expr_ext_cost_custom_hess_0 = []

        ## CONVEX_OVER_NONLINEAR convex-over-nonlinear cost: psi(y(x, u, p) - y_ref; p)
        self.__cost_psi_expr_0 = []
        self.__cost_psi_expr = []
        self.__cost_psi_expr_e = []
        self.__cost_r_in_psi_expr_0 = []
        self.__cost_r_in_psi_expr = []
        self.__cost_r_in_psi_expr_e = []
        self.__cost_conl_custom_outer_hess_0 = []
        self.__cost_conl_custom_outer_hess = []
        self.__cost_conl_custom_outer_hess_e = []
        self.__nu_original = None # TODO: remove? only used by benchmark

        self.__t0 = None # TODO: remove? only used by benchmark
        self.__x_labels = None
        self.__u_labels = None
        self.__t_label = "t"

    @property
    def name(self):
        """The model name is used for code generation. Type: string.
        Default: :code:`None`
        """
        return self.__name

    @name.setter
    def name(self, name):
        self.__name = name

    @property
    def x(self):
        """CasADi variable describing the state of the system;
        Default: :code:`[]`
        """
        return self.__x

    @x.setter
    def x(self, x):
        self.__x = x

    @property
    def xdot(self):
        """CasADi variable describing the derivative of the state wrt time;
        Default: :code:`[]`
        """
        return self.__xdot

    @xdot.setter
    def xdot(self, xdot):
        self.__xdot = xdot

    @property
    def u(self):
        """CasADi variable describing the derivative of the state wrt time;
        Default: :code:`[]`
        """
        return self.__u

    @u.setter
    def u(self, u):
        self.__u = u

    @property
    def z(self):
        """CasADi variable describing the algebraic variables of the DAE;
        Default: :code:`[]`
        """
        return self.__z

    @z.setter
    def z(self, z):
        self.__z = z

    @property
    def p(self):
        """CasADi variable describing stage-wise parameters;
        Default: :code:`[]`
        """
        return self.__p

    @p.setter
    def p(self, p):
        self.__p = p

    @property
    def t(self):
        """CasADi variable representing time t in functions;
        Default: :code:`[]`
        NOTE:
        - For integrators, the start time has to be explicitly set via :py:attr:`acados_template.AcadosSimSolver.set`('t0').
        - For OCPs, the start time is set to 0 on each stage.
        The time dependency can be used within cost formulations and is relevant when cost integration is used.
        Start times of shooting intervals can be added using parameters.
        """
        return self.__t

    @t.setter
    def t(self, t):
        self.__t = t

    @property
    def p_global(self):
        """
        CasADi variable representing global parameters.
        This feature can be used to precompute expensive terms which only depend on these parameters, e.g. spline coefficients, when p_global are underlying data points.
        Only supported for OCP solvers.
        Updating these parameters can be done using :py:attr:`acados_template.acados_ocp_solver.AcadosOcpSolver.set_p_global_and_precompute_dependencies(values)`.
        NOTE: this is only supported with CasADi beta release https://github.com/casadi/casadi/releases/tag/nightly-se
        Default: :code:`[]`
        """
        return self.__p_global

    @p_global.setter
    def p_global(self, p_global):
        self.__p_global = p_global

    @property
    def f_impl_expr(self):
        r"""
        CasADi expression for the implicit dynamics :math:`f_\text{impl}(\dot{x}, x, u, z, p) = 0`.
        Used if :py:attr:`acados_template.acados_ocp_options.AcadosOcpOptions.integrator_type` == 'IRK'.
        Default: :code:`[]`
        """
        return self.__f_impl_expr

    @f_impl_expr.setter
    def f_impl_expr(self, f_impl_expr):
        self.__f_impl_expr = f_impl_expr

    @property
    def f_expl_expr(self):
        r"""
        CasADi expression for the explicit dynamics :math:`\dot{x} = f_\text{expl}(x, u, p)`.
        Used if :py:attr:`acados_template.acados_ocp_options.AcadosOcpOptions.integrator_type` == 'ERK'.
        Default: :code:`[]`
        """
        return self.__f_expl_expr

    @f_expl_expr.setter
    def f_expl_expr(self, f_expl_expr):
        self.__f_expl_expr = f_expl_expr

    @property
    def disc_dyn_expr(self):
        r"""
        CasADi expression for the discrete dynamics :math:`x_{+} = f_\text{disc}(x, u, p)`.
        Used if :py:attr:`acados_template.acados_ocp_options.AcadosOcpOptions.integrator_type` == 'DISCRETE'.
        Default: :code:`[]`
        """
        return self.__disc_dyn_expr

    @disc_dyn_expr.setter
    def disc_dyn_expr(self, disc_dyn_expr):
        self.__disc_dyn_expr = disc_dyn_expr

    @property
    def dyn_ext_fun_type(self):
        """
        Type of external functions for dynamics module; 'casadi' or 'generic';
        Default: 'casadi'
        """
        return self.__dyn_ext_fun_type

    @dyn_ext_fun_type.setter
    def dyn_ext_fun_type(self, dyn_ext_fun_type):
        self.__dyn_ext_fun_type = dyn_ext_fun_type

    @property
    def dyn_generic_source(self):
        """
        Name of source file for discrete dynamics, only relevant if :code:`dyn_ext_fun_type` is :code:`'generic'`;
        Default: :code:`None`
        """
        return self.__dyn_generic_source

    @dyn_generic_source.setter
    def dyn_generic_source(self, dyn_generic_source):
        self.__dyn_generic_source = dyn_generic_source

    @property
    def dyn_disc_fun_jac_hess(self):
        """
        Name of function discrete dynamics + jacobian and hessian, only relevant if :code:`dyn_ext_fun_type` is :code:`'generic'`;
        Default: :code:`None`
        """
        return self.__dyn_disc_fun_jac_hess

    @dyn_disc_fun_jac_hess.setter
    def dyn_disc_fun_jac_hess(self, dyn_disc_fun_jac_hess):
        self.__dyn_disc_fun_jac_hess = dyn_disc_fun_jac_hess


    @property
    def dyn_disc_fun_jac(self):
        """
        Name of function discrete dynamics + jacobian, only relevant if :code:`dyn_ext_fun_type` is :code:`'generic'`;
        Default: :code:`None`
        """
        return self.__dyn_disc_fun_jac

    @dyn_disc_fun_jac.setter
    def dyn_disc_fun_jac(self, dyn_disc_fun_jac):
        self.__dyn_disc_fun_jac = dyn_disc_fun_jac



    @property
    def dyn_disc_fun(self):
        """
        Name of function discrete dynamics, only relevant if :code:`dyn_ext_fun_type` is :code:`'generic'`;
        Default: :code:`None`
        """
        return self.__dyn_disc_fun

    @dyn_disc_fun.setter
    def dyn_disc_fun(self, dyn_disc_fun):
        self.__dyn_disc_fun = dyn_disc_fun

    @property
    def dyn_impl_dae_fun_jac(self):
        """
        Name of source files for implicit DAE function value and jacobian, only relevant if :code:`dyn_ext_fun_type` is :code:`'generic'`;
        Default: :code:`None`
        """
        return self.__dyn_impl_dae_fun_jac

    @dyn_impl_dae_fun_jac.setter
    def dyn_impl_dae_fun_jac(self, dyn_impl_dae_fun_jac):
        self.__dyn_impl_dae_fun_jac = dyn_impl_dae_fun_jac

    @property
    def dyn_impl_dae_jac(self):
        """
        Name of source files for implicit DAE jacobian, only relevant if :code:`dyn_ext_fun_type` is :code:`'generic'`;
        Default: :code:`None`
        """
        return self.__dyn_impl_dae_jac

    @dyn_impl_dae_jac.setter
    def dyn_impl_dae_jac(self, dyn_impl_dae_jac):
        self.__dyn_impl_dae_jac = dyn_impl_dae_jac

    @property
    def dyn_impl_dae_fun(self):
        """
        Name of source files for implicit DAE function value, only relevant if :code:`dyn_ext_fun_type` is :code:`'generic'`;
        Default: :code:`None`
        """
        return self.__dyn_impl_dae_fun

    @dyn_impl_dae_fun.setter
    def dyn_impl_dae_fun(self, dyn_impl_dae_fun):
        self.__dyn_impl_dae_fun = dyn_impl_dae_fun

    @property
    def gnsf_nontrivial_f_LO(self):
        """
        GNSF: Flag indicating whether GNSF stucture has nontrivial f.
        """
        return self.__gnsf_nontrivial_f_LO

    @gnsf_nontrivial_f_LO.setter
    def gnsf_nontrivial_f_LO(self, gnsf_nontrivial_f_LO):
        self.__gnsf_nontrivial_f_LO = gnsf_nontrivial_f_LO

    @property
    def gnsf_purely_linear(self):
        """
        GNSF: Flag indicating whether GNSF stucture is purely linear.
        """
        return self.__gnsf_purely_linear

    @gnsf_purely_linear.setter
    def gnsf_purely_linear(self, gnsf_purely_linear):
        self.__gnsf_purely_linear = gnsf_purely_linear

    @property
    def con_h_expr_0(self):
        r"""
        CasADi expression for the initial constraint :math:`h^0`;
        Default: :code:`[]`
        """
        return self.__con_h_expr_0

    @con_h_expr_0.setter
    def con_h_expr_0(self, con_h_expr_0):
        self.__con_h_expr_0 = con_h_expr_0

    @property
    def con_phi_expr_0(self):
        r"""
        CasADi expression for the outer function of the initial constraint :math:`\phi_0(r_0)`;
        Default: :code:`[]`
        """
        return self.__con_phi_expr_0

    @con_phi_expr_0.setter
    def con_phi_expr_0(self, con_phi_expr_0):
        self.__con_phi_expr_0 = con_phi_expr_0

    @property
    def con_r_expr_0(self):
        r"""
        CasADi expression for the inner function of the initial constraint :math:`\phi_0(r_0)`;
        Default: :code:`[]`
        """
        return self.__con_r_expr_0

    @con_r_expr_0.setter
    def con_r_expr_0(self, con_r_expr_0):
        self.__con_r_expr_0 = con_r_expr_0

    @property
    def con_r_in_phi_0(self):
        r"""
        CasADi variable defining the input to the outer function of the initial constraint :math:`\phi_0(r_0)`;
        Default: :code:`[]`
        """
        return self.__con_r_in_phi_0

    @con_r_in_phi_0.setter
    def con_r_in_phi_0(self, con_r_in_phi_0):
        self.__con_r_in_phi_0 = con_r_in_phi_0

    @property
    def con_h_expr(self):
        """
        CasADi expression for the intermediate constraint :math:`h`;
        Default: :code:`[]`
        """
        return self.__con_h_expr

    @con_h_expr.setter
    def con_h_expr(self, con_h_expr):
        self.__con_h_expr = con_h_expr

    @property
    def con_phi_expr(self):
        r"""
        CasADi expression for the outer function of the intermediate constraint :math:`\phi(r)`;
        Default: :code:`[]`
        """
        return self.__con_phi_expr

    @con_phi_expr.setter
    def con_phi_expr(self, con_phi_expr):
        self.__con_phi_expr = con_phi_expr

    @property
    def con_r_expr(self):
        r"""
        CasADi expression for the inner function of the intermediate constraint :math:`\phi(r)`;
        Default: :code:`[]`
        """
        return self.__con_r_expr

    @con_r_expr.setter
    def con_r_expr(self, con_r_expr):
        self.__con_r_expr = con_r_expr

    @property
    def con_r_in_phi(self):
        r"""
        CasADi variable defining the input to the outer function of the intermediate constraint :math:`\phi(r)`;
        Default: :code:`[]`
        """
        return self.__con_r_in_phi

    @con_r_in_phi.setter
    def con_r_in_phi(self, con_r_in_phi):
        self.__con_r_in_phi = con_r_in_phi

    @property
    def con_h_expr_e(self):
        """
        CasADi expression for the terminal constraint :math:`h_e`;
        Default: :code:`[]`
        """
        return self.__con_h_expr_e

    @con_h_expr_e.setter
    def con_h_expr_e(self, con_h_expr_e):
        self.__con_h_expr_e = con_h_expr_e

    @property
    def con_phi_expr_e(self):
        r"""
        CasADi expression for the outer function of the terminal constraint :math:`\phi_e(r_e)`;
        Default: :code:`[]`
        """
        return self.__con_phi_expr_e

    @con_phi_expr_e.setter
    def con_phi_expr_e(self, con_phi_expr_e):
        self.__con_phi_expr_e = con_phi_expr_e

    @property
    def con_r_expr_e(self):
        r"""
        CasADi expression for the inner function of the terminal constraint :math:`\phi_e(r_e)`;
        Default: :code:`[]`
        """
        return self.__con_r_expr_e

    @con_r_expr_e.setter
    def con_r_expr_e(self, con_r_expr_e):
        self.__con_r_expr_e = con_r_expr_e

    @property
    def con_r_in_phi_e(self):
        r"""
        CasADi variable defining the input to the outer function of the terminal constraint :math:`\phi_e(r_e)`;
        Default: :code:`[]`
        """
        return self.__con_r_in_phi_e

    @con_r_in_phi_e.setter
    def con_r_in_phi_e(self, con_r_in_phi_e):
        self.__con_r_in_phi_e = con_r_in_phi_e

    @property
    def cost_y_expr(self):
        """CasADi expression for nonlinear least squares;
        Default: :code:`[]`
        """
        return self.__cost_y_expr

    @cost_y_expr.setter
    def cost_y_expr(self, cost_y_expr):
        self.__cost_y_expr = cost_y_expr


    @property
    def cost_y_expr_e(self):
        """CasADi expression for nonlinear least squares, terminal;
        Default: :code:`[]`
        """
        return self.__cost_y_expr_e

    @cost_y_expr_e.setter
    def cost_y_expr_e(self, cost_y_expr_e):
        self.__cost_y_expr_e = cost_y_expr_e

    @property
    def cost_y_expr_0(self):
        """CasADi expression for nonlinear least squares, initial;
        Default: :code:`[]`
        """
        return self.__cost_y_expr_0

    @cost_y_expr_0.setter
    def cost_y_expr_0(self, cost_y_expr_0):
        self.__cost_y_expr_0 = cost_y_expr_0

    @property
    def cost_expr_ext_cost(self):
        """CasADi expression for external cost;
        Default: :code:`[]`
        """
        return self.__cost_expr_ext_cost

    @cost_expr_ext_cost.setter
    def cost_expr_ext_cost(self, cost_expr_ext_cost):
        self.__cost_expr_ext_cost = cost_expr_ext_cost

    @property
    def cost_expr_ext_cost_e(self):
        """CasADi expression for external cost, terminal;
        Default: :code:`[]`
        """
        return self.__cost_expr_ext_cost_e

    @cost_expr_ext_cost_e.setter
    def cost_expr_ext_cost_e(self, cost_expr_ext_cost_e):
        self.__cost_expr_ext_cost_e = cost_expr_ext_cost_e

    @property
    def cost_expr_ext_cost_0(self):
        """CasADi expression for external cost, initial;
        Default: :code:`[]`
        """
        return self.__cost_expr_ext_cost_0

    @cost_expr_ext_cost_0.setter
    def cost_expr_ext_cost_0(self, cost_expr_ext_cost_0):
        self.__cost_expr_ext_cost_0 = cost_expr_ext_cost_0


    @property
    def cost_expr_ext_cost_custom_hess(self):
        """CasADi expression for custom hessian (only for external cost);
        Default: :code:`[]`
        """
        return self.__cost_expr_ext_cost_custom_hess

    @cost_expr_ext_cost_custom_hess.setter
    def cost_expr_ext_cost_custom_hess(self, cost_expr_ext_cost_custom_hess):
        self.__cost_expr_ext_cost_custom_hess = cost_expr_ext_cost_custom_hess

    @property
    def cost_expr_ext_cost_custom_hess_e(self):
        """CasADi expression for custom hessian (only for external cost), terminal;
        Default: :code:`[]`
        """
        return self.__cost_expr_ext_cost_custom_hess_e

    @cost_expr_ext_cost_custom_hess_e.setter
    def cost_expr_ext_cost_custom_hess_e(self, cost_expr_ext_cost_custom_hess_e):
        self.__cost_expr_ext_cost_custom_hess_e = cost_expr_ext_cost_custom_hess_e

    @property
    def cost_expr_ext_cost_custom_hess_0(self):
        """CasADi expression for custom hessian (only for external cost), initial;
        Default: :code:`[]`
        """
        return self.__cost_expr_ext_cost_custom_hess_0

    @cost_expr_ext_cost_custom_hess_0.setter
    def cost_expr_ext_cost_custom_hess_0(self, cost_expr_ext_cost_custom_hess_0):
        self.__cost_expr_ext_cost_custom_hess_0 = cost_expr_ext_cost_custom_hess_0

    @property
    def cost_psi_expr_0(self):
        r"""
        CasADi expression for the outer loss function :math:`\psi(r - yref, t, p)`, initial;
        Default: :code:`[]`
        Used if :py:attr:`acados_template.acados_ocp_options.AcadosOcpOptions.cost_type_0` is 'CONVEX_OVER_NONLINEAR'.
        """
        return self.__cost_psi_expr_0

    @cost_psi_expr_0.setter
    def cost_psi_expr_0(self, cost_psi_expr_0):
        self.__cost_psi_expr_0 = cost_psi_expr_0


    @property
    def cost_psi_expr(self):
        r"""
        CasADi expression for the outer loss function :math:`\psi(r - yref, t, p)`;
        Default: :code:`[]`
        Used if :py:attr:`acados_template.acados_ocp_options.AcadosOcpOptions.cost_type` is 'CONVEX_OVER_NONLINEAR'.
        """
        return self.__cost_psi_expr

    @cost_psi_expr.setter
    def cost_psi_expr(self, cost_psi_expr):
        self.__cost_psi_expr = cost_psi_expr

    @property
    def cost_psi_expr_e(self):
        r"""
        CasADi expression for the outer loss function :math:`\psi(r - yref, t, p)`, terminal;
        Default: :code:`[]`
        Used if :py:attr:`acados_template.acados_ocp_options.AcadosOcpOptions.cost_type_e` is 'CONVEX_OVER_NONLINEAR'.
        """
        return self.__cost_psi_expr_e

    @cost_psi_expr_e.setter
    def cost_psi_expr_e(self, cost_psi_expr_e):
        self.__cost_psi_expr_e = cost_psi_expr_e

    @property
    def cost_r_in_psi_expr_0(self):
        r"""
        CasADi symbolic input variable for the argument :math:`r` to the outer loss function :math:`\psi(r, t, p)`, initial;
        Default: :code:`[]`
        Used if :py:attr:`acados_template.acados_ocp_options.AcadosOcpOptions.cost_type_0` is 'CONVEX_OVER_NONLINEAR'.
        """
        return self.__cost_r_in_psi_expr_0

    @cost_r_in_psi_expr_0.setter
    def cost_r_in_psi_expr_0(self, cost_r_in_psi_expr_0):
        self.__cost_r_in_psi_expr_0 = cost_r_in_psi_expr_0

    @property
    def cost_r_in_psi_expr(self):
        r"""
        CasADi symbolic input variable for the argument :math:`r` to the outer loss function :math:`\psi(r, t, p)`;
        Default: :code:`[]`
        Used if :py:attr:`acados_template.acados_ocp_options.AcadosOcpOptions.cost_type` is 'CONVEX_OVER_NONLINEAR'.
        """
        return self.__cost_r_in_psi_expr

    @cost_r_in_psi_expr.setter
    def cost_r_in_psi_expr(self, cost_r_in_psi_expr):
        self.__cost_r_in_psi_expr = cost_r_in_psi_expr


    @property
    def cost_r_in_psi_expr_e(self):
        r"""
        CasADi symbolic input variable for the argument :math:`r` to the outer loss function :math:`\psi(r, t, p)`, terminal;
        Default: :code:`[]`
        Used if :py:attr:`acados_template.acados_ocp_options.AcadosOcpOptions.cost_type_e` is 'CONVEX_OVER_NONLINEAR'.
        """
        return self.__cost_r_in_psi_expr_e

    @cost_r_in_psi_expr_e.setter
    def cost_r_in_psi_expr_e(self, cost_r_in_psi_expr_e):
        self.__cost_r_in_psi_expr_e = cost_r_in_psi_expr_e


    @property
    def cost_conl_custom_outer_hess_0(self):
        """
        CasADi expression for the custom hessian of the outer loss function (only for convex-over-nonlinear cost), initial;
        Default: :code:`[]`
        Used if :py:attr:`acados_template.acados_ocp_options.AcadosOcpOptions.cost_type_0` is 'CONVEX_OVER_NONLINEAR'.
        """
        return self.__cost_conl_custom_outer_hess_0

    @cost_conl_custom_outer_hess_0.setter
    def cost_conl_custom_outer_hess_0(self, cost_conl_custom_outer_hess_0):
        self.__cost_conl_custom_outer_hess_0 = cost_conl_custom_outer_hess_0

    @property
    def cost_conl_custom_outer_hess(self):
        """
        CasADi expression for the custom hessian of the outer loss function (only for convex-over-nonlinear cost);
        Default: :code:`[]`
        Used if :py:attr:`acados_template.acados_ocp_options.AcadosOcpOptions.cost_type` is 'CONVEX_OVER_NONLINEAR'.
        """
        return self.__cost_conl_custom_outer_hess

    @cost_conl_custom_outer_hess.setter
    def cost_conl_custom_outer_hess(self, cost_conl_custom_outer_hess):
        self.__cost_conl_custom_outer_hess = cost_conl_custom_outer_hess


    @property
    def cost_conl_custom_outer_hess_e(self):
        """
        CasADi expression for the custom hessian of the outer loss function (only for convex-over-nonlinear cost), terminal;
        Default: :code:`[]`
        Used if :py:attr:`acados_template.acados_ocp_options.AcadosOcpOptions.cost_type_e` is 'CONVEX_OVER_NONLINEAR'.
        """
        return self.__cost_conl_custom_outer_hess_e

    @cost_conl_custom_outer_hess_e.setter
    def cost_conl_custom_outer_hess_e(self, cost_conl_custom_outer_hess_e):
        self.__cost_conl_custom_outer_hess_e = cost_conl_custom_outer_hess_e


    @property
    def nu_original(self):
        """
        Number of original control inputs (before polynomial control augmentation);
        Default: :code:`None`
        """
        return self.__nu_original

    @nu_original.setter
    def nu_original(self, nu_original):
        self.__nu_original = nu_original


    @property
    def t0(self):
        """
        Only relevant for benchmark, TODO
        """
        return self.__t0

    @t0.setter
    def t0(self, t0):
        self.__t0 = t0


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
            return [f"u{i}" for i in range(self.u.size()[0])]
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
            raise TypeError(f"model.x must be casadi.SX or casadi.MX, got {type(self.x)}")

    def get_casadi_zeros(self):
        if isinstance(self.x, MX):
            return MX.zeros
        elif isinstance(self.x, SX):
            return SX.zeros
        else:
            raise TypeError(f"model.x must be casadi.SX or casadi.MX, got {type(self.x)}")


    def make_consistent(self, dims: Union[AcadosOcpDims, AcadosSimDims]) -> None:

        casadi_symbol = self.get_casadi_symbol()

        # nx
        if is_empty(self.x):
            raise ValueError("model.x must be defined")
        else:
            dims.nx = casadi_length(self.x)

        if is_empty(self.xdot):
            self.xdot = casadi_symbol('xdot', dims.nx, 1)
        else:
            if casadi_length(self.xdot) != dims.nx:
                raise ValueError(f"model.xdot must have length nx = {dims.nx}, got {casadi_length(self.xdot)}")

        # nu
        if is_empty(self.u):
            self.u = casadi_symbol('u', 0, 1)
        dims.nu = casadi_length(self.u)

        # nz
        if is_empty(self.z):
            self.z = casadi_symbol('z', 0, 1)
        dims.nz = casadi_length(self.z)

        # np
        if is_empty(self.p):
            self.p = casadi_symbol('p', 0, 1)
        dims.np = casadi_length(self.p)

        # np_global
        if is_empty(self.p_global):
            self.p_global = casadi_symbol('p_global', 0, 1)
        dims.np_global = casadi_length(self.p_global)

        # sanity checks
        for symbol, name in [(self.x, 'x'), (self.xdot, 'xdot'), (self.u, 'u'), (self.z, 'z'), (self.p, 'p'), (self.p_global, 'p_global')]:
            if not isinstance(symbol, (ca.MX, ca.SX)):
                raise TypeError(f"model.{name} must be casadi.MX, casadi.SX got {type(symbol)}")
            if not symbol.is_valid_input():
                raise ValueError(f"model.{name} must be valid CasADi symbol, got {symbol}")

        # p_global
        if not is_empty(self.p_global):
            if isinstance(dims, AcadosSimDims):
                raise NotImplementedError("model.p_global is only supported for OCPs")
            if any(ca.which_depends(self.p_global, self.p)):
                raise ValueError(f"model.p_global must not depend on model.p, got p_global ={self.p_global}, p = {self.p}")

        # model output dimension nx_next: dimension of the next state
        if isinstance(dims, AcadosOcpDims):
            if not is_empty(self.disc_dyn_expr):
                dims.nx_next = casadi_length(self.disc_dyn_expr)
            else:
                dims.nx_next = casadi_length(self.x)

        if not is_empty(self.f_impl_expr):
            if casadi_length(self.f_impl_expr) != (dims.nx + dims.nz):
                raise ValueError(f"model.f_impl_expr must have length nx + nz = {dims.nx} + {dims.nz}, got {casadi_length(self.f_impl_expr)}")
        if not is_empty(self.f_expl_expr):
            if casadi_length(self.f_expl_expr) != dims.nx:
                raise ValueError(f"model.f_expl_expr must have length nx = {dims.nx}, got {casadi_length(self.f_expl_expr)}")

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
        r"""
        Augment the model with polynomial control.

        Replace the original control input :math:`u` with a polynomial control input :math:`v_{\text{poly}} = \sum_{i=0}^d u_i t^i`
        New controls are :math:`u_0, \dots, u_d`.

        NOTE: bounds on controls are not changed in this function.

        :param degree: degree of the polynomial control
        :type degree: int
        """
        if self.u is None:
            raise ValueError('model.u must be defined')
        if self.nu_original is not None:
            raise ValueError('model.u has already been augmented')

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


    def _has_custom_hess(self) -> bool:
        return not (is_empty(self.cost_expr_ext_cost_custom_hess_0) and
                    is_empty(self.cost_expr_ext_cost_custom_hess) and 
                    is_empty(self.cost_expr_ext_cost_custom_hess_e))
