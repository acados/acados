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

import numpy as np
from .utils import check_if_nparray_and_flatten, check_if_2d_nparray

class AcadosOcpCost:
    r"""
    Class containing the numerical data of the cost:

    NOTE: By default, the Lagrange cost term provided in continuous time is internally integrated using the explicit Euler method, cost_discretization = 'EULER',
    which allows for a seamless OCP discretization with a nonuniform time grid.
    This means that all cost terms, except for the terminal one, are weighted with the corresponding time step.
    :math:`c_\text{total} = \Delta t_0 \cdot l_0(x_0, u_0, z_0, p_0) + ... + \Delta t_{N-1} \cdot l_{N-1}(x_{N-1}, u_{N-1}, z_{N-1}, p_{N-1}) + l_N(x_N, p_N)`.

    If a nonlinear least-squares or convex-over-nonlinear cost is used, the cost can also be integrated using the same integration scheme,
    which is used for the dynamics, cost_discretization = 'INTEGRATOR'.

    In case of LINEAR_LS:
    stage cost is
    :math:`l(x,u,z) = 0.5 \cdot || V_x \, x + V_u \, u + V_z \, z - y_\text{ref}||^2_W`,
    terminal cost is
    :math:`m(x) = 0.5 \cdot || V^e_x \, x - y_\text{ref}^e||^2_{W^e}`

    In case of NONLINEAR_LS:
    stage cost is
    :math:`l(x,u,z,t,p) = 0.5 \cdot || y(x,u,z,t,p) - y_\text{ref}||^2_W`,
    terminal cost is
    :math:`m(x,p) = 0.5 \cdot || y^e(x,p) - y_\text{ref}^e||^2_{W^e}`

    In case of CONVEX_OVER_NONLINEAR:
    stage cost is
    :math:`l(x,u,z,t,p) = \psi(y(x,u,z,t,p) - y_\text{ref}, t, p)`,
    terminal cost is
    :math:`m(x, p) = \psi^e (y^e(x,p) - y_\text{ref}^e, p)`
    """
    def __init__(self):
        # initial stage
        self.__cost_type_0 = None
        self.__W_0 = None
        self.__Vx_0 = None
        self.__Vu_0 = None
        self.__Vz_0 = None
        self.__yref_0 = None
        self.__Zl_0 = None
        self.__Zu_0 = None
        self.__zl_0 = None
        self.__zu_0 = None
        self.__cost_ext_fun_type_0 = 'casadi'
        self.__cost_source_ext_cost_0 = None # TODO add property, only required for generic
        self.__cost_function_ext_cost_0 = None # TODO add property, only required for generic

        # Lagrange term
        self.__cost_type   = 'LINEAR_LS'  # cost type
        self.__W           = np.zeros((0,0))
        self.__Vx          = np.zeros((0,0))
        self.__Vu          = np.zeros((0,0))
        self.__Vz          = np.zeros((0,0))
        self.__yref        = np.array([])
        self.__Zl          = np.array([])
        self.__Zu          = np.array([])
        self.__zl          = np.array([])
        self.__zu          = np.array([])

        # TODO: check how generic works in templates ?!
        self.__cost_ext_fun_type = 'casadi'
        self.__cost_source_ext_cost = None # TODO add property, only required for generic
        self.__cost_function_ext_cost = None # TODO add property, only required for generic

        # Mayer term
        self.__cost_type_e = 'LINEAR_LS'
        self.__W_e         = np.zeros((0,0))
        self.__Vx_e        = np.zeros((0,0))
        self.__yref_e      = np.array([])
        self.__Zl_e        = np.array([])
        self.__Zu_e        = np.array([])
        self.__zl_e        = np.array([])
        self.__zu_e        = np.array([])
        self.__cost_ext_fun_type_e = 'casadi'
        self.__cost_source_ext_cost_e = None # TODO add property, only required for generic
        self.__cost_function_ext_cost_e = None # TODO add property, only required for generic


    # initial stage
    @property
    def cost_type_0(self):
        """Cost type at initial shooting node (0)
        -- string in {EXTERNAL, LINEAR_LS, NONLINEAR_LS, CONVEX_OVER_NONLINEAR} or :code:`None`.
        Default: :code:`None`.

            .. note:: Cost at initial stage is the same as for intermediate shooting nodes if not set differently explicitly.

            .. note:: If :py:attr:`cost_type_0` is set to :code:`None` values in :py:attr:`W_0`, :py:attr:`Vx_0`, :py:attr:`Vu_0`, :py:attr:`Vz_0` and :py:attr:`yref_0` are ignored (set to :code:`None`).
        """
        return self.__cost_type_0

    @property
    def W_0(self):
        """:math:`W_0` - weight matrix at initial shooting node (0).
        Default: :code:`None`.
        """
        return self.__W_0

    @property
    def Vx_0(self):
        """:math:`V_x^0` - x matrix coefficient at initial shooting node (0).
        Default: :code:`None`.
        """
        return self.__Vx_0

    @property
    def Vu_0(self):
        """:math:`V_u^0` - u matrix coefficient at initial shooting node (0).
        Default: :code:`None`.
        """
        return self.__Vu_0

    @property
    def Vz_0(self):
        """:math:`V_z^0` - z matrix coefficient at initial shooting node (0).
        Default: :code:`None`.
        """
        return self.__Vz_0

    @property
    def yref_0(self):
        r""":math:`y_\text{ref}^0` - reference at initial shooting node (0).
        Default: :code:`None`.
        """
        return self.__yref_0

    @property
    def cost_ext_fun_type_0(self):
        """Type of external function for cost at initial shooting node (0)
        -- string in {casadi, generic} or :code:`None`
        Default: :code:'casadi'.

            .. note:: Cost at initial stage is the same as for intermediate shooting nodes if not set differently explicitly.
        """
        return self.__cost_ext_fun_type_0

    @yref_0.setter
    def yref_0(self, yref_0):
        yref_0 = check_if_nparray_and_flatten(yref_0, "yref_0")
        self.__yref_0 = yref_0

    @W_0.setter
    def W_0(self, W_0):
        check_if_2d_nparray(W_0, "W_0")
        self.__W_0 = W_0

    @Vx_0.setter
    def Vx_0(self, Vx_0):
        check_if_2d_nparray(Vx_0, "Vx_0")
        self.__Vx_0 = Vx_0

    @Vu_0.setter
    def Vu_0(self, Vu_0):
        check_if_2d_nparray(Vu_0, "Vu_0")
        self.__Vu_0 = Vu_0

    @Vz_0.setter
    def Vz_0(self, Vz_0):
        check_if_2d_nparray(Vz_0, "Vz_0")
        self.__Vz_0 = Vz_0

    @cost_ext_fun_type_0.setter
    def cost_ext_fun_type_0(self, cost_ext_fun_type_0):
        if cost_ext_fun_type_0 in ['casadi', 'generic']:
            self.__cost_ext_fun_type_0 = cost_ext_fun_type_0
        else:
            raise Exception('Invalid cost_ext_fun_type_0 value, expected numpy array.')

    # Lagrange term
    @property
    def cost_type(self):
        """
        Cost type at intermediate shooting nodes (1 to N-1)
        -- string in {EXTERNAL, LINEAR_LS, NONLINEAR_LS, CONVEX_OVER_NONLINEAR}.
        Default: 'LINEAR_LS'.
        """
        return self.__cost_type

    @property
    def W(self):
        """:math:`W` - weight matrix at intermediate shooting nodes (1 to N-1).
        Default: :code:`np.zeros((0,0))`.
        """
        return self.__W

    @property
    def Vx(self):
        """:math:`V_x` - x matrix coefficient at intermediate shooting nodes (1 to N-1).
        Default: :code:`np.zeros((0,0))`.
        """
        return self.__Vx

    @property
    def Vu(self):
        """:math:`V_u` - u matrix coefficient at intermediate shooting nodes (1 to N-1).
        Default: :code:`np.zeros((0,0))`.
        """
        return self.__Vu

    @property
    def Vz(self):
        """:math:`V_z` - z matrix coefficient at intermediate shooting nodes (1 to N-1).
        Default: :code:`np.zeros((0,0))`.
        """
        return self.__Vz

    @property
    def yref(self):
        r""":math:`y_\text{ref}` - reference at intermediate shooting nodes (1 to N-1).
        Default: :code:`np.array([])`.
        """
        return self.__yref

    @property
    def Zl(self):
        """:math:`Z_l` - diagonal of Hessian wrt lower slack at intermediate shooting nodes (0 to N-1).
        Default: :code:`np.array([])`.
        """
        return self.__Zl

    @property
    def Zu(self):
        """:math:`Z_u` - diagonal of Hessian wrt upper slack at intermediate shooting nodes (0 to N-1).
        Default: :code:`np.array([])`.
        """
        return self.__Zu

    @property
    def zl(self):
        """:math:`z_l` - gradient wrt lower slack at intermediate shooting nodes (0 to N-1).
        Default: :code:`np.array([])`.
        """
        return self.__zl

    @property
    def zu(self):
        """:math:`z_u` - gradient wrt upper slack at intermediate shooting nodes (0 to N-1).
        Default: :code:`np.array([])`.
        """
        return self.__zu

    @property
    def cost_ext_fun_type(self):
        """Type of external function for cost at intermediate shooting nodes (1 to N-1).
        -- string in {casadi, generic}
        Default: :code:'casadi'.
        """
        return self.__cost_ext_fun_type

    @cost_type.setter
    def cost_type(self, cost_type):
        cost_types = ('LINEAR_LS', 'NONLINEAR_LS', 'EXTERNAL', 'CONVEX_OVER_NONLINEAR')
        if cost_type in cost_types:
            self.__cost_type = cost_type
        else:
            raise Exception('Invalid cost_type value.')

    @cost_type_0.setter
    def cost_type_0(self, cost_type_0):
        cost_types = ('LINEAR_LS', 'NONLINEAR_LS', 'EXTERNAL', 'CONVEX_OVER_NONLINEAR')
        if cost_type_0 in cost_types:
            self.__cost_type_0 = cost_type_0
        else:
            raise Exception('Invalid cost_type_0 value.')

    @W.setter
    def W(self, W):
        check_if_2d_nparray(W, "W")
        self.__W = W


    @Vx.setter
    def Vx(self, Vx):
        check_if_2d_nparray(Vx, "Vx")
        self.__Vx = Vx

    @Vu.setter
    def Vu(self, Vu):
        check_if_2d_nparray(Vu, "Vu")
        self.__Vu = Vu

    @Vz.setter
    def Vz(self, Vz):
        check_if_2d_nparray(Vz, "Vz")
        self.__Vz = Vz

    @yref.setter
    def yref(self, yref):
        yref = check_if_nparray_and_flatten(yref, "yref")
        self.__yref = yref

    @Zl.setter
    def Zl(self, Zl):
        if isinstance(Zl, np.ndarray):
            self.__Zl = Zl
        else:
            raise Exception('Invalid Zl value, expected numpy array.')

    @Zu.setter
    def Zu(self, Zu):
        if isinstance(Zu, np.ndarray):
            self.__Zu = Zu
        else:
            raise Exception('Invalid Zu value, expected numpy array.')

    @zl.setter
    def zl(self, zl):
        if isinstance(zl, np.ndarray):
            self.__zl = zl
        else:
            raise Exception('Invalid zl value, expected numpy array.')

    @zu.setter
    def zu(self, zu):
        if isinstance(zu, np.ndarray):
            self.__zu = zu
        else:
            raise Exception('Invalid zu value, expected numpy array.')

    @cost_ext_fun_type.setter
    def cost_ext_fun_type(self, cost_ext_fun_type):
        if cost_ext_fun_type in ['casadi', 'generic']:
            self.__cost_ext_fun_type = cost_ext_fun_type
        else:
            raise Exception("Invalid cost_ext_fun_type value, expected one in ['casadi', 'generic'].")

    # Mayer term
    @property
    def cost_type_e(self):
        """
        Cost type at terminal shooting node (N)
        -- string in {EXTERNAL, LINEAR_LS, NONLINEAR_LS, CONVEX_OVER_NONLINEAR}.
        Default: 'LINEAR_LS'.
        """
        return self.__cost_type_e

    @property
    def W_e(self):
        """:math:`W_e` - weight matrix at terminal shooting node (N).
        Default: :code:`np.zeros((0,0))`.
        """
        return self.__W_e

    @property
    def Vx_e(self):
        """:math:`V_x^e` - x matrix coefficient for cost at terminal shooting node (N).
        Default: :code:`np.zeros((0,0))`.
        """
        return self.__Vx_e

    @property
    def yref_e(self):
        r""":math:`y_\text{ref}^e` - cost reference at terminal shooting node (N).
        Default: :code:`np.array([])`.
        """
        return self.__yref_e

    @property
    def Zl_e(self):
        """:math:`Z_l^e` - diagonal of Hessian wrt lower slack at terminal shooting node (N).
        Default: :code:`np.array([])`.
        """
        return self.__Zl_e

    @property
    def Zu_e(self):
        """:math:`Z_u^e` - diagonal of Hessian wrt upper slack at terminal shooting node (N).
        Default: :code:`np.array([])`.
        """
        return self.__Zu_e

    @property
    def zl_e(self):
        """:math:`z_l^e` - gradient wrt lower slack at terminal shooting node (N).
        Default: :code:`np.array([])`.
        """
        return self.__zl_e

    @property
    def zu_e(self):
        """:math:`z_u^e` - gradient wrt upper slack at terminal shooting node (N).
        Default: :code:`np.array([])`.
        """
        return self.__zu_e


    @property
    def Zl_0(self):
        """:math:`Z_l^0` - diagonal of Hessian wrt lower slack at initial node 0.
        Default: :code:`np.array([])`.
        """
        return self.__Zl_0

    @property
    def Zu_0(self):
        """:math:`Z_u^0` - diagonal of Hessian wrt upper slack at initial node 0.
        Default: :code:`np.array([])`.
        """
        return self.__Zu_0

    @property
    def zl_0(self):
        """:math:`z_l^0` - gradient wrt lower slack at initial node 0.
        Default: :code:`np.array([])`.
        """
        return self.__zl_0

    @property
    def zu_0(self):
        """:math:`z_u^0` - gradient wrt upper slack at initial node 0.
        Default: :code:`np.array([])`.
        """
        return self.__zu_0

    @property
    def cost_ext_fun_type_e(self):
        """Type of external function for cost at terminal shooting node (N).
        -- string in {casadi, generic}
        Default: :code:'casadi'.
        """
        return self.__cost_ext_fun_type_e

    @cost_type_e.setter
    def cost_type_e(self, cost_type_e):
        cost_types = ('LINEAR_LS', 'NONLINEAR_LS', 'EXTERNAL', 'CONVEX_OVER_NONLINEAR')
        if cost_type_e in cost_types:
            self.__cost_type_e = cost_type_e
        else:
            raise Exception('Invalid cost_type_e value.')

    @W_e.setter
    def W_e(self, W_e):
        check_if_2d_nparray(W_e, "W_e")
        self.__W_e = W_e

    @Vx_e.setter
    def Vx_e(self, Vx_e):
        check_if_2d_nparray(Vx_e, "Vx_e")
        self.__Vx_e = Vx_e

    @yref_e.setter
    def yref_e(self, yref_e):
        yref_e = check_if_nparray_and_flatten(yref_e, "yref_e")
        self.__yref_e = yref_e

    @Zl_e.setter
    def Zl_e(self, Zl_e):
        Zl_e = check_if_nparray_and_flatten(Zl_e, "Zl_e")
        self.__Zl_e = Zl_e

    @Zu_e.setter
    def Zu_e(self, Zu_e):
        Zu_e = check_if_nparray_and_flatten(Zu_e, "Zu_e")
        self.__Zu_e = Zu_e

    @zl_e.setter
    def zl_e(self, zl_e):
        zl_e = check_if_nparray_and_flatten(zl_e, "zl_e")
        self.__zl_e = zl_e

    @zu_e.setter
    def zu_e(self, zu_e):
        zu_e = check_if_nparray_and_flatten(zu_e, "zu_e")
        self.__zu_e = zu_e

    @Zl_0.setter
    def Zl_0(self, Zl_0):
        Zl_0 = check_if_nparray_and_flatten(Zl_0, "Zl_0")
        self.__Zl_0 = Zl_0

    @Zu_0.setter
    def Zu_0(self, Zu_0):
        Zu_0 = check_if_nparray_and_flatten(Zu_0, "Zu_0")
        self.__Zu_0 = Zu_0

    @zl_0.setter
    def zl_0(self, zl_0):
        zl_0 = check_if_nparray_and_flatten(zl_0, "zl_0")
        self.__zl_0 = zl_0

    @zu_0.setter
    def zu_0(self, zu_0):
        zu_0 = check_if_nparray_and_flatten(zu_0, "zu_0")
        self.__zu_0 = zu_0

    @cost_ext_fun_type_e.setter
    def cost_ext_fun_type_e(self, cost_ext_fun_type_e):
        if cost_ext_fun_type_e in ['casadi', 'generic']:
            self.__cost_ext_fun_type_e = cost_ext_fun_type_e
        else:
            raise Exception("Invalid cost_ext_fun_type_e value, expected one in ['casadi', 'generic'].")

    def set(self, attr, value):
        setattr(self, attr, value)
