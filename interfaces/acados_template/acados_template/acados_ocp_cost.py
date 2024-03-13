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

class AcadosOcpCost:
    """
    Class containing the numerical data of the cost:

    NOTE: all cost terms, except for the terminal one are weighted with the corresponding time step.
    This means given the time steps are :math:`\Delta t_0,..., \Delta t_N`, the total cost is given by:
    :math:`c_\\text{total} = \Delta t_0 \cdot c_0(x_0, u_0, p_0, z_0) + ... + \Delta t_{N-1} \cdot c_{N-1}(x_0, u_0, p_0, z_0) + c_N(x_N, p_N)`.

    This means the Lagrange cost term is given in continuous time, this makes up for a seeminglessly OCP discretization with a nonuniform time grid.

    In case of LINEAR_LS:
    stage cost is
    :math:`l(x,u,z) = 0.5 \cdot || V_x \, x + V_u \, u + V_z \, z - y_\\text{ref}||^2_W`,
    terminal cost is
    :math:`m(x) = 0.5 \cdot || V^e_x \, x - y_\\text{ref}^e||^2_{W^e}`

    In case of NONLINEAR_LS:
    stage cost is
    :math:`l(x,u,z,p) = 0.5 \cdot || y(x,u,z,p) - y_\\text{ref}||^2_W`,
    terminal cost is
    :math:`m(x,p) = 0.5 \cdot || y^e(x,p) - y_\\text{ref}^e||^2_{W^e}`

    In case of CONVEX_OVER_NONLINEAR:
    stage cost is
    :math:`l(x,u,p) = \psi(y(x,u,p) - y_\\text{ref}, p)`,
    terminal cost is
    :math:`m(x, p) = \psi^e (y^e(x,p) - y_\\text{ref}^e, p)`
    """
    def __init__(self):
        # initial stage
        self.__cost_type_0 = None
        self.__W_0 = None
        self.__Vx_0 = None
        self.__Vu_0 = None
        self.__Vz_0 = None
        self.__yref_0 = None
        self.__cost_ext_fun_type_0 = 'casadi'
        self.__Zl_0 = None
        self.__Zu_0 = None
        self.__zl_0 = None
        self.__zu_0 = None
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
        self.__cost_ext_fun_type = 'casadi'
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
        """:math:`y_\\text{ref}^0` - reference at initial shooting node (0).
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
        if isinstance(yref_0, np.ndarray) and len(yref_0.shape) == 1:
            self.__yref_0 = yref_0
        else:
            raise Exception('Invalid yref_0 value, expected 1-dimensional numpy array.')

    @W_0.setter
    def W_0(self, W_0):
        if isinstance(W_0, np.ndarray) and len(W_0.shape) == 2:
            self.__W_0 = W_0
        else:
            raise Exception('Invalid cost W_0 value. ' \
                + 'Should be 2 dimensional numpy array.')

    @Vx_0.setter
    def Vx_0(self, Vx_0):
        if isinstance(Vx_0, np.ndarray) and len(Vx_0.shape) == 2:
            self.__Vx_0 = Vx_0
        else:
            raise Exception('Invalid cost Vx_0 value. ' \
                + 'Should be 2 dimensional numpy array.')

    @Vu_0.setter
    def Vu_0(self, Vu_0):
        if isinstance(Vu_0, np.ndarray) and len(Vu_0.shape) == 2:
            self.__Vu_0 = Vu_0
        else:
            raise Exception('Invalid cost Vu_0 value. ' \
                + 'Should be 2 dimensional numpy array.')

    @Vz_0.setter
    def Vz_0(self, Vz_0):
        if isinstance(Vz_0, np.ndarray) and len(Vz_0.shape) == 2:
            self.__Vz_0 = Vz_0
        else:
            raise Exception('Invalid cost Vz_0 value. ' \
                + 'Should be 2 dimensional numpy array.')

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
        """:math:`y_\\text{ref}` - reference at intermediate shooting nodes (1 to N-1).
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
        if isinstance(W, np.ndarray) and len(W.shape) == 2:
            self.__W = W
        else:
            raise Exception('Invalid cost W value. ' \
                + 'Should be 2 dimensional numpy array.')


    @Vx.setter
    def Vx(self, Vx):
        if isinstance(Vx, np.ndarray) and len(Vx.shape) == 2:
            self.__Vx = Vx
        else:
            raise Exception('Invalid cost Vx value. ' \
                + 'Should be 2 dimensional numpy array.')

    @Vu.setter
    def Vu(self, Vu):
        if isinstance(Vu, np.ndarray) and len(Vu.shape) == 2:
            self.__Vu = Vu
        else:
            raise Exception('Invalid cost Vu value. ' \
                + 'Should be 2 dimensional numpy array.')

    @Vz.setter
    def Vz(self, Vz):
        if isinstance(Vz, np.ndarray) and len(Vz.shape) == 2:
            self.__Vz = Vz
        else:
            raise Exception('Invalid cost Vz value. ' \
                + 'Should be 2 dimensional numpy array.')

    @yref.setter
    def yref(self, yref):
        if isinstance(yref, np.ndarray) and len(yref.shape) == 1:
            self.__yref = yref
        else:
            raise Exception('Invalid yref value, expected 1-dimensional numpy array.')

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
        """:math:`y_\\text{ref}^e` - cost reference at terminal shooting node (N).
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
        if isinstance(W_e, np.ndarray) and len(W_e.shape) == 2:
            self.__W_e = W_e
        else:
            raise Exception('Invalid cost W_e value. ' \
                + 'Should be 2 dimensional numpy array.')

    @Vx_e.setter
    def Vx_e(self, Vx_e):
        if isinstance(Vx_e, np.ndarray) and len(Vx_e.shape) == 2:
            self.__Vx_e = Vx_e
        else:
            raise Exception('Invalid cost Vx_e value. ' \
                + 'Should be 2 dimensional numpy array.')

    @yref_e.setter
    def yref_e(self, yref_e):
        if isinstance(yref_e, np.ndarray) and len(yref_e.shape) == 1:
            self.__yref_e = yref_e
        else:
            raise Exception('Invalid yref_e value, expected 1-dimensional numpy array.')

    @Zl_e.setter
    def Zl_e(self, Zl_e):
        if isinstance(Zl_e, np.ndarray):
            self.__Zl_e = Zl_e
        else:
            raise Exception('Invalid Zl_e value, expected numpy array.')

    @Zu_e.setter
    def Zu_e(self, Zu_e):
        if isinstance(Zu_e, np.ndarray):
            self.__Zu_e = Zu_e
        else:
            raise Exception('Invalid Zu_e value, expected numpy array.')

    @zl_e.setter
    def zl_e(self, zl_e):
        if isinstance(zl_e, np.ndarray):
            self.__zl_e = zl_e
        else:
            raise Exception('Invalid zl_e value, expected numpy array.')

    @zu_e.setter
    def zu_e(self, zu_e):
        if isinstance(zu_e, np.ndarray):
            self.__zu_e = zu_e
        else:
            raise Exception('Invalid zu_e value, expected numpy array.')

    @Zl_0.setter
    def Zl_0(self, Zl_0):
        if isinstance(Zl_0, np.ndarray):
            self.__Zl_0 = Zl_0
        else:
            raise Exception('Invalid Zl_0 value, expected numpy array.')

    @Zu_0.setter
    def Zu_0(self, Zu_0):
        if isinstance(Zu_0, np.ndarray):
            self.__Zu_0 = Zu_0
        else:
            raise Exception('Invalid Zu_0 value, expected numpy array.')

    @zl_0.setter
    def zl_0(self, zl_0):
        if isinstance(zl_0, np.ndarray):
            self.__zl_0 = zl_0
        else:
            raise Exception('Invalid zl_0 value, expected numpy array.')

    @zu_0.setter
    def zu_0(self, zu_0):
        if isinstance(zu_0, np.ndarray):
            self.__zu_0 = zu_0
        else:
            raise Exception('Invalid zu_0 value, expected numpy array.')

    @cost_ext_fun_type_e.setter
    def cost_ext_fun_type_e(self, cost_ext_fun_type_e):
        if cost_ext_fun_type_e in ['casadi', 'generic']:
            self.__cost_ext_fun_type_e = cost_ext_fun_type_e
        else:
            raise Exception("Invalid cost_ext_fun_type_e value, expected one in ['casadi', 'generic'].")

    def set(self, attr, value):
        setattr(self, attr, value)
