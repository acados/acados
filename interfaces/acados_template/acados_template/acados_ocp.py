# -*- coding: future_fstrings -*-
#
# Copyright 2019 Gianluca Frison, Dimitris Kouzoupis, Robin Verschueren,
# Andrea Zanelli, Niels van Duijkeren, Jonathan Frey, Tommaso Sartor,
# Branimir Novoselnik, Rien Quirynen, Rezart Qelibari, Dang Doan,
# Jonas Koenemann, Yutao Chen, Tobias Schöls, Jonas Schlagenhauf, Moritz Diehl
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
import os
from .acados_model import AcadosModel
from .utils import get_acados_path, J_to_idx, J_to_idx_slack

class AcadosOcpDims:
    """
    class containing the dimensions of the optimal control problem
    """
    def __init__(self):
        self.__nx      = None
        self.__nu      = None
        self.__nz      = 0
        self.__np      = 0
        self.__ny      = 0
        self.__ny_e    = 0
        self.__nr      = 0
        self.__nr_e    = 0
        self.__nh      = 0
        self.__nh_e    = 0
        self.__nphi    = 0
        self.__nphi_e  = 0
        self.__nbx     = 0
        self.__nbx_0   = 0
        self.__nbx_e   = 0
        self.__nbu     = 0
        self.__nsbx    = 0
        self.__nsbx_e  = 0
        self.__nsbu    = 0
        self.__nsh     = 0
        self.__nsh_e   = 0
        self.__nsphi   = 0
        self.__nsphi_e = 0
        self.__ns      = 0
        self.__ns_e    = 0
        self.__ng      = 0
        self.__ng_e    = 0
        self.__nsg     = 0
        self.__nsg_e   = 0
        self.__nbxe_0  = 0
        self.__N       = None


    @property
    def nx(self):
        """:math:`n_x` - number of states"""
        return self.__nx

    @property
    def nz(self):
        """:math:`n_z` - number of algebraic variables"""
        return self.__nz

    @property
    def nu(self):
        """:math:`n_u` - number of inputs"""
        return self.__nu

    @property
    def np(self):
        """:math:`n_p` - number of parameters"""
        return self.__np

    @property
    def ny(self):
        """:math:`n_y` - number of residuals in Lagrange term"""
        return self.__ny

    @property
    def ny_e(self):
        """:math:`n_{y}^e` - number of residuals in Mayer term"""
        return self.__ny_e

    @property
    def nr(self):
        """:math:`n_{\pi}` - dimension of the image of the inner nonlinear function in positive definite constraints"""
        return self.__nr

    @property
    def nr_e(self):
        """:math:`n_{\pi}^e` - dimension of the image of the inner nonlinear function in positive definite constraints"""
        return self.__nr_e

    @property
    def nh(self):
        """:math:`n_h` - number of nonlinear constraints"""
        return self.__nh

    @property
    def nh_e(self):
        """:math:`n_{h}^e` - number of nonlinear constraints at t=T"""
        return self.__nh_e

    @property
    def nphi(self):
        """:math:`n_{\phi}` - number of convex-over-nonlinear constraints"""
        return self.__nphi

    @property
    def nphi_e(self):
        """:math:`n_{\phi}^e` - number of convex-over-nonlinear constraints at t=T"""
        return self.__nphi_e

    @property
    def nbx(self):
        """:math:`n_{b_x}` - number of state bounds"""
        return self.__nbx

    @property
    def nbxe_0(self):
        """:math:`n_{be_{x0}}` - number of state bounds at initial shooting node that are equalities"""
        return self.__nbxe_0

    @property
    def nbx_0(self):
        """:math:`n_{b_{x0}}` - number of state bounds for initial state"""
        return self.__nbx_0

    @property
    def nbx_e(self):
        """:math:`n_{b_x}` - number of state bounds at t=T"""
        return self.__nbx_e

    @property
    def nbu(self):
        """:math:`n_{b_u}` - number of input bounds"""
        return self.__nbu

    @property
    def nsbx(self):
        """:math:`n_{{sb}_x}` - number of soft state bounds"""
        return self.__nsbx

    @property
    def nsbx_e(self):
        """:math:`n_{{sb}^e_{x}}` - number of soft state bounds at t=T"""
        return self.__nsbx_e

    @property
    def nsbu(self):
        """:math:`n_{{sb}_u}` - number of soft input bounds"""
        return self.__nsbu

    @property
    def nsg(self):
        """:math:`n_{{sg}}` - number of soft general linear constraints"""
        return self.__nsg

    @property
    def nsg_e(self):
        """:math:`n_{{sg}^e}` - number of soft general linear constraints at t=T"""
        return self.__nsg_e

    @property
    def nsh(self):
        """:math:`n_{{sh}}` - number of soft nonlinear constraints"""
        return self.__nsh

    @property
    def nsh_e(self):
        """:math:`n_{{sh}}^e` - number of soft nonlinear constraints at t=T"""
        return self.__nsh_e

    @property
    def nsphi(self):
        """:math:`n_{{s\phi}}` - number of soft convex-over-nonlinear constraints"""
        return self.__nsphi

    @property
    def nsphi_e(self):
        """:math:`n_{{s\phi}^e}` - number of soft convex-over-nonlinear constraints at t=T"""
        return self.__nsphi_e

    @property
    def ns(self):
        """:math:`n_{s}` - total number of slacks"""
        return self.__ns

    @property
    def ns_e(self):
        """:math:`n_{s}^e` - total number of slacks at t=T"""
        return self.__ns_e

    @property
    def ng(self):
        """:math:`n_{g}` - number of general polytopic constraints"""
        return self.__ng

    @property
    def ng_e(self):
        """:math:`n_{g}^e` - number of general polytopic constraints at t=T"""
        return self.__ng_e

    @property
    def N(self):
        """:math:`N` - prediction horizon"""
        return self.__N

    @nx.setter
    def nx(self, nx):
        if type(nx) == int and nx > 0:
            self.__nx = nx
        else:
            raise Exception('Invalid nx value, expected positive integer. Exiting.')

    @nz.setter
    def nz(self, nz):
        if type(nz) == int and nz > -1:
            self.__nz = nz
        else:
            raise Exception('Invalid nz value, expected nonnegative integer. Exiting.')

    @nu.setter
    def nu(self, nu):
        if type(nu) == int and nu > -1:
            self.__nu = nu
        else:
            raise Exception('Invalid nu value, expected nonnegative integer. Exiting.')

    @np.setter
    def np(self, np):
        if type(np) == int and np > -1:
            self.__np = np
        else:
            raise Exception('Invalid np value, expected nonnegative integer. Exiting.')

    @ny.setter
    def ny(self, ny):
        if type(ny) == int and ny > -1:
            self.__ny = ny
        else:
            raise Exception('Invalid ny value, expected nonnegative integer. Exiting.')

    @ny_e.setter
    def ny_e(self, ny_e):
        if type(ny_e) == int and ny_e > -1:
            self.__ny_e = ny_e
        else:
            raise Exception('Invalid ny_e value, expected nonnegative integer. Exiting.')

    @nr.setter
    def nr(self, nr):
        if type(nr) == int and nr > -1:
            self.__nr = nr
        else:
            raise Exception('Invalid nr value, expected nonnegative integer. Exiting.')

    @nr_e.setter
    def nr_e(self, nr_e):
        if type(nr_e) == int and nr_e > -1:
            self.__nr_e = nr_e
        else:
            raise Exception('Invalid nr_e value, expected nonnegative integer. Exiting.')

    @nh.setter
    def nh(self, nh):
        if type(nh) == int and nh > -1:
            self.__nh = nh
        else:
            raise Exception('Invalid nh value, expected nonnegative integer. Exiting.')

    @nh_e.setter
    def nh_e(self, nh_e):
        if type(nh_e) == int and nh_e > -1:
            self.__nh_e = nh_e
        else:
            raise Exception('Invalid nh_e value, expected nonnegative integer. Exiting.')

    @nphi.setter
    def nphi(self, nphi):
        if type(nphi) == int and nphi > -1:
            self.__nphi = nphi
        else:
            raise Exception('Invalid nphi value, expected nonnegative integer. Exiting.')

    @nphi_e.setter
    def nphi_e(self, nphi_e):
        if type(nphi_e) == int and nphi_e > -1:
            self.__nphi_e = nphi_e
        else:
            raise Exception('Invalid nphi_e value, expected nonnegative integer. Exiting.')

    @nbx.setter
    def nbx(self, nbx):
        if isinstance(nbx, int) and nbx > -1:
            self.__nbx = nbx
        else:
            raise Exception('Invalid nbx value, expected nonnegative integer. Exiting.')

    @nbxe_0.setter
    def nbxe_0(self, nbxe_0):
        if isinstance(nbxe_0, int) and nbxe_0 > -1:
            self.__nbxe_0 = nbxe_0
        else:
            raise Exception('Invalid nbxe_0 value, expected nonnegative integer. Exiting.')

    @nbx_0.setter
    def nbx_0(self, nbx_0):
        if type(nbx_0) == int and nbx_0 > -1:
            self.__nbx_0 = nbx_0
        else:
            raise Exception('Invalid nbx_0 value, expected nonnegative integer. Exiting.')

    @nbx_e.setter
    def nbx_e(self, nbx_e):
        if type(nbx_e) == int and nbx_e > -1:
            self.__nbx_e = nbx_e
        else:
            raise Exception('Invalid nbx_e value, expected nonnegative integer. Exiting.')

    @nbu.setter
    def nbu(self, nbu):
        if type(nbu) == int and nbu > -1:
            self.__nbu = nbu
        else:
            raise Exception('Invalid nbu value, expected nonnegative integer. Exiting.')

    @nsbx.setter
    def nsbx(self, nsbx):
        if type(nsbx) == int and nsbx > -1:
            self.__nsbx = nsbx
        else:
            raise Exception('Invalid nsbx value, expected nonnegative integer. Exiting.')

    @nsbx_e.setter
    def nsbx_e(self, nsbx_e):
        if type(nsbx_e) == int and nsbx_e > -1:
            self.__nsbx_e = nsbx_e
        else:
            raise Exception('Invalid nsbx_e value, expected nonnegative integer. Exiting.')

    @nsbu.setter
    def nsbu(self, nsbu):
        if type(nsbu) == int and nsbu > -1:
            self.__nsbu = nsbu
        else:
            raise Exception('Invalid nsbu value, expected nonnegative integer. Exiting.')

    @nsg.setter
    def nsg(self, nsg):
        if isinstance(nsg, int) and nsg > -1:
            self.__nsg = nsg
        else:
            raise Exception('Invalid nsg value, expected nonnegative integer. Exiting.')

    @nsg_e.setter
    def nsg_e(self, nsg_e):
        if isinstance(nsg_e, int) and nsg_e > -1:
            self.__nsg_e = nsg_e
        else:
            raise Exception('Invalid nsg_e value, expected nonnegative integer. Exiting.')

    @nsh.setter
    def nsh(self, nsh):
        if isinstance(nsh, int) and nsh > -1:
            self.__nsh = nsh
        else:
            raise Exception('Invalid nsh value, expected nonnegative integer. Exiting.')

    @nsh_e.setter
    def nsh_e(self, nsh_e):
        if isinstance(nsh_e, int) and nsh_e > -1:
            self.__nsh_e = nsh_e
        else:
            raise Exception('Invalid nsh_e value, expected nonnegative integer. Exiting.')

    @nsphi.setter
    def nsphi(self, nsphi):
        if isinstance(nsphi, int) and nsphi > -1:
            self.__nsphi = nsphi
        else:
            raise Exception('Invalid nsphi value, expected nonnegative integer. Exiting.')

    @nsphi_e.setter
    def nsphi_e(self, nsphi_e):
        if isinstance(nsphi_e, int) and nsphi_e > -1:
            self.__nsphi_e = nsphi_e
        else:
            raise Exception('Invalid nsphi_e value, expected nonnegative integer. Exiting.')

    @ns.setter
    def ns(self, ns):
        if isinstance(ns, int) and ns > -1:
            self.__ns = ns
        else:
            raise Exception('Invalid ns value, expected nonnegative integer. Exiting.')

    @ns_e.setter
    def ns_e(self, ns_e):
        if isinstance(ns_e, int) and ns_e > -1:
            self.__ns_e = ns_e
        else:
            raise Exception('Invalid ns_e value, expected nonnegative integer. Exiting.')

    @ng.setter
    def ng(self, ng):
        if isinstance(ng, int) and ng > -1:
            self.__ng = ng
        else:
            raise Exception('Invalid ng value, expected nonnegative integer. Exiting.')

    @ng_e.setter
    def ng_e(self, ng_e):
        if isinstance(ng_e, int) and ng_e > -1:
            self.__ng_e = ng_e
        else:
            raise Exception('Invalid ng_e value, expected nonnegative integer. Exiting.')

    @N.setter
    def N(self, N):
        if isinstance(N, int) and N > 0:
            self.__N = N
        else:
            raise Exception('Invalid N value, expected positive integer. Exiting.')

    def set(self, attr, value):
        setattr(self, attr, value)


class AcadosOcpCost:
    """
    class containing the numerical data of the cost:
    in case of LINEAR_LS:
    stage cost is
    :math:`l(x,u,z) = || V_x x + V_u u + V_z z - y_{\\text{ref}}||^2_W`,
    terminal cost is
    :math:`m(x) = || V^e_x x - y_{\\text{ref}^e}||^2_{W^e}`
    """
    def __init__(self):
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
        # Mayer term
        self.__cost_type_e = 'LINEAR_LS'  # cost type for Mayer term
        self.__W_e         = np.zeros((0,0))
        self.__Vx_e        = np.zeros((0,0))
        self.__yref_e      = np.array([])
        self.__Zl_e        = np.array([])
        self.__Zu_e        = np.array([])
        self.__zl_e        = np.array([])
        self.__zu_e        = np.array([])

    # Lagrange term
    @property
    def cost_type(self):
        """cost type"""
        return self.__cost_type

    @property
    def W(self):
        """:math:`W` - weight matrix"""
        return self.__W

    @property
    def Vx(self):
        """:math:`V_x` - x matrix coefficient"""
        return self.__Vx

    @property
    def Vu(self):
        """:math:`V_u` - u matrix coefficient"""
        return self.__Vu

    @property
    def Vz(self):
        """:math:`V_z` - z matrix coefficient"""
        return self.__Vz

    @property
    def yref(self):
        """:math:`y_{\text{ref}}` - reference"""
        return self.__yref

    @property
    def Zl(self):
        """:math:`Z_l` - diagonal of Hessian wrt lower slack"""
        return self.__Zl

    @property
    def Zu(self):
        """:math:`Z_u` - diagonal of Hessian wrt upper slack"""
        return self.__Zu

    @property
    def zl(self):
        """:math:`z_l` - gradient wrt lower slack"""
        return self.__zl

    @property
    def zu(self):
        """:math:`z_u` - gradient wrt upper slack"""
        return self.__zu

    @cost_type.setter
    def cost_type(self, cost_type):

        cost_types = ('LINEAR_LS', 'NONLINEAR_LS', 'EXTERNAL')

        if type(cost_type) == str and cost_type in cost_types:
            self.__cost_type = cost_type
        else:
            raise Exception('Invalid cost_type value. Exiting.')

    @W.setter
    def W(self, W):
        if isinstance(W, np.ndarray) and len(W.shape) == 2:
            self.__W = W
        else:
            raise Exception('Invalid cost W value. ' \
                + 'Should be 2 dimensional numpy array. Exiting.')


    @Vx.setter
    def Vx(self, Vx):
        if isinstance(Vx, np.ndarray) and len(Vx.shape) == 2:
            self.__Vx = Vx
        else:
            raise Exception('Invalid cost Vx value. ' \
                + 'Should be 2 dimensional numpy array. Exiting.')

    @Vu.setter
    def Vu(self, Vu):
        if isinstance(Vu, np.ndarray) and len(Vu.shape) == 2:
            self.__Vu = Vu
        else:
            raise Exception('Invalid cost Vu value. ' \
                + 'Should be 2 dimensional numpy array. Exiting.')

    @Vz.setter
    def Vz(self, Vz):
        if isinstance(Vz, np.ndarray) and len(Vz.shape) == 2:
            self.__Vz = Vz
        else:
            raise Exception('Invalid cost Vz value. ' \
                + 'Should be 2 dimensional numpy array. Exiting.')

    @yref.setter
    def yref(self, yref):
        if isinstance(yref, np.ndarray):
            self.__yref = yref
        else:
            raise Exception('Invalid yref value, expected numpy array. Exiting.')

    @Zl.setter
    def Zl(self, Zl):
        if isinstance(Zl, np.ndarray):
            self.__Zl = Zl
        else:
            raise Exception('Invalid Zl value, expected numpy array. Exiting.')

    @Zu.setter
    def Zu(self, Zu):
        if isinstance(Zu, np.ndarray):
            self.__Zu = Zu
        else:
            raise Exception('Invalid Zu value, expected numpy array. Exiting.')

    @zl.setter
    def zl(self, zl):
        if isinstance(zl, np.ndarray):
            self.__zl = zl
        else:
            raise Exception('Invalid zl value, expected numpy array. Exiting.')

    @zu.setter
    def zu(self, zu):
        if isinstance(zu, np.ndarray):
            self.__zu = zu
        else:
            raise Exception('Invalid zu value, expected numpy array. Exiting.')

    # Mayer term
    @property
    def cost_type_e(self):
        """cost type for Mayer term, either LINEAR_LS, NONLINEAR_LS, AUTO"""
        return self.__cost_type_e

    @property
    def W_e(self):
        """:math:`W` - weight matrix"""
        return self.__W_e

    @property
    def Vx_e(self):
        """:math:`W^e` - weight matrix for Mayer term"""
        return self.__Vx_e

    @property
    def yref_e(self):
        """:math:`V_x^e` - x matrix coefficient for Mayer term"""
        return self.__yref_e

    @property
    def Zl_e(self):
        """:math:`Z_l^e` - diagonal of Hessian wrt upper slack for Mayer term"""
        return self.__Zl_e

    @property
    def Zu_e(self):
        """:math:`Z_l^e` - diagonal of Hessian wrt upper slack for Mayer term"""
        return self.__Zu_e

    @property
    def zl_e(self):
        """:math:`z_l^e` - gradient wrt lower slack for Mayer term"""
        return self.__zl_e

    @property
    def zu_e(self):
        """:math:`z_u^e` - gradient wrt upper slack for Mayer term"""
        return self.__zu_e

    @cost_type_e.setter
    def cost_type_e(self, cost_type_e):
        cost_types = ('LINEAR_LS', 'NONLINEAR_LS', 'EXTERNAL')

        if type(cost_type_e) == str and cost_type_e in cost_types:
            self.__cost_type_e = cost_type_e
        else:
            raise Exception('Invalid cost_type_e value. Exiting.')

    @W_e.setter
    def W_e(self, W_e):
        if isinstance(W_e, np.ndarray) and len(W_e.shape) == 2:
            self.__W_e = W_e
        else:
            raise Exception('Invalid cost W_e value. ' \
                + 'Should be 2 dimensional numpy array. Exiting.')

    @Vx_e.setter
    def Vx_e(self, Vx_e):
        if isinstance(Vx_e, np.ndarray) and len(Vx_e.shape) == 2:
            self.__Vx_e = Vx_e
        else:
            raise Exception('Invalid cost Vx_e value. ' \
                + 'Should be 2 dimensional numpy array. Exiting.')

    @yref_e.setter
    def yref_e(self, yref_e):
        if isinstance(yref_e, np.ndarray):
            self.__yref_e = yref_e
        else:
            raise Exception('Invalid yref_e value, expected numpy array. Exiting.')

    @Zl_e.setter
    def Zl_e(self, Zl_e):
        if isinstance(Zl_e, np.ndarray):
            self.__Zl_e = Zl_e
        else:
            raise Exception('Invalid Zl_e value, expected numpy array. Exiting.')

    @Zu_e.setter
    def Zu_e(self, Zu_e):
        if isinstance(Zu_e, np.ndarray):
            self.__Zu_e = Zu_e
        else:
            raise Exception('Invalid Zu_e value, expected numpy array. Exiting.')

    @zl_e.setter
    def zl_e(self, zl_e):
        if isinstance(zl_e, np.ndarray):
            self.__zl_e = zl_e
        else:
            raise Exception('Invalid zl_e value, expected numpy array. Exiting.')

    @zu_e.setter
    def zu_e(self, zu_e):
        if isinstance(zu_e, np.ndarray):
            self.__zu_e = zu_e
        else:
            raise Exception('Invalid zu_e value, expected numpy array. Exiting.')

    def set(self, attr, value):
        setattr(self, attr, value)


def print_J_to_idx_note():
    print("NOTE: J* matrix is converted to zero based vector idx* vector, which is returned here.")


class AcadosOcpConstraints:
    """
    class containing the description of the constraints
    """
    def __init__(self):
        self.__constr_type   = 'BGH'
        self.__constr_type_e = 'BGH'
        # initial x
        self.__lbx_0   = np.array([])
        self.__ubx_0   = np.array([])
        self.__idxbx_0 = np.array([])
        self.__idxbxe_0 = np.array([])
        # state bounds
        self.__lbx     = np.array([])
        self.__ubx     = np.array([])
        self.__idxbx   = np.array([])
        # bounds on x at t=T
        self.__lbx_e   = np.array([])
        self.__ubx_e   = np.array([])
        self.__idxbx_e = np.array([])
        # bounds on u
        self.__lbu     = np.array([])
        self.__ubu     = np.array([])
        self.__idxbu   = np.array([])
        # polytopic constraints
        self.__lg      = np.array([])
        self.__ug      = np.array([])
        self.__D       = np.zeros((0,0))
        self.__C       = np.zeros((0,0))
        # polytopic constraints at t=T
        self.__C_e     = np.zeros((0,0))
        self.__lg_e    = np.array([])
        self.__ug_e    = np.array([])
        # nonlinear constraints
        self.__lh      = np.array([])
        self.__uh      = np.array([])
        # nonlinear constraints at t=T
        self.__uh_e    = np.array([])
        self.__lh_e    = np.array([])
        # convex-over-nonlinear constraints
        self.__lphi    = np.array([])
        self.__uphi    = np.array([])
        # nonlinear constraints at t=T
        self.__uphi_e = np.array([])
        self.__lphi_e = np.array([])
        # SLACK BOUNDS
        # soft bounds on x
        self.__lsbx   = np.array([])
        self.__usbx   = np.array([])
        self.__idxsbx = np.array([])
        # soft bounds on u
        self.__lsbu   = np.array([])
        self.__usbu   = np.array([])
        self.__idxsbu = np.array([])
        # soft bounds on x at t=T
        self.__lsbx_e  = np.array([])
        self.__usbx_e  = np.array([])
        self.__idxsbx_e= np.array([])
        # soft bounds on general linear constraints
        self.__lsg    = np.array([])
        self.__usg    = np.array([])
        self.__idxsg  = np.array([])
        # soft bounds on nonlinear constraints
        self.__lsh    = np.array([])
        self.__ush    = np.array([])
        self.__idxsh  = np.array([])
        # soft bounds on nonlinear constraints
        self.__lsphi  = np.array([])
        self.__usphi  = np.array([])
        self.__idxsphi  = np.array([])
        # soft bounds on general linear constraints at t=T
        self.__lsg_e    = np.array([])
        self.__usg_e    = np.array([])
        self.__idxsg_e  = np.array([])
        # soft bounds on nonlinear constraints at t=T
        self.__lsh_e    = np.array([])
        self.__ush_e    = np.array([])
        self.__idxsh_e  = np.array([])
        # soft bounds on nonlinear constraints at t=T
        self.__lsphi_e    = np.array([])
        self.__usphi_e    = np.array([])
        self.__idxsphi_e  = np.array([])


    # types
    @property
    def constr_type(self):
        """Constraints type"""
        return self.__constr_type

    @property
    def constr_type_e(self):
        """Constraints type t=T"""
        return self.__constr_type_e

    # initial bounds on x
    @property
    def lbx_0(self):
        """:math:`\\underline{x_0}` - lower bounds on x0"""
        return self.__lbx_0

    @property
    def ubx_0(self):
        """:math:`\\bar{x_0}` - upper bounds on x0"""
        return self.__ubx_0

    @property
    def idxbx_0(self):
        """indexes of bounds on x0"""
        return self.__idxbx_0

    @property
    def idxbxe_0(self):
        """indexes of bounds on x0 that are equalities (set automatically)"""
        return self.__idxbxe_0

    # bounds on x
    @property
    def lbx(self):
        """:math:`\\underline{x}` - lower bounds on x"""
        return self.__lbx

    @property
    def ubx(self):
        """:math:`\\bar{x}` - upper bounds on x"""
        return self.__ubx

    @property
    def idxbx(self):
        """indexes of bounds on x (defines :math:`J_{bx}`)"""
        return self.__idxbx

    @property
    def Jbx(self):
        """:math:`J_{bx}` - matrix coefficient for bounds on x"""
        print_J_to_idx_note()
        return self.__idxbx

    # bounds on x at t=T
    @property
    def lbx_e(self):
        """:math:`\\underline{x}^e` - lower bounds on x at t=T"""
        return self.__lbx_e

    @property
    def ubx_e(self):
        """:math:`\\bar{x}^e` - upper bounds on x at t=T"""
        return self.__ubx_e

    @property
    def idxbx_e(self):
        """indexes for bounds on x at t=T (defines :math:`J_{bx}^e`)"""
        return self.__idxbx_e

    @property
    def Jbx_e(self):
        """:math:`J_{bx}^e` matrix coefficient for bounds on x at t=T"""
        print_J_to_idx_note()
        return self.__idxbx_e

    # bounds on u
    @property
    def lbu(self):
        """:math:`\\underline{u}` - lower bounds on u"""
        return self.__lbu

    @property
    def ubu(self):
        """:math:`\\bar{u}` - upper bounds on u"""
        return self.__ubu

    @property
    def idxbu(self):
        """indexes of bounds on u (defines :math:`J_{bu}`)"""
        return self.__idxbu

    @property
    def Jbu(self):
        """:math:`J_{bu}` - matrix coefficient for bounds on u"""
        print_J_to_idx_note()
        return self.__idxbu

    # polytopic constraints
    @property
    def C(self):
        """:math:`C` - C matrix in lg <= D * u + C * x <= ug"""
        return self.__C

    @property
    def D(self):
        """:math:`D` - D matrix in lg <= D * u + C * x <= ug"""
        return self.__D

    @property
    def lg(self):
        """:math:`\\underline{g}` - lower bound for general polytopic inequalities"""
        return self.__lg

    @property
    def ug(self):
        """:math:`\\bar{g}` - upper bound for general polytopic inequalities"""
        return self.__ug

    # polytopic constraints at t=T
    @property
    def C_e(self):
        """:math:`C^e` - C matrix at t=T"""
        return self.__C_e

    @property
    def lg_e(self):
        """:math:`\\underline{g}^e` - lower bound on general polytopic inequalities at t=T"""
        return self.__lg_e

    @property
    def ug_e(self):
        """:math:`\\bar{g}^e` - upper bound on general polytopic inequalities at t=T"""
        return self.__ug_e


    # nonlinear constraints
    @property
    def lh(self):
        """:math:`\\underline{h}` - lower bound for nonlinear inequalities"""
        return self.__lh

    @property
    def uh(self):
        """:math:`\\bar{h}` - upper bound for nonlinear inequalities"""
        return self.__uh

    # nonlinear constraints at t=T
    @property
    def lh_e(self):
        """:math:`\\underline{h}^e` - lower bound on nonlinear inequalities at t=T"""
        return self.__lh_e

    @property
    def uh_e(self):
        """:math:`\\bar{h}^e` - upper bound on nonlinear inequalities at t=T"""
        return self.__uh_e

    # convex-over-nonlinear constraints
    @property
    def lphi(self):
        """:math:`\\underline{\phi}` - lower bound for convex-over-nonlinear inequalities"""
        return self.__lphi

    @property
    def uphi(self):
        """:math:`\\bar{\phi}` - upper bound for convex-over-nonlinear inequalities"""
        return self.__uphi

    # convex-over-nonlinear constraints at t=T
    @property
    def lphi_e(self):
        """:math:`\\underline{\phi}^e` - lower bound on convex-over-nonlinear inequalities at t=T"""
        return self.__lphi_e

    @property
    def uphi_e(self):
        """:math:`\\bar{\phi}^e` - upper bound on convex-over-nonlinear inequalities at t=T"""
        return self.__uphi_e


    # SLACK bounds
    # soft bounds on x
    @property
    def lsbx(self):
        """lower bounds on slacks corresponding to soft lower bounds on x"""
        return self.__lsbx

    @property
    def usbx(self):
        """upper bounds on slacks corresponding to soft upper bounds on x"""
        return self.__usbx

    @property
    def idxsbx(self):
        """indexes of soft bounds on x within the indices of bounds on x"""
        return self.__idxsbx

    @property
    def Jsbx(self):
        """:math:`J_{sbx}` - matrix coefficient for soft bounds on x"""
        print_J_to_idx_note()
        return self.__idxsbx

    # soft bounds on u
    @property
    def lsbu(self):
        """lower bounds on slacks corresponding to soft lower bounds on u"""
        return self.__lsbu

    @property
    def usbu(self):
        """upper bounds on slacks corresponding to soft upper bounds on u"""
        return self.__usbu

    @property
    def idxsbu(self):
        """indexes of soft bounds on u within the indices of bounds on u"""
        return self.__idxsbu

    @property
    def Jsbu(self):
        """:math:`J_{sbu}` - matrix coefficient for soft bounds on u"""
        print_J_to_idx_note()
        return self.__idxsbu

    # soft bounds on x at t=T
    @property
    def lsbx_e(self):
        """lower bounds on slacks corresponding to soft lower bounds on x at t=T"""
        return self.__lsbx_e

    @property
    def usbx_e(self):
        """upper bounds on slacks corresponding to soft upper bounds on x at t=T"""
        return self.__usbx_e

    @property
    def idxsbx_e(self):
        """indexes of soft bounds on x at t=T, within the indices of bounds on x at t=T"""
        return self.__idxsbx_e

    @property
    def Jsbx_e(self):
        """:math:`J_{sbx}^e` - matrix coefficient for soft bounds on x at t=T"""
        print_J_to_idx_note()
        return self.__idxsbx_e

    # soft general linear constraints
    @property
    def lsg(self):
        """lower bounds on slacks corresponding to soft lower bounds for general linear constraints"""
        return self.__lsg

    @property
    def usg(self):
        """upper bounds on slacks corresponding to soft upper bounds for general linear constraints"""
        return self.__usg

    @property
    def idxsg(self):
        """indexes of soft general linear constraints within the indices of general linear constraints"""
        return self.__idxsg

    @property
    def Jsg(self):
        """:math:`J_{sg}` - matrix coefficient for soft bounds on general linear constraints"""
        print_J_to_idx_note()
        return self.__idxsg

    # soft nonlinear constraints
    @property
    def lsh(self):
        """lower bounds on slacks corresponding to soft lower bounds for nonlinear constraints"""
        return self.__lsh

    @property
    def ush(self):
        """upper bounds on slacks corresponding to soft upper bounds for nonlinear constraints"""
        return self.__ush

    @property
    def idxsh(self):
        """indexes of soft nonlinear constraints within the indices of nonlinear constraints"""
        return self.__idxsh

    @property
    def Jsh(self):
        """:math:`J_{sh}` - matrix coefficient for soft bounds on nonlinear constraints"""
        print_J_to_idx_note()
        return self.__idxsh

    # soft bounds on convex-over-nonlinear constraints
    @property
    def lsphi(self):
        """lower bounds on slacks corresponding to soft lower bounds for convex-over-nonlinear constraints"""
        return self.__lsphi

    @property
    def usphi(self):
        """upper bounds on slacks corresponding to soft upper bounds for convex-over-nonlinear constraints"""
        return self.__usphi

    @property
    def idxsphi(self):
        """indexes of soft convex-over-nonlinear constraints within the indices of nonlinear constraints"""
        return self.__idxsphi

    @property
    def Jsphi(self):
        """:math:`J_{s, \phi}` - matrix coefficient for soft bounds on convex-over-nonlinear constraints"""
        print_J_to_idx_note()
        return self.__idxsphi


    # soft bounds on general linear constraints at t=T
    @property
    def lsg_e(self):
        """lower bounds on slacks corresponding to soft lower bounds for general linear constraints at t=T"""
        return self.__lsg_e

    @property
    def usg_e(self):
        """upper bounds on slacks corresponding to soft upper bounds for general linear constraints at t=T"""
        return self.__usg_e

    @property
    def idxsg_e(self):
        """indexes of soft general linear constraints at t=T within the indices of general linear constraints at t=T"""
        return self.__idxsg_e

    @property
    def Jsg_e(self):
        """:math:`J_{s,h}^e` - matrix coefficient for soft bounds on general linear constraints at t=T"""
        print_J_to_idx_note()
        return self.__idxsg_e


    # soft bounds on nonlinear constraints at t=T
    @property
    def lsh_e(self):
        """lower bounds on slacks corresponding to soft lower bounds for nonlinear constraints at t=T"""
        return self.__lsh_e

    @property
    def ush_e(self):
        """upper bounds on slacks corresponding to soft upper bounds for nonlinear constraints at t=T"""
        return self.__ush_e

    @property
    def idxsh_e(self):
        """indexes of soft nonlinear constraints at t=T within the indices of nonlinear constraints at t=T"""
        return self.__idxsh_e

    @property
    def Jsh_e(self):
        """:math:`J_{s,h}^e` - matrix coefficient for soft bounds on nonlinear constraints at t=T"""
        print_J_to_idx_note()
        return self.__idxsh_e

    # soft bounds on convex-over-nonlinear constraints at t=T
    @property
    def lsphi_e(self):
        """lower bounds on slacks corresponding to soft lower bounds for convex-over-nonlinear constraints at t=T"""
        return self.__lsphi_e

    @property
    def usphi_e(self):
        """upper bounds on slacks corresponding to soft upper bounds for convex-over-nonlinear constraints at t=T"""
        return self.__usphi_e

    @property
    def idxsphi_e(self):
        """indexes of soft nonlinear constraints at t=T within the indices of nonlinear constraints at t=T"""
        return self.__idxsphi_e

    @property
    def Jsphi_e(self):
        """:math:`J_{sh}^e` - matrix coefficient for soft bounds on convex-over-nonlinear constraints at t=T"""
        print_J_to_idx_note()
        return self.__idxsphi_e

    @property
    def x0(self):
        """:math:`\\bar{x}_0` - initial state"""
        print("x0 is converted to lbx_0, ubx_0, idxbx_0")
        print("idxbx_0: ", self.__idxbx_0)
        print("lbx_0: ", self.__lbx_0)
        print("ubx_0: ", self.__ubx_0)
        print("idxbxe_0: ", self.__idxbxe_0)
        return None

    # SETTERS
    @constr_type.setter
    def constr_type(self, constr_type):
        constr_types = ('BGH', 'BGP')

        if type(constr_type) == str and constr_type in constr_types:
            self.__constr_type = constr_type
        else:
            raise Exception('Invalid constr_type value. Possible values are:\n\n' \
                    + ',\n'.join(constr_types) + '.\n\nYou have: ' + constr_type + '.\n\nExiting.')

    @constr_type_e.setter
    def constr_type_e(self, constr_type_e):
        constr_types = ('BGH', 'BGP')

        if type(constr_type_e) == str and constr_type_e in constr_types:
            self.__constr_type_e = constr_type_e
        else:
            raise Exception('Invalid constr_type_e value. Possible values are:\n\n' \
                    + ',\n'.join(constr_types) + '.\n\nYou have: ' + constr_type_e + '.\n\nExiting.')

    # initial x
    @lbx_0.setter
    def lbx_0(self, lbx_0):
        if type(lbx_0) == np.ndarray:
            self.__lbx_0 = lbx_0
        else:
            raise Exception('Invalid lbx_0 value. Exiting.')

    @ubx_0.setter
    def ubx_0(self, ubx_0):
        if type(ubx_0) == np.ndarray:
            self.__ubx_0 = ubx_0
        else:
            raise Exception('Invalid ubx_0 value. Exiting.')

    @idxbx_0.setter
    def idxbx_0(self, idxbx_0):
        if isinstance(idxbx_0, np.ndarray):
            self.__idxbx_0 = idxbx_0
        else:
            raise Exception('Invalid idxbx_0 value. Exiting.')

    @idxbxe_0.setter
    def idxbxe_0(self, idxbxe_0):
        if isinstance(idxbxe_0, np.ndarray):
            self.__idxbxe_0 = idxbxe_0
        else:
            raise Exception('Invalid idxbxe_0 value. Exiting.')


    @x0.setter
    def x0(self, x0):
        if type(x0) == np.ndarray:
            self.__lbx_0 = x0
            self.__ubx_0 = x0
            self.__idxbx_0 = np.arange(x0.size)
            self.__idxbxe_0 = np.arange(x0.size)
        else:
            raise Exception('Invalid x0 value. Exiting.')

    # bounds on x
    @lbx.setter
    def lbx(self, lbx):
        if type(lbx) == np.ndarray:
            self.__lbx = lbx
        else:
            raise Exception('Invalid lbx value. Exiting.')

    @ubx.setter
    def ubx(self, ubx):
        if type(ubx) == np.ndarray:
            self.__ubx = ubx
        else:
            raise Exception('Invalid ubx value. Exiting.')

    @idxbx.setter
    def idxbx(self, idxbx):
        if type(idxbx) == np.ndarray:
            self.__idxbx = idxbx
        else:
            raise Exception('Invalid idxbx value. Exiting.')

    @Jbx.setter
    def Jbx(self, Jbx):
        if type(Jbx) == np.ndarray:
            self.__idxbx = J_to_idx(Jbx)
        else:
            raise Exception('Invalid Jbx value. Exiting.')

    # bounds on u
    @lbu.setter
    def lbu(self, lbu):
        if type(lbu) == np.ndarray:
            self.__lbu = lbu
        else:
            raise Exception('Invalid lbu value. Exiting.')

    @ubu.setter
    def ubu(self, ubu):
        if type(ubu) == np.ndarray:
            self.__ubu = ubu
        else:
            raise Exception('Invalid ubu value. Exiting.')

    @idxbu.setter
    def idxbu(self, idxbu):
        if type(idxbu) == np.ndarray:
            self.__idxbu = idxbu
        else:
            raise Exception('Invalid idxbu value. Exiting.')

    @Jbu.setter
    def Jbu(self, Jbu):
        if type(Jbu) == np.ndarray:
            self.__idxbu = J_to_idx(Jbu)
        else:
            raise Exception('Invalid Jbu value. Exiting.')

    # bounds on x at t=T
    @lbx_e.setter
    def lbx_e(self, lbx_e):
        if type(lbx_e) == np.ndarray:
            self.__lbx_e = lbx_e
        else:
            raise Exception('Invalid lbx_e value. Exiting.')

    @ubx_e.setter
    def ubx_e(self, ubx_e):
        if type(ubx_e) == np.ndarray:
            self.__ubx_e = ubx_e
        else:
            raise Exception('Invalid ubx_e value. Exiting.')

    @idxbx_e.setter
    def idxbx_e(self, idxbx_e):
        if type(idxbx_e) == np.ndarray:
            self.__idxbx_e = idxbx_e
        else:
            raise Exception('Invalid idxbx_e value. Exiting.')

    @Jbx_e.setter
    def Jbx_e(self, Jbx_e):
        if type(Jbx_e) == np.ndarray:
            self.__idxbx_e = J_to_idx(Jbx_e)
        else:
            raise Exception('Invalid Jbx_e value. Exiting.')

    # polytopic constraints
    @D.setter
    def D(self, D):
        if isinstance(D, np.ndarray) and len(D.shape) == 2:
            self.__D = D
        else:
            raise Exception('Invalid constraint D value.' \
                + 'Should be 2 dimensional numpy array. Exiting.')

    @C.setter
    def C(self, C):
        if isinstance(C, np.ndarray) and len(C.shape) == 2:
            self.__C = C
        else:
            raise Exception('Invalid constraint C value.' \
                + 'Should be 2 dimensional numpy array. Exiting.')

    @lg.setter
    def lg(self, lg):
        if type(lg) == np.ndarray:
            self.__lg = lg
        else:
            raise Exception('Invalid lg value. Exiting.')

    @ug.setter
    def ug(self, ug):
        if type(ug) == np.ndarray:
            self.__ug = ug
        else:
            raise Exception('Invalid ug value. Exiting.')

    # polytopic constraints at t=T
    @C_e.setter
    def C_e(self, C_e):
        if isinstance(C_e, np.ndarray) and len(C_e.shape) == 2:
            self.__C_e = C_e
        else:
            raise Exception('Invalid constraint C_e value.' \
                + 'Should be 2 dimensional numpy array. Exiting.')

    @lg_e.setter
    def lg_e(self, lg_e):
        if type(lg_e) == np.ndarray:
            self.__lg_e = lg_e
        else:
            raise Exception('Invalid lg_e value. Exiting.')

    @ug_e.setter
    def ug_e(self, ug_e):
        if type(ug_e) == np.ndarray:
            self.__ug_e = ug_e
        else:
            raise Exception('Invalid ug_e value. Exiting.')

    # nonlinear constraints
    @lh.setter
    def lh(self, lh):
        if type(lh) == np.ndarray:
            self.__lh = lh
        else:
            raise Exception('Invalid lh value. Exiting.')

    @uh.setter
    def uh(self, uh):
        if type(uh) == np.ndarray:
            self.__uh = uh
        else:
            raise Exception('Invalid uh value. Exiting.')

    # convex-over-nonlinear constraints
    @lphi.setter
    def lphi(self, lphi):
        if type(lphi) == np.ndarray:
            self.__lphi = lphi
        else:
            raise Exception('Invalid lphi value. Exiting.')

    @uphi.setter
    def uphi(self, uphi):
        if type(uphi) == np.ndarray:
            self.__uphi = uphi
        else:
            raise Exception('Invalid uphi value. Exiting.')

    # nonlinear constraints at t=T
    @lh_e.setter
    def lh_e(self, lh_e):
        if type(lh_e) == np.ndarray:
            self.__lh_e = lh_e
        else:
            raise Exception('Invalid lh_e value. Exiting.')

    @uh_e.setter
    def uh_e(self, uh_e):
        if type(uh_e) == np.ndarray:
            self.__uh_e = uh_e
        else:
            raise Exception('Invalid uh_e value. Exiting.')

    # convex-over-nonlinear constraints at t=T
    @lphi_e.setter
    def lphi_e(self, lphi_e):
        if type(lphi_e) == np.ndarray:
            self.__lphi_e = lphi_e
        else:
            raise Exception('Invalid lphi_e value. Exiting.')

    @uphi_e.setter
    def uphi_e(self, uphi_e):
        if type(uphi_e) == np.ndarray:
            self.__uphi_e = uphi_e
        else:
            raise Exception('Invalid uphi_e value. Exiting.')

    # SLACK bounds
    # soft bounds on x
    @lsbx.setter
    def lsbx(self, lsbx):
        if type(lsbx) == np.ndarray:
            self.__lsbx = lsbx
        else:
            raise Exception('Invalid lsbx value. Exiting.')

    @usbx.setter
    def usbx(self, usbx):
        if type(usbx) == np.ndarray:
            self.__usbx = usbx
        else:
            raise Exception('Invalid usbx value. Exiting.')

    @idxsbx.setter
    def idxsbx(self, idxsbx):
        if type(idxsbx) == np.ndarray:
            self.__idxsbx = idxsbx
        else:
            raise Exception('Invalid idxsbx value. Exiting.')

    @Jsbx.setter
    def Jsbx(self, Jsbx):
        if isinstance(Jsbx, np.ndarray):
            self.__idxsbx = J_to_idx_slack(Jsbx)
        else:
            raise Exception('Invalid Jsbx value, expected numpy array. Exiting.')

    # soft bounds on u
    @lsbu.setter
    def lsbu(self, lsbu):
        if type(lsbu) == np.ndarray:
            self.__lsbu = lsbu
        else:
            raise Exception('Invalid lsbu value. Exiting.')

    @usbu.setter
    def usbu(self, usbu):
        if type(usbu) == np.ndarray:
            self.__usbu = usbu
        else:
            raise Exception('Invalid usbu value. Exiting.')

    @idxsbu.setter
    def idxsbu(self, idxsbu):
        if type(idxsbu) == np.ndarray:
            self.__idxsbu = idxsbu
        else:
            raise Exception('Invalid idxsbu value. Exiting.')

    @Jsbu.setter
    def Jsbu(self, Jsbu):
        if type(Jsbu) == np.ndarray:
            self.__idxsbu = J_to_idx_slack(Jsbu)
        else:
            raise Exception('Invalid Jsbu value. Exiting.')

    # soft bounds on x at t=T
    @lsbx_e.setter
    def lsbx_e(self, lsbx_e):
        if type(lsbx_e) == np.ndarray:
            self.__lsbx_e = lsbx_e
        else:
            raise Exception('Invalid lsbx_e value. Exiting.')

    @usbx_e.setter
    def usbx_e(self, usbx_e):
        if type(usbx_e) == np.ndarray:
            self.__usbx_e = usbx_e
        else:
            raise Exception('Invalid usbx_e value. Exiting.')

    @idxsbx_e.setter
    def idxsbx_e(self, idxsbx_e):
        if type(idxsbx_e) == np.ndarray:
            self.__idxsbx_e = idxsbx_e
        else:
            raise Exception('Invalid idxsbx_e value. Exiting.')

    @Jsbx_e.setter
    def Jsbx_e(self, Jsbx_e):
        if type(Jsbx_e) == np.ndarray:
            self.__idxsbx_e = J_to_idx_slack(Jsbx_e)
        else:
            raise Exception('Invalid Jsbx_e value. Exiting.')


    # soft bounds on general linear constraints
    @lsg.setter
    def lsg(self, lsg):
        if isinstance(lsg, np.ndarray):
            self.__lsg = lsg
        else:
            raise Exception('Invalid lsg value. Exiting.')

    @usg.setter
    def usg(self, usg):
        if isinstance(usg, np.ndarray):
            self.__usg = usg
        else:
            raise Exception('Invalid usg value. Exiting.')

    @idxsg.setter
    def idxsg(self, idxsg):
        if isinstance(idxsg, np.ndarray):
            self.__idxsg = idxsg
        else:
            raise Exception('Invalid idxsg value. Exiting.')

    @Jsg.setter
    def Jsg(self, Jsg):
        if isinstance(Jsg, np.ndarray):
            self.__idxsg = J_to_idx_slack(Jsg)
        else:
            raise Exception('Invalid Jsg value, expected numpy array. Exiting.')


    # soft bounds on nonlinear constraints
    @lsh.setter
    def lsh(self, lsh):
        if type(lsh) == np.ndarray:
            self.__lsh = lsh
        else:
            raise Exception('Invalid lsh value. Exiting.')

    @ush.setter
    def ush(self, ush):
        if type(ush) == np.ndarray:
            self.__ush = ush
        else:
            raise Exception('Invalid ush value. Exiting.')

    @idxsh.setter
    def idxsh(self, idxsh):
        if type(idxsh) == np.ndarray:
            self.__idxsh = idxsh
        else:
            raise Exception('Invalid idxsh value. Exiting.')

    # soft bounds on convex-over-nonlinear constraints
    @lsphi.setter
    def lsphi(self, lsphi):
        if type(lsphi) == np.ndarray:
            self.__lsphi = lsphi
        else:
            raise Exception('Invalid lsphi value. Exiting.')

    @usphi.setter
    def usphi(self, usphi):
        if type(usphi) == np.ndarray:
            self.__usphi = usphi
        else:
            raise Exception('Invalid usphi value. Exiting.')

    @idxsphi.setter
    def idxsphi(self, idxsphi):
        if type(idxsphi) == np.ndarray:
            self.__idxsphi = idxsphi
        else:
            raise Exception('Invalid idxsphi value. Exiting.')

    @Jsphi.setter
    def Jsphi(self, Jsphi):
        if isinstance(Jsphi, np.ndarray):
            self.__idxsphi = J_to_idx_slack(Jsphi)
        else:
            raise Exception('Invalid Jsphi value, expected numpy array. Exiting.')

    # soft bounds on general linear constraints at t=T
    @lsg_e.setter
    def lsg_e(self, lsg_e):
        if isinstance(lsg_e, np.ndarray):
            self.__lsg_e = lsg_e
        else:
            raise Exception('Invalid lsg_e value. Exiting.')

    @usg_e.setter
    def usg_e(self, usg_e):
        if isinstance(usg_e, np.ndarray):
            self.__usg_e = usg_e
        else:
            raise Exception('Invalid usg_e value. Exiting.')

    @idxsg_e.setter
    def idxsg_e(self, idxsg_e):
        if isinstance(idxsg_e, np.ndarray):
            self.__idxsg_e = idxsg_e
        else:
            raise Exception('Invalid idxsg_e value. Exiting.')

    @Jsg_e.setter
    def Jsg_e(self, Jsg_e):
        if isinstance(Jsg_e, np.ndarray):
            self.__idxsg_e = J_to_idx_slack(Jsg_e)
        else:
            raise Exception('Invalid Jsg_e value, expected numpy array. Exiting.')

    # soft bounds on nonlinear constraints at t=T
    @lsh_e.setter
    def lsh_e(self, lsh_e):
        if isinstance(lsh_e, np.ndarray):
            self.__lsh_e = lsh_e
        else:
            raise Exception('Invalid lsh_e value. Exiting.')

    @ush_e.setter
    def ush_e(self, ush_e):
        if isinstance(ush_e, np.ndarray):
            self.__ush_e = ush_e
        else:
            raise Exception('Invalid ush_e value. Exiting.')

    @idxsh_e.setter
    def idxsh_e(self, idxsh_e):
        if isinstance(idxsh_e, np.ndarray):
            self.__idxsh_e = idxsh_e
        else:
            raise Exception('Invalid idxsh_e value. Exiting.')

    # soft bounds on convex-over-nonlinear constraints at t=T
    @lsphi_e.setter
    def lsphi_e(self, lsphi_e):
        if isinstance(lsphi_e, np.ndarray):
            self.__lsphi_e = lsphi_e
        else:
            raise Exception('Invalid lsphi_e value. Exiting.')

    @usphi_e.setter
    def usphi_e(self, usphi_e):
        if isinstance(usphi_e, np.ndarray):
            self.__usphi_e = usphi_e
        else:
            raise Exception('Invalid usphi_e value. Exiting.')

    @idxsphi_e.setter
    def idxsphi_e(self, idxsphi_e):
        if isinstance(idxsphi_e, np.ndarray):
            self.__idxsphi_e = idxsphi_e
        else:
            raise Exception('Invalid idxsphi_e value. Exiting.')

    @Jsphi_e.setter
    def Jsphi_e(self, Jsphi_e):
        if isinstance(Jsphi_e, np.ndarray):
            self.__idxsphi_e = J_to_idx_slack(Jsphi_e)
        else:
            raise Exception('Invalid Jsphi_e value. Exiting.')

    def set(self, attr, value):
        setattr(self, attr, value)


class AcadosOcpOptions:
    """
    class containing the description of the solver options
    """
    def __init__(self):
        self.__qp_solver        = 'PARTIAL_CONDENSING_HPIPM'  # qp solver to be used in the NLP solver
        self.__hessian_approx   = 'GAUSS_NEWTON'              # hessian approximation
        self.__integrator_type  = 'ERK'                       # integrator type
        self.__tf               = None                        # prediction horizon
        self.__nlp_solver_type  = 'SQP_RTI'                   # NLP solver
        self.__nlp_solver_step_length = 1.0                   # fixed Newton step length
        self.__levenberg_marquardt = 0.0
        self.__sim_method_num_stages  = 4                     # number of stages in the integrator
        self.__sim_method_num_steps   = 1                     # number of steps in the integrator
        self.__sim_method_newton_iter = 3                     # number of Newton iterations in simulation method
        self.__sim_method_jac_reuse = False
        self.__qp_solver_tol_stat = None                      # QP solver stationarity tolerance
        self.__qp_solver_tol_eq   = None                      # QP solver equality tolerance
        self.__qp_solver_tol_ineq = None                      # QP solver inequality
        self.__qp_solver_tol_comp = None                      # QP solver complementarity
        self.__qp_solver_iter_max = 50                        # QP solver max iter
        self.__qp_solver_cond_N = None                        # QP solver: new horizon after partial condensing
        self.__nlp_solver_tol_stat = 1e-6                     # NLP solver stationarity tolerance
        self.__nlp_solver_tol_eq   = 1e-6                     # NLP solver equality tolerance
        self.__nlp_solver_tol_ineq = 1e-6                     # NLP solver inequality
        self.__nlp_solver_tol_comp = 1e-6                     # NLP solver complementarity
        self.__nlp_solver_max_iter = 100                      # NLP solver maximum number of iterations
        self.__Tsim = None                                    # automatically calculated as tf/N
        self.__print_level = 0                                # print level
        self.__initialize_t_slacks = 0                        # possible values: 0, 1
        self.__model_external_shared_lib_dir   = None         # path to the the .so lib
        self.__model_external_shared_lib_name  = None         # name of the the .so lib
        self.__regularize_method = None
        self.__time_steps = None
        self.__shooting_nodes = None
        self.__exact_hess_cost = 1
        self.__exact_hess_dyn = 1
        self.__exact_hess_constr = 1
        self.__ext_cost_num_hess = 0


    @property
    def qp_solver(self):
        """QP solver to be used in the NLP solver"""
        return self.__qp_solver

    @property
    def hessian_approx(self):
        """Hessian approximation"""
        return self.__hessian_approx

    @property
    def integrator_type(self):
        """Integrator type"""
        return self.__integrator_type

    @property
    def nlp_solver_type(self):
        """NLP solver"""
        return self.__nlp_solver_type

    @property
    def regularize_method(self):
        """Regularization method for the Hessian"""
        return self.__regularize_method

    @property
    def nlp_solver_step_length(self):
        """Fixed Newton step length"""
        return self.__nlp_solver_step_length

    @property
    def levenberg_marquardt(self):
        """Factor for LM regularization"""
        return self.__levenberg_marquardt

    @property
    def sim_method_num_stages(self):
        """Number of stages in the integrator"""
        return self.__sim_method_num_stages

    @property
    def sim_method_num_steps(self):
        """Number of steps in the integrator"""
        return self.__sim_method_num_steps

    @property
    def sim_method_newton_iter(self):
        """Number of Newton iterations in simulation method"""
        return self.__sim_method_newton_iter

    @property
    def sim_method_jac_reuse(self):
        """Boolean determining if jacobians are reused within integrator"""
        return self.__sim_method_jac_reuse

    @property
    def qp_solver_tol_stat(self):
        """QP solver stationarity tolerance"""
        return self.__qp_solver_tol_stat

    @property
    def qp_solver_tol_eq(self):
        """QP solver equality tolerance"""
        return self.__qp_solver_tol_eq

    @property
    def qp_solver_tol_ineq(self):
        """QP solver inequality"""
        return self.__qp_solver_tol_ineq

    @property
    def qp_solver_tol_comp(self):
        """QP solver complementarity"""
        return self.__qp_solver_tol_comp

    @property
    def qp_solver_cond_N(self):
        """QP solver: New horizon after partial condensing"""
        return self.__qp_solver_cond_N

    @property
    def qp_solver_iter_max(self):
        """QP solver: maximum number of iterations"""
        return self.__qp_solver_iter_max

    @property
    def tol(self):
        """NLP solver tolerance"""
        return max([self.__nlp_solver_tol_eq, self.__nlp_solver_tol_ineq,\
                    self.__nlp_solver_tol_comp, self.__nlp_solver_tol_stat])

    @property
    def nlp_solver_tol_stat(self):
        """NLP solver stationarity tolerance"""
        return self.__nlp_solver_tol_stat

    @property
    def nlp_solver_tol_eq(self):
        """NLP solver equality tolerance"""
        return self.__nlp_solver_tol_eq

    @property
    def nlp_solver_tol_ineq(self):
        """NLP solver inequality tolerance"""
        return self.__nlp_solver_tol_ineq

    @property
    def nlp_solver_tol_comp(self):
        """NLP solver complementarity tolerance"""
        return self.__nlp_solver_tol_comp

    @property
    def nlp_solver_max_iter(self):
        """NLP solver maximum number of iterations"""
        return self.__nlp_solver_max_iter

    @property
    def time_steps(self):
        """Vector with time steps between the shooting nodes. Set automatically to uniform discretization if N, tf are provided"""
        return self.__time_steps

    @property
    def shooting_nodes(self):
        """Vector with the shooting nodes, time_steps will be computed from it automatically"""
        return self.__shooting_nodes

    @property
    def tf(self):
        """Prediction horizon"""
        return self.__tf

    @property
    def Tsim(self):
        """Time horizon for one integrator step"""
        return self.__Tsim

    @property
    def print_level(self):
        """Verbosity of printing"""
        return self.__print_level

    @property
    def model_external_shared_lib_dir(self):
        """Path to the .so lib"""
        return self.__model_external_shared_lib_dir

    @property
    def model_external_shared_lib_name(self):
        """Name of the .so lib"""
        return self.__model_external_shared_lib_name

    @property
    def exact_hess_constr(self):
        """Used in case of hessian_approx == 'EXACT'.\n
           Can be used to turn off exact hessian contributions from the constraints module"""
        return self.__exact_hess_constr

    @property
    def exact_hess_cost(self):
        """Used in case of hessian_approx == 'EXACT'.\n
           Can be used to turn off exact hessian contributions from the cost module"""
        return self.__exact_hess_cost

    @property
    def exact_hess_dyn(self):
        """Used in case of hessian_approx == 'EXACT'.\n
           Can be used to turn off exact hessian contributions from the dynamics module"""
        return self.__exact_hess_dyn

    @property
    def ext_cost_num_hess(self):
        """Determines if custom hessian approximation for cost contribution is used (> 0).\n
           Or if hessian contribution is evaluated exactly using CasADi external function (=0 - default)."""
        return self.__ext_cost_num_hess

    @qp_solver.setter
    def qp_solver(self, qp_solver):
        qp_solvers = ('PARTIAL_CONDENSING_HPIPM', 'PARTIAL_CONDENSING_QPOASES', \
                'FULL_CONDENSING_QPOASES', 'FULL_CONDENSING_HPIPM')

        if isinstance(qp_solver, str) and qp_solver in qp_solvers:
            self.__qp_solver = qp_solver
        else:
            raise Exception('Invalid qp_solver value. Possible values are:\n\n' \
                    + ',\n'.join(qp_solvers) + '.\n\nYou have: ' + qp_solver + '.\n\nExiting.')

    @regularize_method.setter
    def regularize_method(self, regularize_method):
        regularize_methods = ('NO_REGULARIZE', 'MIRROR', 'PROJECT', \
                                'PROJECT_REDUC_HESS', 'CONVEXIFY')

        if isinstance(regularize_method, str) and regularize_method in regularize_methods:
            self.__regularize_method = regularize_method
        else:
            raise Exception('Invalid regularize_method value. Possible values are:\n\n' \
                    + ',\n'.join(regularize_methods) + '.\n\nYou have: ' + regularize_method + '.\n\nExiting.')

    @hessian_approx.setter
    def hessian_approx(self, hessian_approx):
        hessian_approxs = ('GAUSS_NEWTON', 'EXACT')

        if type(hessian_approx) == str and hessian_approx in hessian_approxs:
            self.__hessian_approx = hessian_approx
        else:
            raise Exception('Invalid hessian_approx value. Possible values are:\n\n' \
                    + ',\n'.join(hessian_approxs) + '.\n\nYou have: ' + hessian_approx + '.\n\nExiting.')

    @integrator_type.setter
    def integrator_type(self, integrator_type):
        integrator_types = ('ERK', 'IRK', 'GNSF')

        if type(integrator_type) == str and integrator_type in integrator_types:
            self.__integrator_type = integrator_type
        else:
            raise Exception('Invalid integrator_type value. Possible values are:\n\n' \
                    + ',\n'.join(integrator_types) + '.\n\nYou have: ' + integrator_type + '.\n\nExiting.')

    @tf.setter
    def tf(self, tf):
        self.__tf = tf

    @time_steps.setter
    def time_steps(self, time_steps):
        self.__time_steps = time_steps

    @shooting_nodes.setter
    def shooting_nodes(self, shooting_nodes):
        self.__shooting_nodes = shooting_nodes


    @Tsim.setter
    def Tsim(self, Tsim):
        self.__Tsim = Tsim

    @sim_method_num_stages.setter
    def sim_method_num_stages(self, sim_method_num_stages):

        if type(sim_method_num_stages) == int:
            self.__sim_method_num_stages = sim_method_num_stages
        else:
            raise Exception('Invalid sim_method_num_stages value. sim_method_num_stages must be an integer. Exiting.')

    @sim_method_num_steps.setter
    def sim_method_num_steps(self, sim_method_num_steps):

        if type(sim_method_num_steps) == int:
            self.__sim_method_num_steps = sim_method_num_steps
        else:
            raise Exception('Invalid sim_method_num_steps value. sim_method_num_steps must be an integer. Exiting.')

    @sim_method_newton_iter.setter
    def sim_method_newton_iter(self, sim_method_newton_iter):

        if type(sim_method_newton_iter) == int:
            self.__sim_method_newton_iter = sim_method_newton_iter
        else:
            raise Exception('Invalid sim_method_newton_iter value. sim_method_newton_iter must be an integer. Exiting.')

    @sim_method_jac_reuse.setter
    def sim_method_jac_reuse(self, sim_method_jac_reuse):
        if sim_method_jac_reuse in (True, False):
            self.__sim_method_jac_reuse = sim_method_jac_reuse
        else:
            raise Exception('Invalid sim_method_jac_reuse value. sim_method_jac_reuse must be a Boolean.')

    @nlp_solver_type.setter
    def nlp_solver_type(self, nlp_solver_type):
        nlp_solver_types = ('SQP', 'SQP_RTI')

        if type(nlp_solver_type) == str and nlp_solver_type in nlp_solver_types:
            self.__nlp_solver_type = nlp_solver_type
        else:
            raise Exception('Invalid nlp_solver_type value. Possible values are:\n\n' \
                    + ',\n'.join(nlp_solver_types) + '.\n\nYou have: ' + nlp_solver_type + '.\n\nExiting.')

    @nlp_solver_step_length.setter
    def nlp_solver_step_length(self, nlp_solver_step_length):

        if type(nlp_solver_step_length) == float and nlp_solver_step_length > 0:
            self.__nlp_solver_step_length = nlp_solver_step_length
        else:
            raise Exception('Invalid nlp_solver_step_length value. nlp_solver_step_length must be a positive float. Exiting')

    @levenberg_marquardt.setter
    def levenberg_marquardt(self, levenberg_marquardt):
        if isinstance(levenberg_marquardt, float) and levenberg_marquardt >= 0:
            self.__levenberg_marquardt = levenberg_marquardt
        else:
            raise Exception('Invalid levenberg_marquardt value. levenberg_marquardt must be a positive float. Exiting')

    @qp_solver_tol_stat.setter
    def qp_solver_tol_stat(self, qp_solver_tol_stat):

        if type(qp_solver_tol_stat) == float and qp_solver_tol_stat > 0:
            self.__qp_solver_tol_stat = qp_solver_tol_stat
        else:
            raise Exception('Invalid qp_solver_tol_stat value. qp_solver_tol_stat must be a positive float. Exiting')

    @qp_solver_iter_max.setter
    def qp_solver_iter_max(self, qp_solver_iter_max):

        if isinstance(qp_solver_iter_max, int) and qp_solver_iter_max > 0:
            self.__qp_solver_iter_max = qp_solver_iter_max
        else:
            raise Exception('Invalid qp_solver_iter_max value. qp_solver_iter_max must be a positive int. Exiting')

    @qp_solver_cond_N.setter
    def qp_solver_cond_N(self, qp_solver_cond_N):

        if isinstance(qp_solver_cond_N, int) and qp_solver_cond_N > 0:
            self.__qp_solver_cond_N = qp_solver_cond_N
        else:
            raise Exception('Invalid qp_solver_cond_N value. qp_solver_cond_N must be a positive int. Exiting')

    @qp_solver_tol_eq.setter
    def qp_solver_tol_eq(self, qp_solver_tol_eq):

        if type(qp_solver_tol_eq) == float and qp_solver_tol_eq > 0:
            self.__qp_solver_tol_eq = qp_solver_tol_eq
        else:
            raise Exception('Invalid qp_solver_tol_eq value. qp_solver_tol_eq must be a positive float. Exiting')

    @qp_solver_tol_ineq.setter
    def qp_solver_tol_ineq(self, qp_solver_tol_ineq):

        if type(qp_solver_tol_ineq) == float and qp_solver_tol_ineq > 0:
            self.__qp_solver_tol_ineq = qp_solver_tol_ineq
        else:
            raise Exception('Invalid qp_solver_tol_ineq value. qp_solver_tol_ineq must be a positive float. Exiting')

    @qp_solver_tol_comp.setter
    def qp_solver_tol_comp(self, qp_solver_tol_comp):

        if type(qp_solver_tol_comp) == float and qp_solver_tol_comp > 0:
            self.__qp_solver_tol_comp = qp_solver_tol_comp
        else:
            raise Exception('Invalid qp_solver_tol_comp value. qp_solver_tol_comp must be a positive float. Exiting')

    @nlp_solver_tol_stat.setter
    def nlp_solver_tol_stat(self, nlp_solver_tol_stat):

        if type(nlp_solver_tol_stat) == float and nlp_solver_tol_stat > 0:
            self.__nlp_solver_tol_stat = nlp_solver_tol_stat
        else:
            raise Exception('Invalid nlp_solver_tol_stat value. nlp_solver_tol_stat must be a positive float. Exiting')

    @tol.setter
    def tol(self, tol):

        if type(tol) == float and tol > 0:
            self.__nlp_solver_tol_eq = tol
            self.__nlp_solver_tol_ineq = tol
            self.__nlp_solver_tol_stat = tol
            self.__nlp_solver_tol_comp = tol
        else:
            raise Exception('Invalid tol value. tol must be a positive float. Exiting')


    @nlp_solver_tol_eq.setter
    def nlp_solver_tol_eq(self, nlp_solver_tol_eq):

        if type(nlp_solver_tol_eq) == float and nlp_solver_tol_eq > 0:
            self.__nlp_solver_tol_eq = nlp_solver_tol_eq
        else:
            raise Exception('Invalid nlp_solver_tol_eq value. nlp_solver_tol_eq must be a positive float. Exiting')

    @nlp_solver_tol_ineq.setter
    def nlp_solver_tol_ineq(self, nlp_solver_tol_ineq):

        if type(nlp_solver_tol_ineq) == float and nlp_solver_tol_ineq > 0:
            self.__nlp_solver_tol_ineq = nlp_solver_tol_ineq
        else:
            raise Exception('Invalid nlp_solver_tol_ineq value. nlp_solver_tol_ineq must be a positive float. Exiting')

    @nlp_solver_tol_comp.setter
    def nlp_solver_tol_comp(self, nlp_solver_tol_comp):

        if type(nlp_solver_tol_comp) == float and nlp_solver_tol_comp > 0:
            self.__nlp_solver_tol_comp = nlp_solver_tol_comp
        else:
            raise Exception('Invalid nlp_solver_tol_comp value. nlp_solver_tol_comp must be a positive float. Exiting')

    @nlp_solver_max_iter.setter
    def nlp_solver_max_iter(self, nlp_solver_max_iter):

        if type(nlp_solver_max_iter) == int and nlp_solver_max_iter > 0:
            self.__nlp_solver_max_iter = nlp_solver_max_iter
        else:
            raise Exception('Invalid nlp_solver_max_iter value. nlp_solver_max_iter must be a positive int. Exiting')

    @print_level.setter
    def print_level(self, print_level):
        if type(print_level) == int and print_level >= 0:
            self.__print_level = print_level
        else:
            raise Exception('Invalid print_level value. print_level takes one of the values >=0. Exiting')

    @model_external_shared_lib_dir.setter
    def model_external_shared_lib_dir(self, model_external_shared_lib_dir):
        if type(model_external_shared_lib_dir) == str :
            self.__model_external_shared_lib_dir = model_external_shared_lib_dir
        else:
            raise Exception('Invalid model_external_shared_lib_dir value. Str expected.' \
            + '.\n\nYou have: ' + type(model_external_shared_lib_dir) + '.\n\nExiting.')

    @model_external_shared_lib_name.setter
    def model_external_shared_lib_name(self, model_external_shared_lib_name):
        if type(model_external_shared_lib_name) == str :
            if model_external_shared_lib_name[-3:] == '.so' : 
                raise Exception('Invalid model_external_shared_lib_name value. Remove the .so extension.' \
            + '.\n\nYou have: ' + type(model_external_shared_lib_name) + '.\n\nExiting.')
            else :
                self.__model_external_shared_lib_name = model_external_shared_lib_name
        else:
            raise Exception('Invalid model_external_shared_lib_name value. Str expected.' \
            + '.\n\nYou have: ' + type(model_external_shared_lib_name) + '.\n\nExiting.')

    @exact_hess_constr.setter
    def exact_hess_constr(self, exact_hess_constr):
        if exact_hess_constr in [0, 1]:
            self.__exact_hess_constr = exact_hess_constr
        else:
            raise Exception('Invalid exact_hess_constr value. exact_hess_constr takes one of the values 0, 1. Exiting')

    @exact_hess_cost.setter
    def exact_hess_cost(self, exact_hess_cost):
        if exact_hess_cost in [0, 1]:
            self.__exact_hess_cost = exact_hess_cost
        else:
            raise Exception('Invalid exact_hess_cost value. exact_hess_cost takes one of the values 0, 1. Exiting')

    @exact_hess_dyn.setter
    def exact_hess_dyn(self, exact_hess_dyn):
        if exact_hess_dyn in [0, 1]:
            self.__exact_hess_dyn = exact_hess_dyn
        else:
            raise Exception('Invalid exact_hess_dyn value. exact_hess_dyn takes one of the values 0, 1. Exiting')

    @ext_cost_num_hess.setter
    def ext_cost_num_hess(self, ext_cost_num_hess):
        if ext_cost_num_hess in [0, 1]:
            self.__ext_cost_num_hess = ext_cost_num_hess
        else:
            raise Exception('Invalid ext_cost_num_hess value. ext_cost_num_hess takes one of the values 0, 1. Exiting')

    def set(self, attr, value):
        setattr(self, attr, value)


class AcadosOcp:
    """
    class containing the full description of the optimal control problem
    """
    def __init__(self, acados_path=''):
        """
        Keyword arguments:
        acados_path -- path of your acados installation
        """
        if acados_path == '':
            acados_path = get_acados_path()

        self.dims = AcadosOcpDims()
        self.model = AcadosModel()
        self.cost = AcadosOcpCost()
        self.constraints = AcadosOcpConstraints()
        self.solver_options = AcadosOcpOptions()
		
        self.acados_include_path = f'{acados_path}/include'
        self.acados_lib_path = f'{acados_path}/lib'

        self.__parameter_values = np.array([])
        self.__problem_class = 'OCP'

    @property
    def parameter_values(self):
        """:math:`p` - initial values for parameter - can be updated stagewise"""
        return self.__parameter_values

    @parameter_values.setter
    def parameter_values(self, parameter_values):
        if isinstance(parameter_values, np.ndarray):
            self.__parameter_values = parameter_values
        else:
            raise Exception('Invalid parameter_values value. ' +
                            f'Expected numpy array, got {type(parameter_values)}.')

    def set(self, attr, value):
        # tokenize string
        tokens = attr.split('_', 1)
        if len(tokens) > 1:
            setter_to_call = getattr(getattr(self, tokens[0]), 'set')
        else:
            setter_to_call = getattr(self, 'set')

        setter_to_call(tokens[1], value)

        return
