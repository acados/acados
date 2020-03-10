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
        self.__ny      = None
        self.__ny_e    = None
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
            raise Exception('Invalid nx value. Exiting.')

    @nz.setter
    def nz(self, nz):
        if type(nz) == int and nz > -1:
            self.__nz = nz
        else:
            raise Exception('Invalid nz value. Exiting.')

    @nu.setter
    def nu(self, nu):
        if type(nu) == int and nu > -1:
            self.__nu = nu
        else:
            raise Exception('Invalid nu value. Exiting.')

    @np.setter
    def np(self, np):
        if type(np) == int and np > -1:
            self.__np = np
        else:
            raise Exception('Invalid np value. Exiting.')

    @ny.setter
    def ny(self, ny):
        if type(ny) == int and ny > -1:
            self.__ny = ny
        else:
            raise Exception('Invalid ny value. Exiting.')

    @ny_e.setter
    def ny_e(self, ny_e):
        if type(ny_e) == int and ny_e > -1:
            self.__ny_e = ny_e
        else:
            raise Exception('Invalid ny_e value. Exiting.')

    @nr.setter
    def nr(self, nr):
        if type(nr) == int and nr > -1:
            self.__nr = nr
        else:
            raise Exception('Invalid nr value. Exiting.')

    @nr_e.setter
    def nr_e(self, nr_e):
        if type(nr_e) == int and nr_e > -1:
            self.__nr_e = nr_e
        else:
            raise Exception('Invalid nr_e value. Exiting.')

    @nh.setter
    def nh(self, nh):
        if type(nh) == int and nh > -1:
            self.__nh = nh
        else:
            raise Exception('Invalid nh value. Exiting.')

    @nh_e.setter
    def nh_e(self, nh_e):
        if type(nh_e) == int and nh_e > -1:
            self.__nh_e = nh_e
        else:
            raise Exception('Invalid nh_e value. Exiting.')

    @nphi.setter
    def nphi(self, nphi):
        if type(nphi) == int and nphi > -1:
            self.__nphi = nphi
        else:
            raise Exception('Invalid nphi value. Exiting.')

    @nphi_e.setter
    def nphi_e(self, nphi_e):
        if type(nphi_e) == int and nphi_e > -1:
            self.__nphi_e = nphi_e
        else:
            raise Exception('Invalid nphi_e value. Exiting.')

    @nbx.setter
    def nbx(self, nbx):
        if type(nbx) == int and nbx > -1:
            self.__nbx = nbx
        else:
            raise Exception('Invalid nbx value. Exiting.')

    @nbx_0.setter
    def nbx_0(self, nbx_0):
        if type(nbx_0) == int and nbx_0 > -1:
            self.__nbx_0 = nbx_0
        else:
            raise Exception('Invalid nbx_0 value. Exiting.')

    @nbx_e.setter
    def nbx_e(self, nbx_e):
        if type(nbx_e) == int and nbx_e > -1:
            self.__nbx_e = nbx_e
        else:
            raise Exception('Invalid nbx_e value. Exiting.')

    @nbu.setter
    def nbu(self, nbu):
        if type(nbu) == int and nbu > -1:
            self.__nbu = nbu
        else:
            raise Exception('Invalid nbu value. Exiting.')

    @nsbx.setter
    def nsbx(self, nsbx):
        if type(nsbx) == int and nsbx > -1:
            self.__nsbx = nsbx
        else:
            raise Exception('Invalid nsbx value. Exiting.')

    @nsbx_e.setter
    def nsbx_e(self, nsbx_e):
        if type(nsbx_e) == int and nsbx_e > -1:
            self.__nsbx_e = nsbx_e
        else:
            raise Exception('Invalid nsbx_e value. Exiting.')

    @nsbu.setter
    def nsbu(self, nsbu):
        if type(nsbu) == int and nsbu > -1:
            self.__nsbu = nsbu
        else:
            raise Exception('Invalid nsbu value. Exiting.')

    @nsh.setter
    def nsh(self, nsh):
        if type(nsh) == int and nsh > -1:
            self.__nsh = nsh
        else:
            raise Exception('Invalid nsh value. Exiting.')

    @nsh_e.setter
    def nsh_e(self, nsh_e):
        if type(nsh_e) == int and nsh_e > -1:
            self.__nsh_e = nsh_e
        else:
            raise Exception('Invalid nsh_e value. Exiting.')

    @nsphi.setter
    def nsphi(self, nsphi):
        if type(nsphi) == int and nsphi > -1:
            self.__nsphi = nsphi
        else:
            raise Exception('Invalid nsphi value. Exiting.')

    @nsphi_e.setter
    def nsphi_e(self, nsphi_e):
        if type(nsphi_e) == int and nsphi_e > -1:
            self.__nsphi_e = nsphi_e
        else:
            raise Exception('Invalid nsphi_e value. Exiting.')

    @ns.setter
    def ns(self, ns):
        if type(ns) == int and ns > -1:
            self.__ns = ns
        else:
            raise Exception('Invalid ns value. Exiting.')

    @ns_e.setter
    def ns_e(self, ns_e):
        if type(ns_e) == int and ns_e > -1:
            self.__ns_e = ns_e
        else:
            raise Exception('Invalid ns_e value. Exiting.')

    @ng.setter
    def ng(self, ng):
        if type(ng) == int and ng > -1:
            self.__ng = ng
        else:
            raise Exception('Invalid ng value. Exiting.')

    @ng_e.setter
    def ng_e(self, ng_e):
        if type(ng_e) == int and ng_e > -1:
            self.__ng_e = ng_e
        else:
            raise Exception('Invalid ng_e value. Exiting.')

    @N.setter
    def N(self, N):
        if type(N) == int and N > 0:
            self.__N = N
        else:
            raise Exception('Invalid N value. Exiting.')

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
        self.__W           = []           # weight matrix
        self.__Vx          = []           # x matrix coefficient
        self.__Vu          = []           # u matrix coefficient
        self.__Vz          = []           # z matrix coefficient
        self.__yref        = []           # reference
        self.__Zl          = []           # diagonal of Hessian wrt lower slack
        self.__Zu          = []           # diagonal of Hessian wrt upper slack
        self.__zl          = []           # gradient wrt lower slack
        self.__zu          = []           # gradient wrt upper slack
        # Mayer term
        self.__cost_type_e = 'LINEAR_LS'  # cost type for Mayer term
        self.__W_e         = []           # weight matrix for Mayer term
        self.__Vx_e        = []           # x matrix coefficient for Mayer term
        self.__yref_e      = []           # reference for Mayer term
        self.__Zl_e        = []           # diagonal of Hessian wrt lower slack for Mayer term
        self.__Zu_e        = []           # diagonal of Hessian wrt upper slack for Mayer term
        self.__zl_e        = []           # gradient wrt lower slack for Mayer term
        self.__zu_e        = []           # gradient wrt upper slack for Mayer term

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
        if type(W) == np.ndarray:
            self.__W = W
        else:
            raise Exception('Invalid W value. Exiting.')

    @Vx.setter
    def Vx(self, Vx):
        if type(Vx) == np.ndarray:
            self.__Vx = Vx
        else:
            raise Exception('Invalid Vx value. Exiting.')

    @Vu.setter
    def Vu(self, Vu):
        if type(Vu) == np.ndarray:
            self.__Vu = Vu
        else:
            raise Exception('Invalid Vu value. Exiting.')

    @Vz.setter
    def Vz(self, Vz):
        if type(Vz) == np.ndarray:
            self.__Vz = Vz
        else:
            raise Exception('Invalid Vz value. Exiting.')

    @yref.setter
    def yref(self, yref):
        if type(yref) == np.ndarray:
            self.__yref = yref
        else:
            raise Exception('Invalid yref value. Exiting.')

    @Zl.setter
    def Zl(self, Zl):
        if type(Zl) == np.ndarray:
            self.__Zl = Zl
        else:
            raise Exception('Invalid Zl value. Exiting.')

    @Zu.setter
    def Zu(self, Zu):
        if type(Zu) == np.ndarray:
            self.__Zu = Zu
        else:
            raise Exception('Invalid Zu value. Exiting.')

    @zl.setter
    def zl(self, zl):
        if type(zl) == np.ndarray:
            self.__zl = zl
        else:
            raise Exception('Invalid zl value. Exiting.')

    @zu.setter
    def zu(self, zu):
        if type(zu) == np.ndarray:
            self.__zu = zu
        else:
            raise Exception('Invalid zu value. Exiting.')

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
        if type(W_e) == np.ndarray:
            self.__W_e = W_e
        else:
            raise Exception('Invalid W_e value. Exiting.')

    @Vx_e.setter
    def Vx_e(self, Vx_e):
        if type(Vx_e) == np.ndarray:
            self.__Vx_e = Vx_e
        else:
            raise Exception('Invalid Vx_e value. Exiting.')

    @yref_e.setter
    def yref_e(self, yref_e):
        if type(yref_e) == np.ndarray:
            self.__yref_e = yref_e
        else:
            raise Exception('Invalid yref_e value. Exiting.')

    @Zl_e.setter
    def Zl_e(self, Zl_e):
        if type(Zl_e) == np.ndarray:
            self.__Zl_e = Zl_e
        else:
            raise Exception('Invalid Zl_e value. Exiting.')

    @Zu_e.setter
    def Zu_e(self, Zu_e):
        if type(Zu_e) == np.ndarray:
            self.__Zu_e = Zu_e
        else:
            raise Exception('Invalid Zu_e value. Exiting.')

    @zl_e.setter
    def zl_e(self, zl_e):
        if type(zl_e) == np.ndarray:
            self.__zl_e = zl_e
        else:
            raise Exception('Invalid zl_e value. Exiting.')

    @zu_e.setter
    def zu_e(self, zu_e):
        if type(zu_e) == np.ndarray:
            self.__zu_e = zu_e
        else:
            raise Exception('Invalid zu_e value. Exiting.')

    def set(self, attr, value):
        setattr(self, attr, value)

def print_J_to_idx_note():
    print("NOTE: J* matrix is converted to zero based vector idx* vector, which is returned here.")


class AcadosOcpConstraints:
    """
    class containing the description of the constraints
    """
    def __init__(self):
        self.__constr_type   = 'BGH'                  # constraint type
        self.__constr_type_e = 'BGH'                  # constraint type
        # initial x
        self.__lbx_0   = []                           # lower bounds on x for initial state
        self.__ubx_0   = []                           # upper bounds on x for initial state
        self.__idxbx_0 = []                           # indexes for bounds on x0
        # state bounds
        self.__lbx     = []                           # lower bounds on x
        self.__ubx     = []                           # upper bounds on x
        self.__idxbx   = []
        # bounds on x at t=T
        self.__lbx_e   = []                           # lower bounds on x at t=T
        self.__ubx_e   = []                           # upper bounds on x at t=T
        self.__idxbx_e = []
        # bounds on u
        self.__lbu     = []                           # lower bounds on u
        self.__ubu     = []                           # upper bounds on u
        self.__idxbu   = []
        # polytopic constraints
        self.__lg      = []                           # lower bound for general polytopic inequalities
        self.__ug      = []                           # upper bound for general polytopic inequalities
        self.__D       = []                           # D matrix in lg <= D * u + C * x <= ug
        self.__C       = []                           # C matrix in lg <= D * u + C * x <= ug
        # polytopic constraints at t=T
        self.__C_e     = []                           # C matrix at t=T
        self.__lg_e    = []                           # lower bound on general polytopic inequalities at t=T
        self.__ug_e    = []                           # upper bound on general polytopic inequalities at t=T
        # nonlinear constraints
        self.__lh      = []                           # lower bound for nonlinear inequalities
        self.__uh      = []                           # upper bound for nonlinear inequalities
        # nonlinear constraints at t=T
        self.__uh_e    = []                           # upper bound on nonlinear inequalities at t=T
        self.__lh_e    = []                           # lower bound on nonlinear inequalities at t=T
        # convex-over-nonlinear constraints
        self.__lphi    = []                           # lower bound for convex-over-nonlinear inequalities
        self.__uphi    = []                           # upper bound for convex-over-nonlinear inequalities
        # nonlinear constraints at t=T
        self.__uphi_e    = []                         # upper bound on convex-over-nonlinear inequalities at t=T
        self.__lphi_e    = []                         # lower bound on convex-over-nonlinear inequalities at t=T
        # SLACK BOUNDS
        # soft bounds on x
        self.__lsbx   = []                            # lower bounds on slacks corresponding to soft lower bounds on x
        self.__usbx   = []                            # lower bounds on slacks corresponding to soft upper bounds on x
        self.__idxsbx = []                            # indexes of soft bounds on x within the indices of bounds on x
        # soft bounds on u
        self.__lsbu   = []                            # lower bounds on slacks corresponding to soft lower bounds on u
        self.__usbu   = []                            # lower bounds on slacks corresponding to soft upper bounds on u
        self.__idxsbu = []                            # indexes of soft bounds on u within the indices of bounds on u
        # soft bounds on x at t=T
        self.__lsbx_e  = []                           # lower bounds on slacks corresponding to soft lower bounds on x at t=T
        self.__usbx_e  = []                           # lower bounds on slacks corresponding to soft upper bounds on x at t=T
        self.__idxsbx_e= []                           # indexes of soft bounds on x at t=T, within the indices of bounds on x at t=T
        # soft bounds on nonlinear constraints
        self.__lsh    = []                            # lower bounds on slacks corresponding to soft lower bounds for nonlinear constraints
        self.__ush    = []                            # lower bounds on slacks corresponding to soft upper bounds for nonlinear constraints
        self.__idxsh  = []                            # indexes of soft nonlinear constraints within the indices of nonlinear constraints
        # soft bounds on nonlinear constraints
        self.__lsphi  = []                            # lower bounds on slacks corresponding to soft lower bounds for convex-over-nonlinear constraints
        self.__usphi  = []                            # lower bounds on slacks corresponding to soft upper bounds for convex-over-nonlinear constraints
        self.__idxsphi  = []                          # indexes of soft convex-over-nonlinear constraints within the indices of nonlinear constraints
        # soft bounds on nonlinear constraints at t=T
        self.__lsh_e    = []                          # lower bounds on slacks corresponding to soft lower bounds for nonlinear constraints at t=T
        self.__ush_e    = []                          # lower bounds on slacks corresponding to soft upper bounds for nonlinear constraints at t=T
        self.__idxsh_e  = []                          # indexes of soft nonlinear constraints at t=T within the indices of nonlinear constraints at t=T
        # soft bounds on nonlinear constraints at t=T
        self.__lsphi_e    = []                        # lower bounds on slacks corresponding to soft lower bounds for convex-over-nonlinear constraints at t=T
        self.__usphi_e    = []                        # lower bounds on slacks corresponding to soft upper bounds for convex-over-nonlinear constraints at t=T
        self.__idxsphi_e  = []                        # indexes of soft nonlinear constraints at t=T within the indices of nonlinear constraints at t=T


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
        """:math:`\\bar{h}^e` - upper bound on nonlinear inequalities at t=T"""
        return self.__lh_e

    @property
    def uh_e(self):
        """:math:`\\underline{h}^e` - lower bound on nonlinear inequalities at t=T"""
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
        if type(idxbx_0) == np.ndarray:
            self.__idxbx_0 = idxbx_0
        else:
            raise Exception('Invalid idxbx_0 value. Exiting.')

    @x0.setter
    def x0(self, x0):
        if type(x0) == np.ndarray:
            self.__lbx_0 = x0
            self.__ubx_0 = x0
            self.__idxbx_0 = np.arange(x0.size)
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
        if type(D) == np.ndarray:
            self.__D = D
        else:
            raise Exception('Invalid D value. Exiting.')

    @C.setter
    def C(self, C):
        if type(C) == np.ndarray:
            self.__C = C
        else:
            raise Exception('Invalid C value. Exiting.')

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
        if type(C_e) == np.ndarray:
            self.__C_e = C_e
        else:
            raise Exception('Invalid C_e value. Exiting.')

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
        if type(Jsbx) == np.ndarray:
            self.__idxsbx = J_to_idx_slack(Jsbx)
        else:
            raise Exception('Invalid Jsbx value. Exiting.')

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
        if type(Jsphi) == np.ndarray:
            self.__Jsphi = Jsphi
            self.__idxsphi = J_to_idx_slack(Jsphi)
        else:
            raise Exception('Invalid Jsphi value. Exiting.')

    # soft bounds on nonlinear constraints at t=T
    @lsh_e.setter
    def lsh_e(self, lsh_e):
        if type(lsh_e) == np.ndarray:
            self.__lsh_e = lsh_e
        else:
            raise Exception('Invalid lsh_e value. Exiting.')

    @ush_e.setter
    def ush_e(self, ush_e):
        if type(ush_e) == np.ndarray:
            self.__ush_e = ush_e
        else:
            raise Exception('Invalid ush_e value. Exiting.')

    @idxsh_e.setter
    def idxsh_e(self, idxsh_e):
        if type(idxsh_e) == np.ndarray:
            self.__idxsh_e = idxsh_e
        else:
            raise Exception('Invalid idxsh_e value. Exiting.')

    # soft bounds on convex-over-nonlinear constraints at t=T
    @lsphi_e.setter
    def lsphi_e(self, lsphi_e):
        if type(lsphi_e) == np.ndarray:
            self.__lsphi_e = lsphi_e
        else:
            raise Exception('Invalid lsphi_e value. Exiting.')

    @usphi_e.setter
    def usphi_e(self, usphi_e):
        if type(usphi_e) == np.ndarray:
            self.__usphi_e = usphi_e
        else:
            raise Exception('Invalid usphi_e value. Exiting.')

    @idxsphi_e.setter
    def idxsphi_e(self, idxsphi_e):
        if type(idxsphi_e) == np.ndarray:
            self.__idxsphi_e = idxsphi_e
        else:
            raise Exception('Invalid idxsphi_e value. Exiting.')

    @Jsphi_e.setter
    def Jsphi_e(self, Jsphi_e):
        if type(Jsphi_e) == np.ndarray:
            self.__Jsphi_e = Jsphi_e
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
        self.__sim_method_num_stages  = 4                     # number of stages in the integrator
        self.__sim_method_num_steps   = 1                     # number of steps in the integrator
        self.__sim_method_newton_iter = 3                     # number of Newton iterations in simulation method
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
        self.__print_level = 0                                # print level (possible values: 0, 1)
        self.__model_external_shared_lib_dir   = None         # path to the the .so lib
        self.__model_external_shared_lib_name  = None         # name of the the .so lib
        self.__regularize_method = None


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
        return max([self.__nlp_solver_tol_eq, self.__nlp_solver_tol_ineq, self.__nlp_solver_tol_comp, self.__nlp_solver_tol_stat])

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

        self.__parameter_values = []

    @property
    def parameter_values(self):
        """:math:`p` - initial values for parameter"""
        return self.__parameter_values

    @parameter_values.setter
    def parameter_values(self, parameter_values):
        if type(parameter_values) == np.ndarray:
            self.__parameter_values = parameter_values
        else:
            raise Exception('Invalid parameter_values value. Exiting.')

    def set(self, attr, value):
        # tokenize string
        tokens = attr.split('_', 1)
        if len(tokens) > 1:
            setter_to_call = getattr(getattr(self, tokens[0]), 'set')
        else:
            setter_to_call = getattr(self, 'set')

        setter_to_call(tokens[1], value)

        return
