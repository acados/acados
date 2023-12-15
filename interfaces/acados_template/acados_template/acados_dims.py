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

class AcadosSimDims:
    """
    Class containing the dimensions of the model to be simulated.
    """
    def __init__(self):
        self.__nx = None
        self.__nu = None
        self.__nz = 0
        self.__np = 0

    @property
    def nx(self):
        """:math:`n_x` - number of states. Type: int > 0"""
        return self.__nx

    @property
    def nz(self):
        """:math:`n_z` - number of algebraic variables. Type: int >= 0"""
        return self.__nz

    @property
    def nu(self):
        """:math:`n_u` - number of inputs. Type: int >= 0"""
        return self.__nu

    @property
    def np(self):
        """:math:`n_p` - number of parameters. Type: int >= 0"""
        return self.__np

    @nx.setter
    def nx(self, nx):
        if isinstance(nx, int) and nx > 0:
            self.__nx = nx
        else:
            raise Exception('Invalid nx value, expected positive integer.')

    @nz.setter
    def nz(self, nz):
        if isinstance(nz, int) and nz > -1:
            self.__nz = nz
        else:
            raise Exception('Invalid nz value, expected nonnegative integer.')

    @nu.setter
    def nu(self, nu):
        if isinstance(nu, int) and nu > -1:
            self.__nu = nu
        else:
            raise Exception('Invalid nu value, expected nonnegative integer.')

    @np.setter
    def np(self, np):
        if isinstance(np, int) and np > -1:
            self.__np = np
        else:
            raise Exception('Invalid np value, expected nonnegative integer.')


class AcadosOcpDims:
    """
    Class containing the dimensions of the optimal control problem.
    """
    def __init__(self):
        # number of shooting intervals
        self.__N = None
        # model
        self.__nx = None
        self.__nu = None
        self.__nz = 0
        self.__np = 0
        # cost
        self.__ny = 0
        self.__ny_e = 0
        self.__ny_0 = 0
        # bounds
        self.__nbu = 0
        self.__nbx = 0
        self.__nbx_0 = 0
        self.__nbx_e = 0
        # nonlinear constraints
        self.__nh = 0
        self.__nh_0 = 0
        self.__nh_e = 0
        self.__nr = 0
        self.__nr_0 = 0
        self.__nr_e = 0
        self.__nphi = 0
        self.__nphi_0 = 0
        self.__nphi_e = 0
        # general linear constraints
        self.__ng = 0
        self.__ng_e = 0
        # soft constraints
        self.__nsbx = 0
        self.__nsbx_e = 0
        self.__nsbu = 0
        self.__nsh = 0
        self.__nsh_e = 0
        self.__nsh_0 = 0
        self.__nsphi = 0
        self.__nsphi_e = 0
        self.__nsphi_0 = 0
        self.__ns_0 = 0
        self.__ns = 0
        self.__ns_e = 0
        self.__nsg = 0
        self.__nsg_e = 0
        # equalities within x bounds
        self.__nbxe_0 = None


    @property
    def nx(self):
        """:math:`n_x` - number of states.
        Type: int; default: None"""
        return self.__nx

    @property
    def nz(self):
        """:math:`n_z` - number of algebraic variables.
        Type: int; default: 0"""
        return self.__nz

    @property
    def nu(self):
        """:math:`n_u` - number of inputs.
        Type: int; default: None"""
        return self.__nu

    @property
    def np(self):
        """:math:`n_p` - number of parameters.
        Type: int; default: 0"""
        return self.__np

    @property
    def ny(self):
        """:math:`n_y` - number of residuals in Lagrange term.
        Type: int; default: 0"""
        return self.__ny

    @property
    def ny_0(self):
        """:math:`n_{y}^0` - number of residuals in Mayer term.
        Type: int; default: 0"""
        return self.__ny_0

    @property
    def ny_e(self):
        """:math:`n_{y}^e` - number of residuals in Mayer term.
        Type: int; default: 0"""
        return self.__ny_e

    @property
    def nh(self):
        """:math:`n_h` - number of nonlinear constraints.
        Type: int; default: 0"""
        return self.__nh

    @property
    def nh_0(self):
        """:math:`n_{h}^e` - number of nonlinear constraints at initial shooting node.
        Type: int; default: 0"""
        return self.__nh_0

    @property
    def nh_e(self):
        """:math:`n_{h}^e` - number of nonlinear constraints at terminal shooting node N.
        Type: int; default: 0"""
        return self.__nh_e

    @property
    def nr(self):
        """:math:`n_{\pi}` - dimension of the image of the inner nonlinear function in positive definite constraints.
        Type: int; default: 0"""
        return self.__nr

    @property
    def nr_e(self):
        """:math:`n_{\pi}^e` - dimension of the image of the inner nonlinear function in positive definite constraints.
        Type: int; default: 0"""
        return self.__nr_e

    @property
    def nr_0(self):
        """:math:`n_{\pi}^0` - dimension of the image of the inner nonlinear function in positive definite constraints.
        Type: int; default: 0"""
        return self.__nr_0

    @property
    def nphi(self):
        """:math:`n_{\phi}` - number of convex-over-nonlinear constraints.
        Type: int; default: 0"""
        return self.__nphi

    @property
    def nphi_0(self):
        """:math:`n_{\phi}^0` - number of convex-over-nonlinear constraints at initial shooting node 0.
        Type: int; default: 0"""
        return self.__nphi_0

    @property
    def nphi_e(self):
        """:math:`n_{\phi}^e` - number of convex-over-nonlinear constraints at terminal shooting node N.
        Type: int; default: 0"""
        return self.__nphi_e

    @property
    def nbx(self):
        """:math:`n_{b_x}` - number of state bounds.
        Type: int; default: 0"""
        return self.__nbx

    @property
    def nbxe_0(self):
        """:math:`n_{be_{x0}}` - number of state bounds at initial shooting node that are equalities.
        Type: int; default: None"""
        return self.__nbxe_0

    @property
    def nbx_0(self):
        """:math:`n_{b_{x0}}` - number of state bounds for initial state.
        Type: int; default: 0"""
        return self.__nbx_0

    @property
    def nbx_e(self):
        """:math:`n_{b_x}` - number of state bounds at terminal shooting node N.
        Type: int; default: 0"""
        return self.__nbx_e

    @property
    def nbu(self):
        """:math:`n_{b_u}` - number of input bounds.
        Type: int; default: 0"""
        return self.__nbu

    @property
    def nsbx(self):
        """:math:`n_{{sb}_x}` - number of soft state bounds.
        Type: int; default: 0"""
        return self.__nsbx

    @property
    def nsbx_e(self):
        """:math:`n_{{sb}^e_{x}}` - number of soft state bounds at terminal shooting node N.
        Type: int; default: 0"""
        return self.__nsbx_e

    @property
    def nsbu(self):
        """:math:`n_{{sb}_u}` - number of soft input bounds.
        Type: int; default: 0"""
        return self.__nsbu

    @property
    def nsg(self):
        """:math:`n_{{sg}}` - number of soft general linear constraints.
        Type: int; default: 0"""
        return self.__nsg

    @property
    def nsg_e(self):
        """:math:`n_{{sg}^e}` - number of soft general linear constraints at terminal shooting node N.
        Type: int; default: 0"""
        return self.__nsg_e

    @property
    def nsh_0(self):
        """:math:`n_{{sh}}^0` - number of soft nonlinear constraints at shooting node 0.
        Type: int; default: 0"""
        return self.__nsh_0

    @property
    def nsh(self):
        """:math:`n_{{sh}}` - number of soft nonlinear constraints.
        Type: int; default: 0"""
        return self.__nsh

    @property
    def nsh_e(self):
        """:math:`n_{{sh}}^e` - number of soft nonlinear constraints at terminal shooting node N.
        Type: int; default: 0"""
        return self.__nsh_e

    @property
    def nsphi_0(self):
        """:math:`n_{{s\phi}^0}` - number of soft convex-over-nonlinear constraints at shooting node 0.
        Type: int; default: 0"""
        return self.__nsphi_0

    @property
    def nsphi(self):
        """:math:`n_{{s\phi}}` - number of soft convex-over-nonlinear constraints.
        Type: int; default: 0"""
        return self.__nsphi

    @property
    def nsphi_e(self):
        """:math:`n_{{s\phi}^e}` - number of soft convex-over-nonlinear constraints at terminal shooting node N.
        Type: int; default: 0"""
        return self.__nsphi_e

    @property
    def ns_0(self):
        """:math:`n_{s}^0` - total number of slacks at shooting node 0.
        Type: int; default: 0"""
        return self.__ns_0

    @property
    def ns(self):
        """:math:`n_{s}` - total number of slacks at stages (1, N-1).
        Type: int; default: 0"""
        return self.__ns

    @property
    def ns_e(self):
        """:math:`n_{s}^e` - total number of slacks at terminal shooting node N.
        Type: int; default: 0"""
        return self.__ns_e

    @property
    def ng(self):
        """:math:`n_{g}` - number of general polytopic constraints.
        Type: int; default: 0"""
        return self.__ng

    @property
    def ng_e(self):
        """:math:`n_{g}^e` - number of general polytopic constraints at terminal shooting node N.
        Type: int; default: 0"""
        return self.__ng_e

    @property
    def N(self):
        """:math:`N` - prediction horizon.
        Type: int; default: None"""
        return self.__N

    @nx.setter
    def nx(self, nx):
        if isinstance(nx, int) and nx > 0:
            self.__nx = nx
        else:
            raise Exception('Invalid nx value, expected positive integer.')

    @nz.setter
    def nz(self, nz):
        if isinstance(nz, int) and nz > -1:
            self.__nz = nz
        else:
            raise Exception('Invalid nz value, expected nonnegative integer.')

    @nu.setter
    def nu(self, nu):
        if isinstance(nu, int) and nu > -1:
            self.__nu = nu
        else:
            raise Exception('Invalid nu value, expected nonnegative integer.')

    @np.setter
    def np(self, np):
        if isinstance(np, int) and np > -1:
            self.__np = np
        else:
            raise Exception('Invalid np value, expected nonnegative integer.')

    @ny_0.setter
    def ny_0(self, ny_0):
        if isinstance(ny_0, int) and ny_0 > -1:
            self.__ny_0 = ny_0
        else:
            raise Exception('Invalid ny_0 value, expected nonnegative integer.')

    @ny.setter
    def ny(self, ny):
        if isinstance(ny, int) and ny > -1:
            self.__ny = ny
        else:
            raise Exception('Invalid ny value, expected nonnegative integer.')

    @ny_e.setter
    def ny_e(self, ny_e):
        if isinstance(ny_e, int) and ny_e > -1:
            self.__ny_e = ny_e
        else:
            raise Exception('Invalid ny_e value, expected nonnegative integer.')

    @nh.setter
    def nh(self, nh):
        if isinstance(nh, int) and nh > -1:
            self.__nh = nh
        else:
            raise Exception('Invalid nh value, expected nonnegative integer.')

    @nh_0.setter
    def nh_0(self, nh_0):
        if isinstance(nh_0, int) and nh_0 > -1:
            self.__nh_0 = nh_0
        else:
            raise Exception('Invalid nh_0 value, expected nonnegative integer.')

    @nh_e.setter
    def nh_e(self, nh_e):
        if isinstance(nh_e, int) and nh_e > -1:
            self.__nh_e = nh_e
        else:
            raise Exception('Invalid nh_e value, expected nonnegative integer.')

    @nphi_0.setter
    def nphi_0(self, nphi_0):
        if isinstance(nphi_0, int) and nphi_0 > -1:
            self.__nphi_0 = nphi_0
        else:
            raise Exception('Invalid nphi_0 value, expected nonnegative integer.')

    @nphi.setter
    def nphi(self, nphi):
        if isinstance(nphi, int) and nphi > -1:
            self.__nphi = nphi
        else:
            raise Exception('Invalid nphi value, expected nonnegative integer.')

    @nphi_e.setter
    def nphi_e(self, nphi_e):
        if isinstance(nphi_e, int) and nphi_e > -1:
            self.__nphi_e = nphi_e
        else:
            raise Exception('Invalid nphi_e value, expected nonnegative integer.')

    @nr_0.setter
    def nr_0(self, nr_0):
        if isinstance(nr_0, int) and nr_0 > -1:
            self.__nr_0 = nr_0
        else:
            raise Exception('Invalid nr_0 value, expected nonnegative integer.')

    @nr.setter
    def nr(self, nr):
        if isinstance(nr, int) and nr > -1:
            self.__nr = nr
        else:
            raise Exception('Invalid nr value, expected nonnegative integer.')

    @nr_e.setter
    def nr_e(self, nr_e):
        if isinstance(nr_e, int) and nr_e > -1:
            self.__nr_e = nr_e
        else:
            raise Exception('Invalid nr_e value, expected nonnegative integer.')

    @nbx.setter
    def nbx(self, nbx):
        if isinstance(nbx, int) and nbx > -1:
            self.__nbx = nbx
        else:
            raise Exception('Invalid nbx value, expected nonnegative integer.')

    @nbxe_0.setter
    def nbxe_0(self, nbxe_0):
        if isinstance(nbxe_0, int) and nbxe_0 > -1:
            self.__nbxe_0 = nbxe_0
        else:
            raise Exception('Invalid nbxe_0 value, expected nonnegative integer.')

    @nbx_0.setter
    def nbx_0(self, nbx_0):
        if isinstance(nbx_0, int) and nbx_0 > -1:
            self.__nbx_0 = nbx_0
        else:
            raise Exception('Invalid nbx_0 value, expected nonnegative integer.')

    @nbx_e.setter
    def nbx_e(self, nbx_e):
        if isinstance(nbx_e, int) and nbx_e > -1:
            self.__nbx_e = nbx_e
        else:
            raise Exception('Invalid nbx_e value, expected nonnegative integer.')

    @nbu.setter
    def nbu(self, nbu):
        if isinstance(nbu, int) and nbu > -1:
            self.__nbu = nbu
        else:
            raise Exception('Invalid nbu value, expected nonnegative integer.')

    @nsbx.setter
    def nsbx(self, nsbx):
        if isinstance(nsbx, int) and nsbx > -1:
            self.__nsbx = nsbx
        else:
            raise Exception('Invalid nsbx value, expected nonnegative integer.')

    @nsbx_e.setter
    def nsbx_e(self, nsbx_e):
        if isinstance(nsbx_e, int) and nsbx_e > -1:
            self.__nsbx_e = nsbx_e
        else:
            raise Exception('Invalid nsbx_e value, expected nonnegative integer.')

    @nsbu.setter
    def nsbu(self, nsbu):
        if isinstance(nsbu, int) and nsbu > -1:
            self.__nsbu = nsbu
        else:
            raise Exception('Invalid nsbu value, expected nonnegative integer.')

    @nsg.setter
    def nsg(self, nsg):
        if isinstance(nsg, int) and nsg > -1:
            self.__nsg = nsg
        else:
            raise Exception('Invalid nsg value, expected nonnegative integer.')

    @nsg_e.setter
    def nsg_e(self, nsg_e):
        if isinstance(nsg_e, int) and nsg_e > -1:
            self.__nsg_e = nsg_e
        else:
            raise Exception('Invalid nsg_e value, expected nonnegative integer.')

    @nsh_0.setter
    def nsh_0(self, nsh_0):
        if isinstance(nsh_0, int) and nsh_0 > -1:
            self.__nsh_0 = nsh_0
        else:
            raise Exception('Invalid nsh_0 value, expected nonnegative integer.')

    @nsh.setter
    def nsh(self, nsh):
        if isinstance(nsh, int) and nsh > -1:
            self.__nsh = nsh
        else:
            raise Exception('Invalid nsh value, expected nonnegative integer.')

    @nsh_e.setter
    def nsh_e(self, nsh_e):
        if isinstance(nsh_e, int) and nsh_e > -1:
            self.__nsh_e = nsh_e
        else:
            raise Exception('Invalid nsh_e value, expected nonnegative integer.')

    @nsphi_0.setter
    def nsphi_0(self, nsphi_0):
        if isinstance(nsphi_0, int) and nsphi_0 > -1:
            self.__nsphi_0 = nsphi_0
        else:
            raise Exception('Invalid nsphi_0 value, expected nonnegative integer.')

    @nsphi.setter
    def nsphi(self, nsphi):
        if isinstance(nsphi, int) and nsphi > -1:
            self.__nsphi = nsphi
        else:
            raise Exception('Invalid nsphi value, expected nonnegative integer.')

    @nsphi_e.setter
    def nsphi_e(self, nsphi_e):
        if isinstance(nsphi_e, int) and nsphi_e > -1:
            self.__nsphi_e = nsphi_e
        else:
            raise Exception('Invalid nsphi_e value, expected nonnegative integer.')

    @ns_0.setter
    def ns_0(self, ns_0):
        if isinstance(ns_0, int) and ns_0 > -1:
            self.__ns_0 = ns_0
        else:
            raise Exception('Invalid ns_0 value, expected nonnegative integer.')

    @ns.setter
    def ns(self, ns):
        if isinstance(ns, int) and ns > -1:
            self.__ns = ns
        else:
            raise Exception('Invalid ns value, expected nonnegative integer.')

    @ns_e.setter
    def ns_e(self, ns_e):
        if isinstance(ns_e, int) and ns_e > -1:
            self.__ns_e = ns_e
        else:
            raise Exception('Invalid ns_e value, expected nonnegative integer.')

    @ng.setter
    def ng(self, ng):
        if isinstance(ng, int) and ng > -1:
            self.__ng = ng
        else:
            raise Exception('Invalid ng value, expected nonnegative integer.')

    @ng_e.setter
    def ng_e(self, ng_e):
        if isinstance(ng_e, int) and ng_e > -1:
            self.__ng_e = ng_e
        else:
            raise Exception('Invalid ng_e value, expected nonnegative integer.')

    @N.setter
    def N(self, N):
        if isinstance(N, int) and N > 0:
            self.__N = N
        else:
            raise Exception('Invalid N value, expected positive integer.')
