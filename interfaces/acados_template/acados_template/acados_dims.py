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

def check_int_value(name, value, *, positive=False, nonnegative=False):
    if not isinstance(value, int):
        raise TypeError(f"Invalid {name} value: expected an integer, got {type(value).__name__}.")
    if positive and value <= 0:
        raise ValueError(f"Invalid {name} value: expected a positive integer, got {value}.")
    if nonnegative and value < 0:
        raise ValueError(f"Invalid {name} value: expected a nonnegative integer, got {value}.")


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
        check_int_value("nx", nx, positive=True)
        self.__nx = nx

    @nz.setter
    def nz(self, nz):
        check_int_value("nz", nz, nonnegative=True)
        self.__nz = nz

    @nu.setter
    def nu(self, nu):
        check_int_value("nu", nu, nonnegative=True)
        self.__nu = nu

    @np.setter
    def np(self, np):
        check_int_value("np", np, nonnegative=True)
        self.__np = np



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
        self.__nx_next = None
        # cost
        self.__ny = 0
        self.__ny_0 = 0
        self.__ny_e = 0
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
        self.__nsh_0 = 0
        self.__nsh_e = 0
        self.__nsphi = 0
        self.__nsphi_0 = 0
        self.__nsphi_e = 0
        self.__ns = 0
        self.__ns_0 = 0
        self.__ns_e = 0
        self.__nsg = 0
        self.__nsg_e = 0
        # equalities within x bounds
        self.__nbxe_0 = None
        # global parameters
        self.__np_global = 0
        self.__n_global_data = 0


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
    def nx_next(self):
        r""":math:`n_{x, \text{next}}` - state dimension of next state.
        Type: int; default: None"""
        return self.__nx_next

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
        r""":math:`n_{\pi}` - dimension of the image of the inner nonlinear function in positive definite constraints.
        Type: int; default: 0"""
        return self.__nr

    @property
    def nr_e(self):
        r""":math:`n_{\pi}^e` - dimension of the image of the inner nonlinear function in positive definite constraints.
        Type: int; default: 0"""
        return self.__nr_e

    @property
    def nr_0(self):
        r""":math:`n_{\pi}^0` - dimension of the image of the inner nonlinear function in positive definite constraints.
        Type: int; default: 0"""
        return self.__nr_0

    @property
    def nphi(self):
        r""":math:`n_{\phi}` - number of convex-over-nonlinear constraints.
        Type: int; default: 0"""
        return self.__nphi

    @property
    def nphi_0(self):
        r""":math:`n_{\phi}^0` - number of convex-over-nonlinear constraints at initial shooting node 0.
        Type: int; default: 0"""
        return self.__nphi_0

    @property
    def nphi_e(self):
        r""":math:`n_{\phi}^e` - number of convex-over-nonlinear constraints at terminal shooting node N.
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
        r""":math:`n_{{s\phi}^0}` - number of soft convex-over-nonlinear constraints at shooting node 0.
        Type: int; default: 0"""
        return self.__nsphi_0

    @property
    def nsphi(self):
        r""":math:`n_{{s\phi}}` - number of soft convex-over-nonlinear constraints.
        Type: int; default: 0"""
        return self.__nsphi

    @property
    def nsphi_e(self):
        r""":math:`n_{{s\phi}^e}` - number of soft convex-over-nonlinear constraints at terminal shooting node N.
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
    def np_global(self):
        """number of global parameters p_global; default: 0"""
        return self.__np_global

    @property
    def n_global_data(self):
        """size of global_data; expressions that only depend on p_global; detected automatically during code generation.
        Type: int; default: 0"""
        return self.__n_global_data

    @property
    def N(self):
        """
        :math:`N` - Number of shooting intervals.
        DEPRECATED: use ocp.solver_options.N instead.

        Type: int; default: None"""
        return self.__N

    @nx.setter
    def nx(self, nx):
        check_int_value("nx", nx, positive=True)
        self.__nx = nx

    @nz.setter
    def nz(self, nz):
        check_int_value("nz", nz, nonnegative=True)
        self.__nz = nz

    @nu.setter
    def nu(self, nu):
        check_int_value("nu", nu, nonnegative=True)
        self.__nu = nu

    @np.setter
    def np(self, np):
        check_int_value("np", np, nonnegative=True)
        self.__np = np

    @nx_next.setter
    def nx_next(self, nx_next):
        check_int_value("nx_next", nx_next, positive=True)
        self.__nx_next = nx_next

    @np_global.setter
    def np_global(self, np_global):
        check_int_value("np_global", np_global, nonnegative=True)
        self.__np_global = np_global

    @n_global_data.setter
    def n_global_data(self, n_global_data):
        check_int_value("n_global_data", n_global_data, nonnegative=True)
        self.__n_global_data = n_global_data

    @ny_0.setter
    def ny_0(self, ny_0):
        check_int_value("ny_0", ny_0, nonnegative=True)
        self.__ny_0 = ny_0

    @ny.setter
    def ny(self, ny):
        check_int_value("ny", ny, nonnegative=True)
        self.__ny = ny

    @ny_e.setter
    def ny_e(self, ny_e):
        check_int_value("ny_e", ny_e, nonnegative=True)
        self.__ny_e = ny_e

    @nh.setter
    def nh(self, nh):
        check_int_value("nh", nh, nonnegative=True)
        self.__nh = nh

    @nh_0.setter
    def nh_0(self, nh_0):
        check_int_value("nh_0", nh_0, nonnegative=True)
        self.__nh_0 = nh_0

    @nh_e.setter
    def nh_e(self, nh_e):
        check_int_value("nh_e", nh_e, nonnegative=True)
        self.__nh_e = nh_e

    @nphi_0.setter
    def nphi_0(self, nphi_0):
        check_int_value("nphi_0", nphi_0, nonnegative=True)
        self.__nphi_0 = nphi_0

    @nphi.setter
    def nphi(self, nphi):
        check_int_value("nphi", nphi, nonnegative=True)
        self.__nphi = nphi

    @nphi_e.setter
    def nphi_e(self, nphi_e):
        check_int_value("nphi_e", nphi_e, nonnegative=True)
        self.__nphi_e = nphi_e

    @nr_0.setter
    def nr_0(self, nr_0):
        check_int_value("nr_0", nr_0, nonnegative=True)
        self.__nr_0 = nr_0

    @nr.setter
    def nr(self, nr):
        check_int_value("nr", nr, nonnegative=True)
        self.__nr = nr

    @nr_e.setter
    def nr_e(self, nr_e):
        check_int_value("nr_e", nr_e, nonnegative=True)
        self.__nr_e = nr_e

    @nbx.setter
    def nbx(self, nbx):
        check_int_value("nbx", nbx, nonnegative=True)
        self.__nbx = nbx

    @nbxe_0.setter
    def nbxe_0(self, nbxe_0):
        check_int_value("nbxe_0", nbxe_0, nonnegative=True)
        self.__nbxe_0 = nbxe_0

    @nbx_0.setter
    def nbx_0(self, nbx_0):
        check_int_value("nbx_0", nbx_0, nonnegative=True)
        self.__nbx_0 = nbx_0

    @nbx_e.setter
    def nbx_e(self, nbx_e):
        check_int_value("nbx_e", nbx_e, nonnegative=True)
        self.__nbx_e = nbx_e

    @nbu.setter
    def nbu(self, nbu):
        check_int_value("nbu", nbu, nonnegative=True)
        self.__nbu = nbu

    @nsbx.setter
    def nsbx(self, nsbx):
        check_int_value("nsbx", nsbx, nonnegative=True)
        self.__nsbx = nsbx

    @nsbx_e.setter
    def nsbx_e(self, nsbx_e):
        check_int_value("nsbx_e", nsbx_e, nonnegative=True)
        self.__nsbx_e = nsbx_e

    @nsbu.setter
    def nsbu(self, nsbu):
        check_int_value("nsbu", nsbu, nonnegative=True)
        self.__nsbu = nsbu

    @nsg.setter
    def nsg(self, nsg):
        check_int_value("nsg", nsg, nonnegative=True)
        self.__nsg = nsg

    @nsg_e.setter
    def nsg_e(self, nsg_e):
        check_int_value("nsg_e", nsg_e, nonnegative=True)
        self.__nsg_e = nsg_e

    @nsh_0.setter
    def nsh_0(self, nsh_0):
        check_int_value("nsh_0", nsh_0, nonnegative=True)
        self.__nsh_0 = nsh_0

    @nsh.setter
    def nsh(self, nsh):
        check_int_value("nsh", nsh, nonnegative=True)
        self.__nsh = nsh

    @nsh_e.setter
    def nsh_e(self, nsh_e):
        check_int_value("nsh_e", nsh_e, nonnegative=True)
        self.__nsh_e = nsh_e

    @nsphi_0.setter
    def nsphi_0(self, nsphi_0):
        check_int_value("nsphi_0", nsphi_0, nonnegative=True)
        self.__nsphi_0 = nsphi_0

    @nsphi.setter
    def nsphi(self, nsphi):
        check_int_value("nsphi", nsphi, nonnegative=True)
        self.__nsphi = nsphi

    @nsphi_e.setter
    def nsphi_e(self, nsphi_e):
        check_int_value("nsphi_e", nsphi_e, nonnegative=True)
        self.__nsphi_e = nsphi_e

    @ns_0.setter
    def ns_0(self, ns_0):
        check_int_value("ns_0", ns_0, nonnegative=True)
        self.__ns_0 = ns_0

    @ns.setter
    def ns(self, ns):
        check_int_value("ns", ns, nonnegative=True)
        self.__ns = ns

    @ns_e.setter
    def ns_e(self, ns_e):
        check_int_value("ns_e", ns_e, nonnegative=True)
        self.__ns_e = ns_e

    @ng.setter
    def ng(self, ng):
        check_int_value("ng", ng, nonnegative=True)
        self.__ng = ng

    @ng_e.setter
    def ng_e(self, ng_e):
        check_int_value("ng_e", ng_e, nonnegative=True)
        self.__ng_e = ng_e

    @N.setter
    def N(self, N):
        check_int_value("N", N, nonnegative=True)
        self.__N = N
