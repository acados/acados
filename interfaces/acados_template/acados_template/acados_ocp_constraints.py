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
from .utils import J_to_idx, print_J_to_idx_note, J_to_idx_slack

class AcadosOcpConstraints:
    """
    class containing the description of the constraints
    """
    def __init__(self):
        self.__constr_type_0 = 'BGH'
        self.__constr_type = 'BGH'
        self.__constr_type_e = 'BGH'
        # initial x
        self.__lbx_0   = np.array([])
        self.__ubx_0   = np.array([])
        self.__idxbx_0 = np.array([])
        self.__idxbxe_0 = np.array([])
        self.__has_x0 = False
        # state bounds
        self.__lbx     = np.array([])
        self.__ubx     = np.array([])
        self.__idxbx   = np.array([])
        # bounds on x at shooting node N
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
        # polytopic constraints at shooting node N
        self.__C_e     = np.zeros((0,0))
        self.__lg_e    = np.array([])
        self.__ug_e    = np.array([])
        # nonlinear constraints
        self.__lh      = np.array([])
        self.__uh      = np.array([])
        # nonlinear constraints at initial shooting node
        self.__uh_0    = np.array([])
        self.__lh_0    = np.array([])
        # nonlinear constraints at shooting node N
        self.__uh_e    = np.array([])
        self.__lh_e    = np.array([])
        # convex-over-nonlinear constraints
        self.__lphi    = np.array([])
        self.__uphi    = np.array([])
        self.__uphi_e = np.array([])
        self.__lphi_e = np.array([])
        self.__uphi_0 = np.array([])
        self.__lphi_0 = np.array([])
        # SLACK BOUNDS
        # soft bounds on x
        self.__lsbx   = np.array([])
        self.__usbx   = np.array([])
        self.__idxsbx = np.array([])
        # soft bounds on u
        self.__lsbu   = np.array([])
        self.__usbu   = np.array([])
        self.__idxsbu = np.array([])
        # soft bounds on x at shooting node N
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
        # soft bounds on general linear constraints at shooting node N
        self.__lsg_e    = np.array([])
        self.__usg_e    = np.array([])
        self.__idxsg_e  = np.array([])
        # soft bounds on nonlinear constraints at shooting node N
        self.__lsh_e    = np.array([])
        self.__ush_e    = np.array([])
        self.__idxsh_e  = np.array([])
        # soft bounds on nonlinear constraints at shooting node N
        self.__lsphi_e    = np.array([])
        self.__usphi_e    = np.array([])
        self.__idxsphi_e  = np.array([])
        # soft bounds on nonlinear constraints at shooting node 0
        self.__lsh_0    = np.array([])
        self.__ush_0    = np.array([])
        self.__idxsh_0  = np.array([])
        # soft bounds on nonlinear constraints at shooting node 0
        self.__lsphi_0    = np.array([])
        self.__usphi_0    = np.array([])
        self.__idxsphi_0  = np.array([])


    # types
    @property
    def constr_type(self):
        """Constraints type for shooting nodes (0 to N-1). string in {BGH, BGP}.
        Detected automatically from model.
        Default: BGH; BGP is for convex over nonlinear.
        """
        return self.__constr_type

    @property
    def constr_type_0(self):
        """Constraints type for initial shooting node. string in {BGH, BGP}.
        Detected automatically from model.
        Default: BGH; BGP is for convex over nonlinear.
        """
        return self.__constr_type_0

    @property
    def constr_type_e(self):
        """Constraints type for terminal shooting node N. string in {BGH, BGP}.
        Detected automatically from model.
        Default: BGH; BGP is for convex over nonlinear.
        """
        return self.__constr_type_e

    # initial bounds on x
    @property
    def lbx_0(self):
        """:math:`\\underline{x_0}` - lower bounds on x at initial stage 0.
        Type: :code:`np.ndarray`; default: :code:`np.array([])`."""
        return self.__lbx_0

    @property
    def ubx_0(self):
        """:math:`\\bar{x_0}` - upper bounds on x at initial stage 0.
        Type: :code:`np.ndarray`; default: :code:`np.array([])`"""
        return self.__ubx_0

    @property
    def Jbx_0(self):
        """:math:`J_{bx,0}` - matrix coefficient for bounds on x at initial stage 0.
        Translated internally to :py:attr:`idxbx_0`"""
        print_J_to_idx_note()
        return self.__idxbx_0

    @property
    def idxbx_0(self):
        """Indices of bounds on x at initial stage 0
        -- can be set automatically via x0.
        Can be set by using :py:attr:`Jbx_0`.
        Type: :code:`np.ndarray`; default: :code:`np.array([])`"""
        return self.__idxbx_0

    @property
    def idxbxe_0(self):
        """Indices of bounds on x0 that are equalities -- can be set automatically via :py:attr:`x0`.
        Type: :code:`np.ndarray`; default: :code:`np.array([])`"""
        return self.__idxbxe_0

    # bounds on x
    @property
    def lbx(self):
        """:math:`\\underline{x}` - lower bounds on x at intermediate shooting nodes (1 to N-1).
        Type: :code:`np.ndarray`; default: :code:`np.array([])`"""
        return self.__lbx

    @property
    def ubx(self):
        """:math:`\\bar{x}` - upper bounds on x at intermediate shooting nodes (1 to N-1).
        Type: :code:`np.ndarray`; default: :code:`np.array([])`"""
        return self.__ubx

    @property
    def idxbx(self):
        """indices of bounds on x (defines :math:`J_{bx}`) at intermediate shooting nodes (1 to N-1).
        Can be set by using :py:attr:`Jbx`.
        Type: :code:`np.ndarray`; default: :code:`np.array([])`"""
        return self.__idxbx

    @property
    def Jbx(self):
        """:math:`J_{bx}` - matrix coefficient for bounds on x
        at intermediate shooting nodes (1 to N-1).
        Translated internally into :py:attr:`idxbx`."""
        print_J_to_idx_note()
        return self.__idxbx

    # bounds on x at shooting node N
    @property
    def lbx_e(self):
        """:math:`\\underline{x}^e` - lower bounds on x at terminal shooting node N.
        Type: :code:`np.ndarray`; default: :code:`np.array([])`"""
        return self.__lbx_e

    @property
    def ubx_e(self):
        """:math:`\\bar{x}^e` - upper bounds on x at terminal shooting node N.
        Type: :code:`np.ndarray`; default: :code:`np.array([])`"""
        return self.__ubx_e

    @property
    def idxbx_e(self):
        """Indices for bounds on x at terminal shooting node N (defines :math:`J_{bx}^e`).
        Can be set by using :py:attr:`Jbx_e`.
        Type: :code:`np.ndarray`; default: :code:`np.array([])`"""
        return self.__idxbx_e

    @property
    def Jbx_e(self):
        """:math:`J_{bx}^e` matrix coefficient for bounds on x at terminal shooting node N.
        Translated internally into :py:attr:`idxbx_e`."""
        print_J_to_idx_note()
        return self.__idxbx_e

    # bounds on u
    @property
    def lbu(self):
        """:math:`\\underline{u}` - lower bounds on u at shooting nodes (0 to N-1).
        Type: :code:`np.ndarray`; default: :code:`np.array([])`
        """
        return self.__lbu

    @property
    def ubu(self):
        """:math:`\\bar{u}` - upper bounds on u at shooting nodes (0 to N-1).
        Type: :code:`np.ndarray`; default: :code:`np.array([])`
        """
        return self.__ubu

    @property
    def idxbu(self):
        """Indices of bounds on u (defines :math:`J_{bu}`) at shooting nodes (0 to N-1).
        Can be set by using :py:attr:`Jbu`.
        Type: :code:`np.ndarray`; default: :code:`np.array([])`
        """
        return self.__idxbu

    @property
    def Jbu(self):
        """:math:`J_{bu}` - matrix coefficient for bounds on u at shooting nodes (0 to N-1).
        Translated internally to :py:attr:`idxbu`.
        """
        print_J_to_idx_note()
        return self.__idxbu

    # polytopic constraints
    @property
    def C(self):
        """:math:`C` - C matrix in :math:`\\underline{g} \\leq D \, u + C \, x \\leq \\bar{g}`
        at shooting nodes (0 to N-1).
        Type: :code:`np.ndarray`; default: :code:`np.array((0,0))`.
        """
        return self.__C

    @property
    def D(self):
        """:math:`D` - D matrix in :math:`\\underline{g} \\leq D \, u + C \, x \\leq \\bar{g}`
        at shooting nodes (0 to N-1).
        Type: :code:`np.ndarray`; default: :code:`np.array((0,0))`
        """
        return self.__D

    @property
    def lg(self):
        """:math:`\\underline{g}` - lower bound for general polytopic inequalities
        at shooting nodes (0 to N-1).
        Type: :code:`np.ndarray`; default: :code:`np.array([])`
        """
        return self.__lg

    @property
    def ug(self):
        """:math:`\\bar{g}` - upper bound for general polytopic inequalities
        at shooting nodes (0 to N-1).
        Type: :code:`np.ndarray`; default: :code:`np.array([])`.
        """
        return self.__ug

    # polytopic constraints at shooting node N
    @property
    def C_e(self):
        """:math:`C^e` - C matrix at terminal shooting node N.
        Type: :code:`np.ndarray`; default: :code:`np.array((0,0))`.
        """
        return self.__C_e

    @property
    def lg_e(self):
        """:math:`\\underline{g}^e` - lower bound on general polytopic inequalities
        at terminal shooting node N.
        Type: :code:`np.ndarray`; default: :code:`np.array([])`.
        """
        return self.__lg_e

    @property
    def ug_e(self):
        """:math:`\\bar{g}^e` - upper bound on general polytopic inequalities
        at terminal shooting node N.
        Type: :code:`np.ndarray`; default: :code:`np.array([])`.
        """
        return self.__ug_e


    # nonlinear constraints
    @property
    def lh(self):
        """:math:`\\underline{h}` - lower bound for nonlinear inequalities
        at intermediate shooting nodes (1 to N-1).
        Type: :code:`np.ndarray`; default: :code:`np.array([])`.
        """
        return self.__lh

    @property
    def uh(self):
        """:math:`\\bar{h}` - upper bound for nonlinear inequalities
        at intermediate shooting nodes (1 to N-1).
        Type: :code:`np.ndarray`; default: :code:`np.array([])`.
        """
        return self.__uh

    # nonlinear constraints at initial shooting node
    @property
    def lh_0(self):
        """:math:`\\underline{h}^0` - lower bound on nonlinear inequalities
        at initial shooting node (0).
        Type: :code:`np.ndarray`; default: :code:`np.array([])`.
        """
        return self.__lh_0

    @property
    def uh_0(self):
        """:math:`\\bar{h}^0` - upper bound on nonlinear inequalities
        at initial shooting node (0).
        Type: :code:`np.ndarray`; default: :code:`np.array([])`.
        """
        return self.__uh_0

    # nonlinear constraints at shooting node N
    @property
    def lh_e(self):
        """:math:`\\underline{h}^e` - lower bound on nonlinear inequalities
        at terminal shooting node N.
        Type: :code:`np.ndarray`; default: :code:`np.array([])`.
        """
        return self.__lh_e

    @property
    def uh_e(self):
        """:math:`\\bar{h}^e` - upper bound on nonlinear inequalities
        at terminal shooting node N.
        Type: :code:`np.ndarray`; default: :code:`np.array([])`.
        """
        return self.__uh_e

    # convex-over-nonlinear constraints
    @property
    def lphi(self):
        """:math:`\\underline{\phi}` - lower bound for convex-over-nonlinear inequalities
        at shooting nodes (0 to N-1).
        Type: :code:`np.ndarray`; default: :code:`np.array([])`.
        """
        return self.__lphi

    @property
    def uphi(self):
        """:math:`\\bar{\phi}` - upper bound for convex-over-nonlinear inequalities
        at shooting nodes (0 to N-1).
        Type: :code:`np.ndarray`; default: :code:`np.array([])`.
        """
        return self.__uphi

    # convex-over-nonlinear constraints at shooting node N
    @property
    def lphi_e(self):
        """:math:`\\underline{\phi}^e` - lower bound on convex-over-nonlinear inequalities
        at terminal shooting node N.
        Type: :code:`np.ndarray`; default: :code:`np.array([])`.
        """
        return self.__lphi_e

    @property
    def uphi_e(self):
        """:math:`\\bar{\phi}^e` - upper bound on convex-over-nonlinear inequalities
        at terminal shooting node N.
        Type: :code:`np.ndarray`; default: :code:`np.array([])`.
        """
        return self.__uphi_e

    @property
    def lphi_0(self):
        """:math:`\\underline{\phi}^0` - lower bound on convex-over-nonlinear inequalities
        at shooting node 0.
        Type: :code:`np.ndarray`; default: :code:`np.array([])`.
        """
        return self.__lphi_0

    @property
    def uphi_0(self):
        """:math:`\\bar{\phi}^0` - upper bound on convex-over-nonlinear inequalities
        at shooting node 0.
        Type: :code:`np.ndarray`; default: :code:`np.array([])`.
        """
        return self.__uphi_0


    # SLACK bounds
    # soft bounds on x
    @property
    def lsbx(self):
        """Lower bounds on slacks corresponding to soft lower bounds on x
        at stages (1 to N-1);
        not required - zeros by default"""
        return self.__lsbx

    @property
    def usbx(self):
        """Lower bounds on slacks corresponding to soft upper bounds on x
        at stages (1 to N-1);
        not required - zeros by default"""
        return self.__usbx

    @property
    def idxsbx(self):
        """Indices of soft bounds on x within the indices of bounds on x
        at stages (1 to N-1).
        Can be set by using :py:attr:`Jsbx`.
        Type: :code:`np.ndarray`; default: :code:`np.array([])`"""
        return self.__idxsbx

    @property
    def Jsbx(self):
        """:math:`J_{sbx}` - matrix coefficient for soft bounds on x
        at stages (1 to N-1);
        Translated internally into :py:attr:`idxsbx`."""
        print_J_to_idx_note()
        return self.__idxsbx

    # soft bounds on u
    @property
    def lsbu(self):
        """Lower bounds on slacks corresponding to soft lower bounds on u
        at stages (0 to N-1).
        Not required - zeros by default."""
        return self.__lsbu

    @property
    def usbu(self):
        """Lower bounds on slacks corresponding to soft upper bounds on u
        at stages (0 to N-1);
        not required - zeros by default"""
        return self.__usbu

    @property
    def idxsbu(self):
        """Indices of soft bounds on u within the indices of bounds on u
        at stages (0 to N-1).
        Can be set by using :py:attr:`Jsbu`.
        Type: :code:`np.ndarray`; default: :code:`np.array([])`"""
        return self.__idxsbu

    @property
    def Jsbu(self):
        """:math:`J_{sbu}` - matrix coefficient for soft bounds on u
        at stages (0 to N-1);
        internally translated into :py:attr:`idxsbu`"""
        print_J_to_idx_note()
        return self.__idxsbu

    # soft bounds on x at shooting node N
    @property
    def lsbx_e(self):
        """Lower bounds on slacks corresponding to soft lower bounds on x at shooting node N.
        Not required - zeros by default"""
        return self.__lsbx_e

    @property
    def usbx_e(self):
        """Lower bounds on slacks corresponding to soft upper bounds on x at shooting node N.
        Not required - zeros by default"""
        return self.__usbx_e

    @property
    def idxsbx_e(self):
        """Indices of soft bounds on x at shooting node N, within the indices of bounds on x at terminal shooting node N.
        Can be set by using :py:attr:`Jsbx_e`.
        Type: :code:`np.ndarray`; default: :code:`np.array([])`"""
        return self.__idxsbx_e

    @property
    def Jsbx_e(self):
        """:math:`J_{sbx}^e` - matrix coefficient for soft bounds on x at terminal shooting node N.
        Translated internally to :py:attr:`idxsbx_e`"""
        print_J_to_idx_note()
        return self.__idxsbx_e

    # soft general linear constraints
    @property
    def lsg(self):
        """Lower bounds on slacks corresponding to soft lower bounds for general linear constraints
        at stages (0 to N-1).
        Type: :code:`np.ndarray`; default: :code:`np.array([])`
        """
        return self.__lsg

    @property
    def usg(self):
        """Lower bounds on slacks corresponding to soft upper bounds for general linear constraints.
        Not required - zeros by default"""
        return self.__usg

    @property
    def idxsg(self):
        """Indices of soft general linear constraints within the indices of general linear constraints.
        Can be set by using :py:attr:`Jsg`.
        Type: :code:`np.ndarray`; default: :code:`np.array([])`"""
        return self.__idxsg

    @property
    def Jsg(self):
        """:math:`J_{sg}` - matrix coefficient for soft bounds on general linear constraints.
        Translated internally to :py:attr:`idxsg`"""
        print_J_to_idx_note()
        return self.__idxsg

    # soft nonlinear constraints
    @property
    def lsh(self):
        """Lower bounds on slacks corresponding to soft lower bounds for nonlinear constraints.
        Not required - zeros by default"""
        return self.__lsh

    @property
    def ush(self):
        """Lower bounds on slacks corresponding to soft upper bounds for nonlinear constraints.
        Not required - zeros by default"""
        return self.__ush

    @property
    def idxsh(self):
        """Indices of soft nonlinear constraints within the indices of nonlinear constraints.
        Can be set by using :py:attr:`Jbx`.
        Type: :code:`np.ndarray`; default: :code:`np.array([])`"""
        return self.__idxsh

    @property
    def Jsh(self):
        """:math:`J_{sh}` - matrix coefficient for soft bounds on nonlinear constraints.
        Translated internally to :py:attr:`idxsh`"""
        print_J_to_idx_note()
        return self.__idxsh

    # soft bounds on convex-over-nonlinear constraints
    @property
    def lsphi(self):
        """Lower bounds on slacks corresponding to soft lower bounds for convex-over-nonlinear constraints.
        Not required - zeros by default"""
        return self.__lsphi

    @property
    def usphi(self):
        """Lower bounds on slacks corresponding to soft upper bounds for convex-over-nonlinear constraints.
        Not required - zeros by default"""
        return self.__usphi

    @property
    def idxsphi(self):
        """Indices of soft convex-over-nonlinear constraints within the indices of nonlinear constraints.
        Can be set by using :py:attr:`Jsphi`.
        Type: :code:`np.ndarray`; default: :code:`np.array([])`"""
        return self.__idxsphi

    @property
    def Jsphi(self):
        """:math:`J_{s, \phi}` - matrix coefficient for soft bounds on convex-over-nonlinear constraints.
        Translated internally into :py:attr:`idxsphi`."""
        print_J_to_idx_note()
        return self.__idxsphi


    # soft bounds on general linear constraints at shooting node N
    @property
    def lsg_e(self):
        """Lower bounds on slacks corresponding to soft lower bounds for general linear constraints at shooting node N.
        Not required - zeros by default"""
        return self.__lsg_e

    @property
    def usg_e(self):
        """Lower bounds on slacks corresponding to soft upper bounds for general linear constraints at shooting node N.
        Not required - zeros by default"""
        return self.__usg_e

    @property
    def idxsg_e(self):
        """Indices of soft general linear constraints at shooting node N within the indices of general linear constraints at shooting node N.
        Can be set by using :py:attr:`Jsg_e`."""
        return self.__idxsg_e

    @property
    def Jsg_e(self):
        """:math:`J_{s,h}^e` - matrix coefficient for soft bounds on general linear constraints at terminal shooting node N.
        Translated internally to :py:attr:`idxsg_e`"""
        print_J_to_idx_note()
        return self.__idxsg_e


    # soft bounds on nonlinear constraints at shooting node 0
    @property
    def lsh_0(self):
        """Lower bounds on slacks corresponding to soft lower bounds for nonlinear constraints at initial shooting node 0.
        Not required - zeros by default"""
        return self.__lsh_0

    @property
    def ush_0(self):
        """Lower bounds on slacks corresponding to soft upper bounds for nonlinear constraints at initial shooting node 0.
        Not required - zeros by default"""
        return self.__ush_0

    @property
    def idxsh_0(self):
        """Indices of soft nonlinear constraints at shooting node N within the indices of nonlinear constraints at initial shooting node 0.
        Can be set by using :py:attr:`Jsh_0`."""
        return self.__idxsh_0

    @property
    def Jsh_0(self):
        """:math:`J_{s,h}^0` - matrix coefficient for soft bounds on nonlinear constraints at initial shooting node 0; fills :py:attr:`idxsh_0`"""
        print_J_to_idx_note()
        return self.__idxsh_0

    # soft bounds on convex-over-nonlinear constraints at shooting node N
    @property
    def lsphi_0(self):
        """Lower bounds on slacks corresponding to soft lower bounds for convex-over-nonlinear constraints at initial shooting node 0.
        Not required - zeros by default"""
        return self.__lsphi_0

    @property
    def usphi_0(self):
        """Lower bounds on slacks corresponding to soft upper bounds for convex-over-nonlinear constraints at initial shooting node 0.
        Not required - zeros by default"""
        return self.__usphi_0

    @property
    def idxsphi_0(self):
        """Indices of soft nonlinear constraints at shooting node N within the indices of nonlinear constraints at initial shooting node 0.
        Can be set by using :py:attr:`Jsphi_0`.
        Type: :code:`np.ndarray`; default: :code:`np.array([])`"""
        return self.__idxsphi_0

    @property
    def Jsphi_0(self):
        """:math:`J_{sh}^0` - matrix coefficient for soft bounds on convex-over-nonlinear constraints at shooting node N.
        Translated internally to :py:attr:`idxsphi_0`"""
        print_J_to_idx_note()
        return self.__idxsphi_0


    # soft bounds on nonlinear constraints at shooting node N
    @property
    def lsh_e(self):
        """Lower bounds on slacks corresponding to soft lower bounds for nonlinear constraints at terminal shooting node N.
        Not required - zeros by default"""
        return self.__lsh_e

    @property
    def ush_e(self):
        """Lower bounds on slacks corresponding to soft upper bounds for nonlinear constraints at terminal shooting node N.
        Not required - zeros by default"""
        return self.__ush_e

    @property
    def idxsh_e(self):
        """Indices of soft nonlinear constraints at shooting node N within the indices of nonlinear constraints at terminal shooting node N.
        Can be set by using :py:attr:`Jsh_e`."""
        return self.__idxsh_e

    @property
    def Jsh_e(self):
        """:math:`J_{s,h}^e` - matrix coefficient for soft bounds on nonlinear constraints at terminal shooting node N; fills :py:attr:`idxsh_e`"""
        print_J_to_idx_note()
        return self.__idxsh_e

    # soft bounds on convex-over-nonlinear constraints at shooting node N
    @property
    def lsphi_e(self):
        """Lower bounds on slacks corresponding to soft lower bounds for convex-over-nonlinear constraints at terminal shooting node N.
        Not required - zeros by default"""
        return self.__lsphi_e

    @property
    def usphi_e(self):
        """Lower bounds on slacks corresponding to soft upper bounds for convex-over-nonlinear constraints at terminal shooting node N.
        Not required - zeros by default"""
        return self.__usphi_e

    @property
    def idxsphi_e(self):
        """Indices of soft nonlinear constraints at shooting node N within the indices of nonlinear constraints at terminal shooting node N.
        Can be set by using :py:attr:`Jsphi_e`.
        Type: :code:`np.ndarray`; default: :code:`np.array([])`"""
        return self.__idxsphi_e

    @property
    def Jsphi_e(self):
        """:math:`J_{sh}^e` - matrix coefficient for soft bounds on convex-over-nonlinear constraints at shooting node N.
        Translated internally to :py:attr:`idxsphi_e`"""
        print_J_to_idx_note()
        return self.__idxsphi_e

    @property
    def x0(self):
        """:math:`x_0 \\in \mathbb{R}^{n_x}` - initial state --
        Translated internally to :py:attr:`idxbx_0`, :py:attr:`lbx_0`, :py:attr:`ubx_0`, :py:attr:`idxbxe_0` """
        print("x0 is converted to lbx_0, ubx_0, idxbx_0")
        print("idxbx_0: ", self.__idxbx_0)
        print("lbx_0: ", self.__lbx_0)
        print("ubx_0: ", self.__ubx_0)
        print("idxbxe_0: ", self.__idxbxe_0)
        return None

    @property
    def has_x0(self):
        return self.__has_x0

    # SETTERS
    @constr_type.setter
    def constr_type(self, constr_type):
        constr_types = ('BGH', 'BGP')
        if constr_type in constr_types:
            self.__constr_type = constr_type
        else:
            raise Exception('Invalid constr_type value. Possible values are:\n\n' \
                    + ',\n'.join(constr_types) + '.\n\nYou have: ' + constr_type + '.\n\n')

    @constr_type_0.setter
    def constr_type_0(self, constr_type_0):
        constr_types = ('BGH', 'BGP')
        if constr_type_0 in constr_types:
            self.__constr_type_0 = constr_type_0
        else:
            raise Exception('Invalid constr_type_0 value. Possible values are:\n\n' \
                    + ',\n'.join(constr_types) + '.\n\nYou have: ' + constr_type_0 + '.\n\n')

    @constr_type_e.setter
    def constr_type_e(self, constr_type_e):
        constr_types = ('BGH', 'BGP')
        if constr_type_e in constr_types:
            self.__constr_type_e = constr_type_e
        else:
            raise Exception('Invalid constr_type_e value. Possible values are:\n\n' \
                    + ',\n'.join(constr_types) + '.\n\nYou have: ' + constr_type_e + '.\n\n')

    # initial x
    @lbx_0.setter
    def lbx_0(self, lbx_0):
        if isinstance(lbx_0, np.ndarray):
            self.__lbx_0 = lbx_0
        else:
            raise Exception('Invalid lbx_0 value.')

    @ubx_0.setter
    def ubx_0(self, ubx_0):
        if isinstance(ubx_0, np.ndarray):
            self.__ubx_0 = ubx_0
        else:
            raise Exception('Invalid ubx_0 value.')

    @idxbx_0.setter
    def idxbx_0(self, idxbx_0):
        if isinstance(idxbx_0, np.ndarray):
            self.__idxbx_0 = idxbx_0
        else:
            raise Exception('Invalid idxbx_0 value.')

    @Jbx_0.setter
    def Jbx_0(self, Jbx_0):
        if isinstance(Jbx_0, np.ndarray):
            self.__idxbx_0 = J_to_idx(Jbx_0)
        else:
            raise Exception('Invalid Jbx_0 value.')

    @idxbxe_0.setter
    def idxbxe_0(self, idxbxe_0):
        if isinstance(idxbxe_0, np.ndarray):
            self.__idxbxe_0 = idxbxe_0
        else:
            raise Exception('Invalid idxbxe_0 value.')

    @x0.setter
    def x0(self, x0):
        if isinstance(x0, np.ndarray):
            self.__lbx_0 = x0
            self.__ubx_0 = x0
            self.__idxbx_0 = np.arange(x0.size)
            self.__idxbxe_0 = np.arange(x0.size)
            self.__has_x0 = True
        else:
            raise Exception('Invalid x0 value.')

    # bounds on x
    @lbx.setter
    def lbx(self, lbx):
        if isinstance(lbx, np.ndarray):
            self.__lbx = lbx
        else:
            raise Exception('Invalid lbx value.')

    @ubx.setter
    def ubx(self, ubx):
        if isinstance(ubx, np.ndarray):
            self.__ubx = ubx
        else:
            raise Exception('Invalid ubx value.')

    @idxbx.setter
    def idxbx(self, idxbx):
        if isinstance(idxbx, np.ndarray):
            self.__idxbx = idxbx
        else:
            raise Exception('Invalid idxbx value.')

    @Jbx.setter
    def Jbx(self, Jbx):
        if isinstance(Jbx, np.ndarray):
            self.__idxbx = J_to_idx(Jbx)
        else:
            raise Exception('Invalid Jbx value.')

    # bounds on u
    @lbu.setter
    def lbu(self, lbu):
        if isinstance(lbu, np.ndarray):
            self.__lbu = lbu
        else:
            raise Exception('Invalid lbu value.')

    @ubu.setter
    def ubu(self, ubu):
        if isinstance(ubu, np.ndarray):
            self.__ubu = ubu
        else:
            raise Exception('Invalid ubu value.')

    @idxbu.setter
    def idxbu(self, idxbu):
        if isinstance(idxbu, np.ndarray):
            self.__idxbu = idxbu
        else:
            raise Exception('Invalid idxbu value.')

    @Jbu.setter
    def Jbu(self, Jbu):
        if isinstance(Jbu, np.ndarray):
            self.__idxbu = J_to_idx(Jbu)
        else:
            raise Exception('Invalid Jbu value.')

    # bounds on x at shooting node N
    @lbx_e.setter
    def lbx_e(self, lbx_e):
        if isinstance(lbx_e, np.ndarray):
            self.__lbx_e = lbx_e
        else:
            raise Exception('Invalid lbx_e value.')

    @ubx_e.setter
    def ubx_e(self, ubx_e):
        if isinstance(ubx_e, np.ndarray):
            self.__ubx_e = ubx_e
        else:
            raise Exception('Invalid ubx_e value.')

    @idxbx_e.setter
    def idxbx_e(self, idxbx_e):
        if isinstance(idxbx_e, np.ndarray):
            self.__idxbx_e = idxbx_e
        else:
            raise Exception('Invalid idxbx_e value.')

    @Jbx_e.setter
    def Jbx_e(self, Jbx_e):
        if isinstance(Jbx_e, np.ndarray):
            self.__idxbx_e = J_to_idx(Jbx_e)
        else:
            raise Exception('Invalid Jbx_e value.')

    # polytopic constraints
    @D.setter
    def D(self, D):
        if isinstance(D, np.ndarray) and len(D.shape) == 2:
            self.__D = D
        else:
            raise Exception('Invalid constraint D value.' \
                + 'Should be 2 dimensional numpy array.')

    @C.setter
    def C(self, C):
        if isinstance(C, np.ndarray) and len(C.shape) == 2:
            self.__C = C
        else:
            raise Exception('Invalid constraint C value.' \
                + 'Should be 2 dimensional numpy array.')

    @lg.setter
    def lg(self, lg):
        if isinstance(lg, np.ndarray):
            self.__lg = lg
        else:
            raise Exception('Invalid lg value.')

    @ug.setter
    def ug(self, ug):
        if isinstance(ug, np.ndarray):
            self.__ug = ug
        else:
            raise Exception('Invalid ug value.')

    # polytopic constraints at shooting node N
    @C_e.setter
    def C_e(self, C_e):
        if isinstance(C_e, np.ndarray) and len(C_e.shape) == 2:
            self.__C_e = C_e
        else:
            raise Exception('Invalid constraint C_e value.' \
                + 'Should be 2 dimensional numpy array.')

    @lg_e.setter
    def lg_e(self, lg_e):
        if isinstance(lg_e, np.ndarray):
            self.__lg_e = lg_e
        else:
            raise Exception('Invalid lg_e value.')

    @ug_e.setter
    def ug_e(self, ug_e):
        if isinstance(ug_e, np.ndarray):
            self.__ug_e = ug_e
        else:
            raise Exception('Invalid ug_e value.')

    # nonlinear constraints
    @lh.setter
    def lh(self, lh):
        if isinstance(lh, np.ndarray):
            self.__lh = lh
        else:
            raise Exception('Invalid lh value.')

    @uh.setter
    def uh(self, uh):
        if isinstance(uh, np.ndarray):
            self.__uh = uh
        else:
            raise Exception('Invalid uh value.')


    # nonlinear constraints at initial shooting node
    @lh_0.setter
    def lh_0(self, lh_0):
        if isinstance(lh_0, np.ndarray):
            self.__lh_0 = lh_0
        else:
            raise Exception('Invalid lh_0 value.')

    @uh_0.setter
    def uh_0(self, uh_0):
        if isinstance(uh_0, np.ndarray):
            self.__uh_0 = uh_0
        else:
            raise Exception('Invalid uh_0 value.')

    @lphi_0.setter
    def lphi_0(self, lphi_0):
        if isinstance(lphi_0, np.ndarray):
            self.__lphi_0 = lphi_0
        else:
            raise Exception('Invalid lphi_0 value.')

    @uphi_0.setter
    def uphi_0(self, uphi_0):
        if isinstance(uphi_0, np.ndarray):
            self.__uphi_0 = uphi_0
        else:
            raise Exception('Invalid uphi_0 value.')

    # nonlinear constraints at shooting node N
    @lh_e.setter
    def lh_e(self, lh_e):
        if isinstance(lh_e, np.ndarray):
            self.__lh_e = lh_e
        else:
            raise Exception('Invalid lh_e value.')

    @uh_e.setter
    def uh_e(self, uh_e):
        if isinstance(uh_e, np.ndarray):
            self.__uh_e = uh_e
        else:
            raise Exception('Invalid uh_e value.')

    # convex-over-nonlinear constraints
    @lphi.setter
    def lphi(self, lphi):
        if isinstance(lphi, np.ndarray):
            self.__lphi = lphi
        else:
            raise Exception('Invalid lphi value.')

    @uphi.setter
    def uphi(self, uphi):
        if isinstance(uphi, np.ndarray):
            self.__uphi = uphi
        else:
            raise Exception('Invalid uphi value.')

    # convex-over-nonlinear constraints at shooting node N
    @lphi_e.setter
    def lphi_e(self, lphi_e):
        if isinstance(lphi_e, np.ndarray):
            self.__lphi_e = lphi_e
        else:
            raise Exception('Invalid lphi_e value.')

    @uphi_e.setter
    def uphi_e(self, uphi_e):
        if isinstance(uphi_e, np.ndarray):
            self.__uphi_e = uphi_e
        else:
            raise Exception('Invalid uphi_e value.')

    # SLACK bounds
    # soft bounds on x
    @lsbx.setter
    def lsbx(self, lsbx):
        if isinstance(lsbx, np.ndarray):
            self.__lsbx = lsbx
        else:
            raise Exception('Invalid lsbx value.')

    @usbx.setter
    def usbx(self, usbx):
        if isinstance(usbx, np.ndarray):
            self.__usbx = usbx
        else:
            raise Exception('Invalid usbx value.')

    @idxsbx.setter
    def idxsbx(self, idxsbx):
        if isinstance(idxsbx, np.ndarray):
            self.__idxsbx = idxsbx
        else:
            raise Exception('Invalid idxsbx value.')

    @Jsbx.setter
    def Jsbx(self, Jsbx):
        if isinstance(Jsbx, np.ndarray):
            self.__idxsbx = J_to_idx_slack(Jsbx)
        else:
            raise Exception('Invalid Jsbx value, expected numpy array.')

    # soft bounds on u
    @lsbu.setter
    def lsbu(self, lsbu):
        if isinstance(lsbu, np.ndarray):
            self.__lsbu = lsbu
        else:
            raise Exception('Invalid lsbu value.')

    @usbu.setter
    def usbu(self, usbu):
        if isinstance(usbu, np.ndarray):
            self.__usbu = usbu
        else:
            raise Exception('Invalid usbu value.')

    @idxsbu.setter
    def idxsbu(self, idxsbu):
        if isinstance(idxsbu, np.ndarray):
            self.__idxsbu = idxsbu
        else:
            raise Exception('Invalid idxsbu value.')

    @Jsbu.setter
    def Jsbu(self, Jsbu):
        if isinstance(Jsbu, np.ndarray):
            self.__idxsbu = J_to_idx_slack(Jsbu)
        else:
            raise Exception('Invalid Jsbu value.')

    # soft bounds on x at shooting node N
    @lsbx_e.setter
    def lsbx_e(self, lsbx_e):
        if isinstance(lsbx_e, np.ndarray):
            self.__lsbx_e = lsbx_e
        else:
            raise Exception('Invalid lsbx_e value.')

    @usbx_e.setter
    def usbx_e(self, usbx_e):
        if isinstance(usbx_e, np.ndarray):
            self.__usbx_e = usbx_e
        else:
            raise Exception('Invalid usbx_e value.')

    @idxsbx_e.setter
    def idxsbx_e(self, idxsbx_e):
        if isinstance(idxsbx_e, np.ndarray):
            self.__idxsbx_e = idxsbx_e
        else:
            raise Exception('Invalid idxsbx_e value.')

    @Jsbx_e.setter
    def Jsbx_e(self, Jsbx_e):
        if isinstance(Jsbx_e, np.ndarray):
            self.__idxsbx_e = J_to_idx_slack(Jsbx_e)
        else:
            raise Exception('Invalid Jsbx_e value.')


    # soft bounds on general linear constraints
    @lsg.setter
    def lsg(self, lsg):
        if isinstance(lsg, np.ndarray):
            self.__lsg = lsg
        else:
            raise Exception('Invalid lsg value.')

    @usg.setter
    def usg(self, usg):
        if isinstance(usg, np.ndarray):
            self.__usg = usg
        else:
            raise Exception('Invalid usg value.')

    @idxsg.setter
    def idxsg(self, idxsg):
        if isinstance(idxsg, np.ndarray):
            self.__idxsg = idxsg
        else:
            raise Exception('Invalid idxsg value.')

    @Jsg.setter
    def Jsg(self, Jsg):
        if isinstance(Jsg, np.ndarray):
            self.__idxsg = J_to_idx_slack(Jsg)
        else:
            raise Exception('Invalid Jsg value, expected numpy array.')


    # soft bounds on nonlinear constraints
    @lsh.setter
    def lsh(self, lsh):
        if isinstance(lsh, np.ndarray):
            self.__lsh = lsh
        else:
            raise Exception('Invalid lsh value.')

    @ush.setter
    def ush(self, ush):
        if isinstance(ush, np.ndarray):
            self.__ush = ush
        else:
            raise Exception('Invalid ush value.')

    @idxsh.setter
    def idxsh(self, idxsh):
        if isinstance(idxsh, np.ndarray):
            self.__idxsh = idxsh
        else:
            raise Exception('Invalid idxsh value.')


    @Jsh.setter
    def Jsh(self, Jsh):
        if isinstance(Jsh, np.ndarray):
            self.__idxsh = J_to_idx_slack(Jsh)
        else:
            raise Exception('Invalid Jsh value, expected numpy array.')

    # soft bounds on convex-over-nonlinear constraints
    @lsphi.setter
    def lsphi(self, lsphi):
        if isinstance(lsphi, np.ndarray):
            self.__lsphi = lsphi
        else:
            raise Exception('Invalid lsphi value.')

    @usphi.setter
    def usphi(self, usphi):
        if isinstance(usphi, np.ndarray):
            self.__usphi = usphi
        else:
            raise Exception('Invalid usphi value.')

    @idxsphi.setter
    def idxsphi(self, idxsphi):
        if isinstance(idxsphi, np.ndarray):
            self.__idxsphi = idxsphi
        else:
            raise Exception('Invalid idxsphi value.')

    @Jsphi.setter
    def Jsphi(self, Jsphi):
        if isinstance(Jsphi, np.ndarray):
            self.__idxsphi = J_to_idx_slack(Jsphi)
        else:
            raise Exception('Invalid Jsphi value, expected numpy array.')

    # soft bounds on general linear constraints at shooting node N
    @lsg_e.setter
    def lsg_e(self, lsg_e):
        if isinstance(lsg_e, np.ndarray):
            self.__lsg_e = lsg_e
        else:
            raise Exception('Invalid lsg_e value.')

    @usg_e.setter
    def usg_e(self, usg_e):
        if isinstance(usg_e, np.ndarray):
            self.__usg_e = usg_e
        else:
            raise Exception('Invalid usg_e value.')

    @idxsg_e.setter
    def idxsg_e(self, idxsg_e):
        if isinstance(idxsg_e, np.ndarray):
            self.__idxsg_e = idxsg_e
        else:
            raise Exception('Invalid idxsg_e value.')

    @Jsg_e.setter
    def Jsg_e(self, Jsg_e):
        if isinstance(Jsg_e, np.ndarray):
            self.__idxsg_e = J_to_idx_slack(Jsg_e)
        else:
            raise Exception('Invalid Jsg_e value, expected numpy array.')

    # soft bounds on nonlinear constraints at shooting node N
    @lsh_e.setter
    def lsh_e(self, lsh_e):
        if isinstance(lsh_e, np.ndarray):
            self.__lsh_e = lsh_e
        else:
            raise Exception('Invalid lsh_e value.')

    @ush_e.setter
    def ush_e(self, ush_e):
        if isinstance(ush_e, np.ndarray):
            self.__ush_e = ush_e
        else:
            raise Exception('Invalid ush_e value.')

    @idxsh_e.setter
    def idxsh_e(self, idxsh_e):
        if isinstance(idxsh_e, np.ndarray):
            self.__idxsh_e = idxsh_e
        else:
            raise Exception('Invalid idxsh_e value.')

    @Jsh_e.setter
    def Jsh_e(self, Jsh_e):
        if isinstance(Jsh_e, np.ndarray):
            self.__idxsh_e = J_to_idx_slack(Jsh_e)
        else:
            raise Exception('Invalid Jsh_e value, expected numpy array.')


    # soft bounds on convex-over-nonlinear constraints at shooting node N
    @lsphi_e.setter
    def lsphi_e(self, lsphi_e):
        if isinstance(lsphi_e, np.ndarray):
            self.__lsphi_e = lsphi_e
        else:
            raise Exception('Invalid lsphi_e value.')

    @usphi_e.setter
    def usphi_e(self, usphi_e):
        if isinstance(usphi_e, np.ndarray):
            self.__usphi_e = usphi_e
        else:
            raise Exception('Invalid usphi_e value.')

    @idxsphi_e.setter
    def idxsphi_e(self, idxsphi_e):
        if isinstance(idxsphi_e, np.ndarray):
            self.__idxsphi_e = idxsphi_e
        else:
            raise Exception('Invalid idxsphi_e value.')

    @Jsphi_e.setter
    def Jsphi_e(self, Jsphi_e):
        if isinstance(Jsphi_e, np.ndarray):
            self.__idxsphi_e = J_to_idx_slack(Jsphi_e)
        else:
            raise Exception('Invalid Jsphi_e value.')

    # soft constraints at shooting node 0
    @lsh_0.setter
    def lsh_0(self, lsh_0):
        if isinstance(lsh_0, np.ndarray):
            self.__lsh_0 = lsh_0
        else:
            raise Exception('Invalid lsh_0 value.')

    @ush_0.setter
    def ush_0(self, ush_0):
        if isinstance(ush_0, np.ndarray):
            self.__ush_0 = ush_0
        else:
            raise Exception('Invalid ush_0 value.')

    @idxsh_0.setter
    def idxsh_0(self, idxsh_0):
        if isinstance(idxsh_0, np.ndarray):
            self.__idxsh_0 = idxsh_0
        else:
            raise Exception('Invalid idxsh_0 value.')

    @Jsh_0.setter
    def Jsh_0(self, Jsh_0):
        if isinstance(Jsh_0, np.ndarray):
            self.__idxsh_0 = J_to_idx_slack(Jsh_0)
        else:
            raise Exception('Invalid Jsh_0 value, expected numpy array.')

    @lsphi_0.setter
    def lsphi_0(self, lsphi_0):
        if isinstance(lsphi_0, np.ndarray):
            self.__lsphi_0 = lsphi_0
        else:
            raise Exception('Invalid lsphi_0 value.')

    @usphi_0.setter
    def usphi_0(self, usphi_0):
        if isinstance(usphi_0, np.ndarray):
            self.__usphi_0 = usphi_0
        else:
            raise Exception('Invalid usphi_0 value.')

    @idxsphi_0.setter
    def idxsphi_0(self, idxsphi_0):
        if isinstance(idxsphi_0, np.ndarray):
            self.__idxsphi_0 = idxsphi_0
        else:
            raise Exception('Invalid idxsphi_0 value.')

    @Jsphi_0.setter
    def Jsphi_0(self, Jsphi_0):
        if isinstance(Jsphi_0, np.ndarray):
            self.__idxsphi_0 = J_to_idx_slack(Jsphi_0)
        else:
            raise Exception('Invalid Jsphi_0 value.')

    def set(self, attr, value):
        setattr(self, attr, value)
