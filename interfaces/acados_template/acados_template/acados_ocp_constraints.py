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
from .utils import J_to_idx, print_J_to_idx_note, J_to_idx_slack, cast_to_1d_nparray, cast_to_2d_nparray, is_empty

class AcadosOcpConstraints:
    """
    Class containing the description of the constraints.

    Soft constraints can be formulated in two ways:
    1) via idxsbu, idxsbx, idxsg, idxsh, idxsphi, lsbu, usbu, lsbx, usbx, lsg, usg, lsh, ush, lsphi, usphi
    2) via idxs_rev, ls, us and *_0, *_e variants

    Option 1) is what was implemented in acados <= 0.5.1
    Option 2) is the new way of formulating soft constraints, which is more flexible, as one slack variable can be used for multiple constraints.
    """
    def __init__(self):
        self.__constr_types = ('BGH', 'BGP')

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
        self.__lg_e    = np.array([])
        self.__ug_e    = np.array([])
        self.__C_e     = np.zeros((0,0))
        # nonlinear constraints at initial shooting node
        self.__lh_0    = np.array([])
        self.__uh_0    = np.array([])
        # nonlinear constraints
        self.__lh      = np.array([])
        self.__uh      = np.array([])
        # nonlinear constraints at shooting node N
        self.__lh_e    = np.array([])
        self.__uh_e    = np.array([])
        # convex-over-nonlinear constraints
        self.__lphi_0 = np.array([])
        self.__uphi_0 = np.array([])
        self.__lphi    = np.array([])
        self.__uphi    = np.array([])
        self.__lphi_e = np.array([])
        self.__uphi_e = np.array([])

        # idxs_rev slack formulation
        self.__idxs_rev_0 = np.array([])
        self.__idxs_rev = np.array([])
        self.__idxs_rev_e = np.array([])

        self.__ls_0 = np.array([])
        self.__ls = np.array([])
        self.__ls_e = np.array([])

        self.__us_0 = np.array([])
        self.__us = np.array([])
        self.__us_e = np.array([])

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
        self.__idxsbx_e = np.array([])
        # soft bounds on general linear constraints
        self.__lsg    = np.array([])
        self.__usg    = np.array([])
        self.__idxsg  = np.array([])
        # soft bounds on general linear constraints at shooting node N
        self.__lsg_e    = np.array([])
        self.__usg_e    = np.array([])
        self.__idxsg_e  = np.array([])
        # soft bounds on nonlinear constraints at shooting node 0
        self.__lsh_0    = np.array([])
        self.__ush_0    = np.array([])
        self.__idxsh_0  = np.array([])
        # soft bounds on nonlinear constraints
        self.__lsh    = np.array([])
        self.__ush    = np.array([])
        self.__idxsh  = np.array([])
        # soft bounds on nonlinear constraints at shooting node N
        self.__lsh_e    = np.array([])
        self.__ush_e    = np.array([])
        self.__idxsh_e  = np.array([])
        # soft bounds on convex-over-nonlinear constraint (BGP) at shooting node 0
        self.__lsphi_0    = np.array([])
        self.__usphi_0    = np.array([])
        self.__idxsphi_0  = np.array([])
        # soft bounds on convex-over-nonlinear constraint (BGP)
        self.__lsphi  = np.array([])
        self.__usphi  = np.array([])
        self.__idxsphi  = np.array([])
        # soft bounds on convex-over-nonlinear constraint (BGP) at shooting node N
        self.__lsphi_e    = np.array([])
        self.__usphi_e    = np.array([])
        self.__idxsphi_e  = np.array([])


    # types
    @property
    def constr_type(self):
        """Constraints type for shooting nodes (0 to N-1). string in {BGH, BGP}.
        Detected automatically from model.
        Default: BGH; BGP is for convex over nonlinear.
        """
        return self.__constr_type

    @constr_type.setter
    def constr_type(self, constr_type):
        if constr_type in self.__constr_types:
            self.__constr_type = constr_type
        else:
            raise ValueError('Invalid constr_type value. Possible values are:\n\n' \
                    + ',\n'.join(self.__constr_types) + '.\n\nYou have: ' + constr_type + '.\n\n')

    @property
    def constr_type_0(self):
        """Constraints type for initial shooting node. string in {BGH, BGP}.
        Detected automatically from model.
        Default: BGH; BGP is for convex over nonlinear.
        """
        return self.__constr_type_0

    @constr_type_0.setter
    def constr_type_0(self, constr_type_0):
        if constr_type_0 in self.__constr_types:
            self.__constr_type_0 = constr_type_0
        else:
            raise ValueError('Invalid constr_type_0 value. Possible values are:\n\n' \
                    + ',\n'.join(self.__constr_types) + '.\n\nYou have: ' + constr_type_0 + '.\n\n')

    @property
    def constr_type_e(self):
        """Constraints type for terminal shooting node N. string in {BGH, BGP}.
        Detected automatically from model.
        Default: BGH; BGP is for convex over nonlinear.
        """
        return self.__constr_type_e

    @constr_type_e.setter
    def constr_type_e(self, constr_type_e):
        if constr_type_e in self.__constr_types:
            self.__constr_type_e = constr_type_e
        else:
            raise ValueError('Invalid constr_type_e value. Possible values are:\n\n' \
                    + ',\n'.join(self.__constr_types) + '.\n\nYou have: ' + constr_type_e + '.\n\n')

    @property
    def lbx_0(self):
        r""":math:`\underline{x_0}` - lower bounds on x at initial stage 0.
        Type: :code:`np.ndarray`; default: :code:`np.array([])`."""
        return self.__lbx_0

    @lbx_0.setter
    def lbx_0(self, lbx_0):
        self.__lbx_0 = cast_to_1d_nparray(lbx_0, "lbx_0")

    @property
    def ubx_0(self):
        r""":math:`\bar{x_0}` - upper bounds on x at initial stage 0.
        Type: :code:`np.ndarray`; default: :code:`np.array([])`"""
        return self.__ubx_0

    @ubx_0.setter
    def ubx_0(self, ubx_0):
        self.__ubx_0 = cast_to_1d_nparray(ubx_0, "ubx_0")

    @property
    def Jbx_0(self):
        """:math:`J_{bx,0}` - matrix coefficient for bounds on x at initial stage 0.
        Translated internally to :py:attr:`idxbx_0`"""
        print_J_to_idx_note()
        return self.__idxbx_0

    @Jbx_0.setter
    def Jbx_0(self, Jbx_0):
        Jbx_0 = cast_to_2d_nparray(Jbx_0, "Jbx_0")
        self.__idxbx_0 = J_to_idx(Jbx_0)

    @property
    def idxbx_0(self):
        """Indices of bounds on x at initial stage 0
        -- can be set automatically via x0.
        Can be set by using :py:attr:`Jbx_0`.
        Type: :code:`np.ndarray`; default: :code:`np.array([])`"""
        return self.__idxbx_0

    @idxbx_0.setter
    def idxbx_0(self, idxbx_0):
        self.__idxbx_0 = cast_to_1d_nparray(idxbx_0, "idxbx_0")

    @property
    def idxbxe_0(self):
        """Indices of bounds on x0 that are equalities -- can be set automatically via :py:attr:`x0`.
        Type: :code:`np.ndarray`; default: :code:`np.array([])`"""
        return self.__idxbxe_0

    @idxbxe_0.setter
    def idxbxe_0(self, idxbxe_0):
        self.__idxbxe_0 = cast_to_1d_nparray(idxbxe_0, "idxbxe_0")

    @property
    def lbx(self):
        r""":math:`\underline{x}` - lower bounds on x at intermediate shooting nodes (1 to N-1).
        Type: :code:`np.ndarray`; default: :code:`np.array([])`"""
        return self.__lbx

    @lbx.setter
    def lbx(self, lbx):
        self.__lbx = cast_to_1d_nparray(lbx, "lbx")

    @property
    def ubx(self):
        r""":math:`\bar{x}` - upper bounds on x at intermediate shooting nodes (1 to N-1).
        Type: :code:`np.ndarray`; default: :code:`np.array([])`"""
        return self.__ubx

    @ubx.setter
    def ubx(self, ubx):
        self.__ubx = cast_to_1d_nparray(ubx, "ubx")

    @property
    def idxbx(self):
        """indices of bounds on x (defines :math:`J_{bx}`) at intermediate shooting nodes (1 to N-1).
        Can be set by using :py:attr:`Jbx`.
        Type: :code:`np.ndarray`; default: :code:`np.array([])`"""
        return self.__idxbx

    @idxbx.setter
    def idxbx(self, idxbx):
        self.__idxbx = cast_to_1d_nparray(idxbx, "idxbx")

    @property
    def Jbx(self):
        """:math:`J_{bx}` - matrix coefficient for bounds on x
        at intermediate shooting nodes (1 to N-1).
        Translated internally into :py:attr:`idxbx`."""
        print_J_to_idx_note()
        return self.__idxbx

    @Jbx.setter
    def Jbx(self, Jbx):
        Jbx = cast_to_2d_nparray(Jbx, "Jbx")
        self.__idxbx = J_to_idx(Jbx)

    @property
    def lbx_e(self):
        r""":math:`\underline{x}^e` - lower bounds on x at terminal shooting node N.
        Type: :code:`np.ndarray`; default: :code:`np.array([])`"""
        return self.__lbx_e

    @lbx_e.setter
    def lbx_e(self, lbx_e):
        self.__lbx_e = cast_to_1d_nparray(lbx_e, "lbx_e")

    @property
    def ubx_e(self):
        r""":math:`\bar{x}^e` - upper bounds on x at terminal shooting node N.
        Type: :code:`np.ndarray`; default: :code:`np.array([])`"""
        return self.__ubx_e

    @ubx_e.setter
    def ubx_e(self, ubx_e):
        self.__ubx_e = cast_to_1d_nparray(ubx_e, "ubx_e")

    @property
    def idxbx_e(self):
        """Indices for bounds on x at terminal shooting node N (defines :math:`J_{bx}^e`).
        Can be set by using :py:attr:`Jbx_e`.
        Type: :code:`np.ndarray`; default: :code:`np.array([])`"""
        return self.__idxbx_e

    @idxbx_e.setter
    def idxbx_e(self, idxbx_e):
        self.__idxbx_e = cast_to_1d_nparray(idxbx_e, "idxbx_e")

    @property
    def Jbx_e(self):
        """:math:`J_{bx}^e` matrix coefficient for bounds on x at terminal shooting node N.
        Translated internally into :py:attr:`idxbx_e`."""
        print_J_to_idx_note()
        return self.__idxbx_e

    @Jbx_e.setter
    def Jbx_e(self, Jbx_e):
        Jbx_e = cast_to_2d_nparray(Jbx_e, "Jbx_e")
        self.__idxbx_e = J_to_idx(Jbx_e)

    @property
    def lbu(self):
        r""":math:`\underline{u}` - lower bounds on u at shooting nodes (0 to N-1).
        Type: :code:`np.ndarray`; default: :code:`np.array([])`
        """
        return self.__lbu

    @lbu.setter
    def lbu(self, lbu):
        self.__lbu = cast_to_1d_nparray(lbu, "lbu")

    @property
    def ubu(self):
        r""":math:`\bar{u}` - upper bounds on u at shooting nodes (0 to N-1).
        Type: :code:`np.ndarray`; default: :code:`np.array([])`
        """
        return self.__ubu

    @ubu.setter
    def ubu(self, ubu):
        self.__ubu = cast_to_1d_nparray(ubu, "ubu")

    @property
    def idxbu(self):
        """Indices of bounds on u (defines :math:`J_{bu}`) at shooting nodes (0 to N-1).
        Can be set by using :py:attr:`Jbu`.
        Type: :code:`np.ndarray`; default: :code:`np.array([])`
        """
        return self.__idxbu

    @idxbu.setter
    def idxbu(self, idxbu):
        self.__idxbu = cast_to_1d_nparray(idxbu, "idxbu")

    @property
    def Jbu(self):
        """:math:`J_{bu}` - matrix coefficient for bounds on u at shooting nodes (0 to N-1).
        Translated internally to :py:attr:`idxbu`.
        """
        print_J_to_idx_note()
        return self.__idxbu

    @Jbu.setter
    def Jbu(self, Jbu):
        Jbu = cast_to_2d_nparray(Jbu, "Jbu")
        self.__idxbu = J_to_idx(Jbu)

    @property
    def C(self):
        r""":math:`C` - C matrix in :math:`\underline{g} \leq D \, u + C \, x \leq \bar{g}`
        at shooting nodes (0 to N-1).
        Type: :code:`np.ndarray`; default: :code:`np.array((0,0))`.
        """
        return self.__C

    @C.setter
    def C(self, C):
        self.__C = cast_to_2d_nparray(C, "C")

    @property
    def D(self):
        r""":math:`D` - D matrix in :math:`\underline{g} \leq D \, u + C \, x \leq \bar{g}`
        at shooting nodes (0 to N-1).
        Type: :code:`np.ndarray`; default: :code:`np.array((0,0))`
        """
        return self.__D

    @D.setter
    def D(self, D):
        self.__D = cast_to_2d_nparray(D, "D")

    @property
    def lg(self):
        r""":math:`\underline{g}` - lower bound for general polytopic inequalities
        at shooting nodes (0 to N-1).
        Type: :code:`np.ndarray`; default: :code:`np.array([])`
        """
        return self.__lg

    @lg.setter
    def lg(self, value):
        self.__lg = cast_to_1d_nparray(value, 'lg')

    @property
    def ug(self):
        r""":math:`\bar{g}` - upper bound for general polytopic inequalities
        at shooting nodes (0 to N-1).
        Type: :code:`np.ndarray`; default: :code:`np.array([])`.
        """
        return self.__ug

    @ug.setter
    def ug(self, value):
        self.__ug = cast_to_1d_nparray(value, 'ug')

    @property
    def C_e(self):
        """:math:`C^e` - C matrix at terminal shooting node N.
        Type: :code:`np.ndarray`; default: :code:`np.array((0,0))`.
        """
        return self.__C_e

    @C_e.setter
    def C_e(self, C_e):
        self.__C_e = cast_to_2d_nparray(C_e, "C_e")

    @property
    def lg_e(self):
        r""":math:`\underline{g}^e` - lower bound on general polytopic inequalities
        at terminal shooting node N.
        Type: :code:`np.ndarray`; default: :code:`np.array([])`.
        """
        return self.__lg_e

    @lg_e.setter
    def lg_e(self, value):
        self.__lg_e = cast_to_1d_nparray(value, 'lg_e')

    @property
    def ug_e(self):
        r""":math:`\bar{g}^e` - upper bound on general polytopic inequalities
        at terminal shooting node N.
        Type: :code:`np.ndarray`; default: :code:`np.array([])`.
        """
        return self.__ug_e


    @ug_e.setter
    def ug_e(self, value):
        self.__ug_e = cast_to_1d_nparray(value, 'ug_e')

    @property
    def lh(self):
        r""":math:`\underline{h}` - lower bound for nonlinear inequalities
        at intermediate shooting nodes (1 to N-1).
        Type: :code:`np.ndarray`; default: :code:`np.array([])`.
        """
        return self.__lh

    @lh.setter
    def lh(self, value):
        self.__lh = cast_to_1d_nparray(value, 'lh')

    @property
    def uh(self):
        r""":math:`\bar{h}` - upper bound for nonlinear inequalities
        at intermediate shooting nodes (1 to N-1).
        Type: :code:`np.ndarray`; default: :code:`np.array([])`.
        """
        return self.__uh

    @uh.setter
    def uh(self, value):
        self.__uh = cast_to_1d_nparray(value, 'uh')

    @property
    def lh_0(self):
        r""":math:`\underline{h}^0` - lower bound on nonlinear inequalities
        at initial shooting node (0).
        Type: :code:`np.ndarray`; default: :code:`np.array([])`.
        """
        return self.__lh_0

    @lh_0.setter
    def lh_0(self, value):
        self.__lh_0 = cast_to_1d_nparray(value, 'lh_0')

    @property
    def uh_0(self):
        r""":math:`\bar{h}^0` - upper bound on nonlinear inequalities
        at initial shooting node (0).
        Type: :code:`np.ndarray`; default: :code:`np.array([])`.
        """
        return self.__uh_0

    @uh_0.setter
    def uh_0(self, value):
        self.__uh_0 = cast_to_1d_nparray(value, 'uh_0')

    @property
    def lh_e(self):
        r""":math:`\underline{h}^e` - lower bound on nonlinear inequalities
        at terminal shooting node N.
        Type: :code:`np.ndarray`; default: :code:`np.array([])`.
        """
        return self.__lh_e

    @lh_e.setter
    def lh_e(self, value):
        self.__lh_e = cast_to_1d_nparray(value, 'lh_e')

    @property
    def uh_e(self):
        r""":math:`\bar{h}^e` - upper bound on nonlinear inequalities
        at terminal shooting node N.
        Type: :code:`np.ndarray`; default: :code:`np.array([])`.
        """
        return self.__uh_e

    @uh_e.setter
    def uh_e(self, value):
        self.__uh_e = cast_to_1d_nparray(value, 'uh_e')

    @property
    def lphi(self):
        r""":math:`\underline{\phi}` - lower bound for convex-over-nonlinear inequalities
        at shooting nodes (0 to N-1).
        Type: :code:`np.ndarray`; default: :code:`np.array([])`.
        """
        return self.__lphi

    @lphi.setter
    def lphi(self, value):
        self.__lphi = cast_to_1d_nparray(value, 'lphi')

    @property
    def uphi(self):
        r""":math:`\bar{\phi}` - upper bound for convex-over-nonlinear inequalities
        at shooting nodes (0 to N-1).
        Type: :code:`np.ndarray`; default: :code:`np.array([])`.
        """
        return self.__uphi

    @uphi.setter
    def uphi(self, value):
        self.__uphi = cast_to_1d_nparray(value, 'uphi')

    @property
    def lphi_e(self):
        r""":math:`\underline{\phi}^e` - lower bound on convex-over-nonlinear inequalities
        at terminal shooting node N.
        Type: :code:`np.ndarray`; default: :code:`np.array([])`.
        """
        return self.__lphi_e

    @lphi_e.setter
    def lphi_e(self, value):
        self.__lphi_e = cast_to_1d_nparray(value, 'lphi_e')

    @property
    def uphi_e(self):
        r""":math:`\bar{\phi}^e` - upper bound on convex-over-nonlinear inequalities
        at terminal shooting node N.
        Type: :code:`np.ndarray`; default: :code:`np.array([])`.
        """
        return self.__uphi_e

    @uphi_e.setter
    def uphi_e(self, value):
        self.__uphi_e = cast_to_1d_nparray(value, 'uphi_e')

    @property
    def lphi_0(self):
        r""":math:`\underline{\phi}^0` - lower bound on convex-over-nonlinear inequalities
        at shooting node 0.
        Type: :code:`np.ndarray`; default: :code:`np.array([])`.
        """
        return self.__lphi_0

    @lphi_0.setter
    def lphi_0(self, value):
        self.__lphi_0 = cast_to_1d_nparray(value, 'lphi_0')

    @property
    def uphi_0(self):
        r""":math:`\bar{\phi}^0` - upper bound on convex-over-nonlinear inequalities
        at shooting node 0.
        Type: :code:`np.ndarray`; default: :code:`np.array([])`.
        """
        return self.__uphi_0


    @uphi_0.setter
    def uphi_0(self, value):
        self.__uphi_0 = cast_to_1d_nparray(value, 'uphi_0')

    @property
    def idxs_rev_0(self):
        """Indices of slack variables associated with each constraint at initial shooting node 0, zero-based.
        Type: :code:`np.ndarray`; default: :code:`np.array([])`.
        """
        return self.__idxs_rev_0

    @idxs_rev_0.setter
    def idxs_rev_0(self, idxs_rev_0):
        self.__idxs_rev_0 = cast_to_1d_nparray(idxs_rev_0, "idxs_rev_0")

    @property
    def idxs_rev(self):
        """Indices of slack variables associated with each constraint at shooting nodes (1 to N-1), zero-based.
        Type: :code:`np.ndarray`; default: :code:`np.array([])`.
        """
        return self.__idxs_rev

    @idxs_rev.setter
    def idxs_rev(self, idxs_rev):
        self.__idxs_rev = cast_to_1d_nparray(idxs_rev, "idxs_rev")

    @property
    def idxs_rev_e(self):
        """Indices of slack variables associated with each constraint at terminal shooting node N, zero-based.
        Type: :code:`np.ndarray`; default: :code:`np.array([])`.
        """
        return self.__idxs_rev_e

    @idxs_rev_e.setter
    def idxs_rev_e(self, idxs_rev_e):
        self.__idxs_rev_e = cast_to_1d_nparray(idxs_rev_e, "idxs_rev_e")

    @property
    def ls_0(self):
        """Lower bounds on slacks associated with lower bound constraints at initial shooting node 0.
        Type: :code:`np.ndarray`; default: :code:`np.array([])`.
        """
        return self.__ls_0

    @ls_0.setter
    def ls_0(self, ls_0):
        self.__ls_0 = cast_to_1d_nparray(ls_0, "ls_0")

    @property
    def ls(self):
        """Lower bounds on slacks associated with lower bound constraints at shooting nodes (1 to N-1).
        Type: :code:`np.ndarray`; default: :code:`np.array([])`.
        """
        return self.__ls

    @ls.setter
    def ls(self, ls):
        self.__ls = cast_to_1d_nparray(ls, "ls")

    @property
    def ls_e(self):
        """Lower bounds on slacks associated with lower bound constraints at terminal shooting node N.
        Type: :code:`np.ndarray`; default: :code:`np.array([])`.
        """
        return self.__ls_e

    @ls_e.setter
    def ls_e(self, ls_e):
        self.__ls_e = cast_to_1d_nparray(ls_e, "ls_e")

    @property
    def us_0(self):
        """Lower bounds on slacks associated with upper bound constraints at initial shooting node 0.
        Type: :code:`np.ndarray`; default: :code:`np.array([])`.
        """
        return self.__us_0

    @us_0.setter
    def us_0(self, us_0):
        self.__us_0 = cast_to_1d_nparray(us_0, "us_0")

    @property
    def us(self):
        """Lower bounds on slacks associated with upper bound constraints at shooting nodes (1 to N-1).
        Type: :code:`np.ndarray`; default: :code:`np.array([])`.
        """
        return self.__us

    @us.setter
    def us(self, us):
        self.__us = cast_to_1d_nparray(us, "us")

    @property
    def us_e(self):
        """Lower bounds on slacks associated with upper bound constraints at terminal shooting node N.
        Type: :code:`np.ndarray`; default: :code:`np.array([])`.
        """
        return self.__us_e

    @us_e.setter
    def us_e(self, us_e):
        self.__us_e = cast_to_1d_nparray(us_e, "us_e")

    @property
    def lsbx(self):
        """Lower bounds on slacks corresponding to soft lower bounds on x
        at stages (1 to N-1);
        not required - zeros by default"""
        return self.__lsbx

    @lsbx.setter
    def lsbx(self, value):
        self.__lsbx = cast_to_1d_nparray(value, 'lsbx')

    @property
    def usbx(self):
        """Lower bounds on slacks corresponding to soft upper bounds on x
        at stages (1 to N-1);
        not required - zeros by default"""
        return self.__usbx

    @usbx.setter
    def usbx(self, value):
        self.__usbx = cast_to_1d_nparray(value, 'usbx')

    @property
    def idxsbx(self):
        """Indices of soft bounds on x within the indices of bounds on x
        at stages (1 to N-1).
        Can be set by using :py:attr:`Jsbx`.
        Type: :code:`np.ndarray`; default: :code:`np.array([])`"""
        return self.__idxsbx

    @idxsbx.setter
    def idxsbx(self, idxsbx):
        self.__idxsbx = cast_to_1d_nparray(idxsbx, "idxsbx")

    @property
    def Jsbx(self):
        """:math:`J_{sbx}` - matrix coefficient for soft bounds on x
        at stages (1 to N-1);
        Translated internally into :py:attr:`idxsbx`."""
        print_J_to_idx_note()
        return self.__idxsbx

    @Jsbx.setter
    def Jsbx(self, Jsbx):
        Jsbx = cast_to_2d_nparray(Jsbx, "Jsbx")
        self.__idxsbx = J_to_idx_slack(Jsbx)

    @property
    def lsbu(self):
        """Lower bounds on slacks corresponding to soft lower bounds on u
        at stages (0 to N-1).
        Not required - zeros by default."""
        return self.__lsbu

    @lsbu.setter
    def lsbu(self, value):
        self.__lsbu = cast_to_1d_nparray(value, 'lsbu')

    @property
    def usbu(self):
        """Lower bounds on slacks corresponding to soft upper bounds on u
        at stages (0 to N-1);
        not required - zeros by default"""
        return self.__usbu

    @usbu.setter
    def usbu(self, value):
        self.__usbu = cast_to_1d_nparray(value, 'usbu')

    @property
    def idxsbu(self):
        """Indices of soft bounds on u within the indices of bounds on u
        at stages (0 to N-1).
        Can be set by using :py:attr:`Jsbu`.
        Type: :code:`np.ndarray`; default: :code:`np.array([])`"""
        return self.__idxsbu

    @idxsbu.setter
    def idxsbu(self, idxsbu):
        self.__idxsbu = cast_to_1d_nparray(idxsbu, "idxsbu")

    @property
    def Jsbu(self):
        """:math:`J_{sbu}` - matrix coefficient for soft bounds on u
        at stages (0 to N-1);
        internally translated into :py:attr:`idxsbu`"""
        print_J_to_idx_note()
        return self.__idxsbu

    @Jsbu.setter
    def Jsbu(self, Jsbu):
        Jsbu = cast_to_2d_nparray(Jsbu, "Jsbu")
        self.__idxsbu = J_to_idx_slack(Jsbu)

    @property
    def lsbx_e(self):
        """Lower bounds on slacks corresponding to soft lower bounds on x at shooting node N.
        Not required - zeros by default"""
        return self.__lsbx_e

    @lsbx_e.setter
    def lsbx_e(self, value):
        self.__lsbx_e = cast_to_1d_nparray(value, 'lsbx_e')

    @property
    def usbx_e(self):
        """Lower bounds on slacks corresponding to soft upper bounds on x at shooting node N.
        Not required - zeros by default"""
        return self.__usbx_e

    @usbx_e.setter
    def usbx_e(self, value):
        self.__usbx_e = cast_to_1d_nparray(value, 'usbx_e')

    @property
    def idxsbx_e(self):
        """Indices of soft bounds on x at shooting node N, within the indices of bounds on x at terminal shooting node N.
        Can be set by using :py:attr:`Jsbx_e`.
        Type: :code:`np.ndarray`; default: :code:`np.array([])`"""
        return self.__idxsbx_e

    @idxsbx_e.setter
    def idxsbx_e(self, idxsbx_e):
        self.__idxsbx_e = cast_to_1d_nparray(idxsbx_e, "idxsbx_e")

    @property
    def Jsbx_e(self):
        """:math:`J_{sbx}^e` - matrix coefficient for soft bounds on x at terminal shooting node N.
        Translated internally to :py:attr:`idxsbx_e`"""
        print_J_to_idx_note()
        return self.__idxsbx_e

    @Jsbx_e.setter
    def Jsbx_e(self, Jsbx_e):
        Jsbx_e = cast_to_2d_nparray(Jsbx_e, "Jsbx_e")
        self.__idxsbx_e = J_to_idx_slack(Jsbx_e)

    @property
    def lsg(self):
        """Lower bounds on slacks corresponding to soft lower bounds for general linear constraints
        at stages (0 to N-1).
        Type: :code:`np.ndarray`; default: :code:`np.array([])`
        """
        return self.__lsg

    @lsg.setter
    def lsg(self, value):
        self.__lsg = cast_to_1d_nparray(value, 'lsg')

    @property
    def usg(self):
        """Lower bounds on slacks corresponding to soft upper bounds for general linear constraints.
        Not required - zeros by default"""
        return self.__usg

    @usg.setter
    def usg(self, value):
        self.__usg = cast_to_1d_nparray(value, 'usg')

    @property
    def idxsg(self):
        """Indices of soft general linear constraints within the indices of general linear constraints.
        Can be set by using :py:attr:`Jsg`.
        Type: :code:`np.ndarray`; default: :code:`np.array([])`"""
        return self.__idxsg

    @idxsg.setter
    def idxsg(self, value):
        self.__idxsg = cast_to_1d_nparray(value, 'idxsg')

    @property
    def Jsg(self):
        """:math:`J_{sg}` - matrix coefficient for soft bounds on general linear constraints.
        Translated internally to :py:attr:`idxsg`"""
        print_J_to_idx_note()
        return self.__idxsg

    @Jsg.setter
    def Jsg(self, Jsg):
        Jsg = cast_to_2d_nparray(Jsg, "Jsg")
        self.__idxsg = J_to_idx_slack(Jsg)

    @property
    def lsh(self):
        """Lower bounds on slacks corresponding to soft lower bounds for nonlinear constraints.
        Not required - zeros by default"""
        return self.__lsh

    @lsh.setter
    def lsh(self, value):
        self.__lsh = cast_to_1d_nparray(value, 'lsh')

    @property
    def ush(self):
        """Lower bounds on slacks corresponding to soft upper bounds for nonlinear constraints.
        Not required - zeros by default"""
        return self.__ush

    @ush.setter
    def ush(self, value):
        self.__ush = cast_to_1d_nparray(value, 'ush')

    @property
    def idxsh(self):
        """Indices of soft nonlinear constraints within the indices of nonlinear constraints.
        Can be set by using :py:attr:`Jbx`.
        Type: :code:`np.ndarray`; default: :code:`np.array([])`"""
        return self.__idxsh

    @idxsh.setter
    def idxsh(self, value):
        self.__idxsh = cast_to_1d_nparray(value, 'idxsh')


    @property
    def Jsh(self):
        """:math:`J_{sh}` - matrix coefficient for soft bounds on nonlinear constraints.
        Translated internally to :py:attr:`idxsh`"""
        print_J_to_idx_note()
        return self.__idxsh

    @Jsh.setter
    def Jsh(self, Jsh):
        Jsh = cast_to_2d_nparray(Jsh, "Jsh")
        self.__idxsh = J_to_idx_slack(Jsh)

    @property
    def lsphi(self):
        """Lower bounds on slacks corresponding to soft lower bounds for convex-over-nonlinear constraints.
        Not required - zeros by default"""
        return self.__lsphi

    @lsphi.setter
    def lsphi(self, value):
        self.__lsphi = cast_to_1d_nparray(value, 'lsphi')

    @property
    def usphi(self):
        """Lower bounds on slacks corresponding to soft upper bounds for convex-over-nonlinear constraints.
        Not required - zeros by default"""
        return self.__usphi

    @usphi.setter
    def usphi(self, value):
        self.__usphi = cast_to_1d_nparray(value, 'usphi')

    @property
    def idxsphi(self):
        """Indices of soft convex-over-nonlinear constraints within the indices of nonlinear constraints.
        Can be set by using :py:attr:`Jsphi`.
        Type: :code:`np.ndarray`; default: :code:`np.array([])`"""
        return self.__idxsphi

    @idxsphi.setter
    def idxsphi(self, value):
        self.__idxsphi = cast_to_1d_nparray(value, 'idxsphi')

    @property
    def Jsphi(self):
        r""":math:`J_{s, \phi}` - matrix coefficient for soft bounds on convex-over-nonlinear constraints.
        Translated internally into :py:attr:`idxsphi`."""
        print_J_to_idx_note()
        return self.__idxsphi


    @Jsphi.setter
    def Jsphi(self, Jsphi):
        Jsphi = cast_to_2d_nparray(Jsphi, "Jsphi")
        self.__idxsphi = J_to_idx_slack(Jsphi)

    @property
    def lsg_e(self):
        """Lower bounds on slacks corresponding to soft lower bounds for general linear constraints at shooting node N.
        Not required - zeros by default"""
        return self.__lsg_e

    @lsg_e.setter
    def lsg_e(self, value):
        self.__lsg_e = cast_to_1d_nparray(value, 'lsg_e')

    @property
    def usg_e(self):
        """Lower bounds on slacks corresponding to soft upper bounds for general linear constraints at shooting node N.
        Not required - zeros by default"""
        return self.__usg_e

    @usg_e.setter
    def usg_e(self, value):
        self.__usg_e = cast_to_1d_nparray(value, 'usg_e')

    @property
    def idxsg_e(self):
        """Indices of soft general linear constraints at shooting node N within the indices of general linear constraints at shooting node N.
        Can be set by using :py:attr:`Jsg_e`."""
        return self.__idxsg_e

    @idxsg_e.setter
    def idxsg_e(self, value):
        self.__idxsg_e = cast_to_1d_nparray(value, 'idxsg_e')

    @property
    def Jsg_e(self):
        """:math:`J_{s,h}^e` - matrix coefficient for soft bounds on general linear constraints at terminal shooting node N.
        Translated internally to :py:attr:`idxsg_e`"""
        print_J_to_idx_note()
        return self.__idxsg_e


    @Jsg_e.setter
    def Jsg_e(self, Jsg_e):
        Jsg_e = cast_to_2d_nparray(Jsg_e, "Jsg_e")
        self.__idxsg_e = J_to_idx_slack(Jsg_e)

    @property
    def lsh_0(self):
        """Lower bounds on slacks corresponding to soft lower bounds for nonlinear constraints at initial shooting node 0.
        Not required - zeros by default"""
        return self.__lsh_0

    @lsh_0.setter
    def lsh_0(self, value):
        self.__lsh_0 = cast_to_1d_nparray(value, 'lsh_0')

    @property
    def ush_0(self):
        """Lower bounds on slacks corresponding to soft upper bounds for nonlinear constraints at initial shooting node 0.
        Not required - zeros by default"""
        return self.__ush_0

    @ush_0.setter
    def ush_0(self, value):
        self.__ush_0 = cast_to_1d_nparray(value, 'ush_0')

    @property
    def idxsh_0(self):
        """Indices of soft nonlinear constraints at shooting node N within the indices of nonlinear constraints at initial shooting node 0.
        Can be set by using :py:attr:`Jsh_0`."""
        return self.__idxsh_0

    @idxsh_0.setter
    def idxsh_0(self, value):
        self.__idxsh_0 = cast_to_1d_nparray(value, 'idxsh_0')

    @property
    def Jsh_0(self):
        """:math:`J_{s,h}^0` - matrix coefficient for soft bounds on nonlinear constraints at initial shooting node 0; fills :py:attr:`idxsh_0`"""
        print_J_to_idx_note()
        return self.__idxsh_0

    @Jsh_0.setter
    def Jsh_0(self, Jsh_0):
        Jsh_0 = cast_to_2d_nparray(Jsh_0, "Jsh_0")
        self.__idxsh_0 = J_to_idx_slack(Jsh_0)

    @property
    def lsphi_0(self):
        """Lower bounds on slacks corresponding to soft lower bounds for convex-over-nonlinear constraints at initial shooting node 0.
        Not required - zeros by default"""
        return self.__lsphi_0

    @lsphi_0.setter
    def lsphi_0(self, value):
        self.__lsphi_0 = cast_to_1d_nparray(value, 'lsphi_0')

    @property
    def usphi_0(self):
        """Lower bounds on slacks corresponding to soft upper bounds for convex-over-nonlinear constraints at initial shooting node 0.
        Not required - zeros by default"""
        return self.__usphi_0

    @usphi_0.setter
    def usphi_0(self, value):
        self.__usphi_0 = cast_to_1d_nparray(value, 'usphi_0')

    @property
    def idxsphi_0(self):
        """Indices of soft nonlinear constraints at shooting node N within the indices of nonlinear constraints at initial shooting node 0.
        Can be set by using :py:attr:`Jsphi_0`.
        Type: :code:`np.ndarray`; default: :code:`np.array([])`"""
        return self.__idxsphi_0

    @idxsphi_0.setter
    def idxsphi_0(self, value):
        self.__idxsphi_0 = cast_to_1d_nparray(value, 'idxsphi_0')

    @property
    def Jsphi_0(self):
        """:math:`J_{sh}^0` - matrix coefficient for soft bounds on convex-over-nonlinear constraints at shooting node N.
        Translated internally to :py:attr:`idxsphi_0`"""
        print_J_to_idx_note()
        return self.__idxsphi_0


    @Jsphi_0.setter
    def Jsphi_0(self, Jsphi_0):
        Jsphi_0 = cast_to_2d_nparray(Jsphi_0, "Jsphi_0")
        self.__idxsphi_0 = J_to_idx_slack(Jsphi_0)

    @property
    def lsh_e(self):
        """Lower bounds on slacks corresponding to soft lower bounds for nonlinear constraints at terminal shooting node N.
        Not required - zeros by default"""
        return self.__lsh_e

    @lsh_e.setter
    def lsh_e(self, value):
        self.__lsh_e = cast_to_1d_nparray(value, 'lsh_e')

    @property
    def ush_e(self):
        """Lower bounds on slacks corresponding to soft upper bounds for nonlinear constraints at terminal shooting node N.
        Not required - zeros by default"""
        return self.__ush_e

    @ush_e.setter
    def ush_e(self, value):
        self.__ush_e = cast_to_1d_nparray(value, 'ush_e')

    @property
    def idxsh_e(self):
        """Indices of soft nonlinear constraints at shooting node N within the indices of nonlinear constraints at terminal shooting node N.
        Can be set by using :py:attr:`Jsh_e`."""
        return self.__idxsh_e

    @idxsh_e.setter
    def idxsh_e(self, value):
        self.__idxsh_e = cast_to_1d_nparray(value, 'idxsh_e')

    @property
    def Jsh_e(self):
        """:math:`J_{s,h}^e` - matrix coefficient for soft bounds on nonlinear constraints at terminal shooting node N; fills :py:attr:`idxsh_e`"""
        print_J_to_idx_note()
        return self.__idxsh_e

    @Jsh_e.setter
    def Jsh_e(self, Jsh_e):
        Jsh_e = cast_to_2d_nparray(Jsh_e, "Jsh_e")
        self.__idxsh_e = J_to_idx_slack(Jsh_e)


    @property
    def lsphi_e(self):
        """Lower bounds on slacks corresponding to soft lower bounds for convex-over-nonlinear constraints at terminal shooting node N.
        Not required - zeros by default"""
        return self.__lsphi_e

    @lsphi_e.setter
    def lsphi_e(self, value):
        self.__lsphi_e = cast_to_1d_nparray(value, 'lsphi_e')

    @property
    def usphi_e(self):
        """Lower bounds on slacks corresponding to soft upper bounds for convex-over-nonlinear constraints at terminal shooting node N.
        Not required - zeros by default"""
        return self.__usphi_e

    @usphi_e.setter
    def usphi_e(self, value):
        self.__usphi_e = cast_to_1d_nparray(value, 'usphi_e')

    @property
    def idxsphi_e(self):
        """Indices of soft nonlinear constraints at shooting node N within the indices of nonlinear constraints at terminal shooting node N.
        Can be set by using :py:attr:`Jsphi_e`.
        Type: :code:`np.ndarray`; default: :code:`np.array([])`"""
        return self.__idxsphi_e

    @idxsphi_e.setter
    def idxsphi_e(self, value):
        self.__idxsphi_e = cast_to_1d_nparray(value, 'idxsphi_e')

    @property
    def Jsphi_e(self):
        """:math:`J_{sh}^e` - matrix coefficient for soft bounds on convex-over-nonlinear constraints at shooting node N.
        Translated internally to :py:attr:`idxsphi_e`"""
        print_J_to_idx_note()
        return self.__idxsphi_e

    @Jsphi_e.setter
    def Jsphi_e(self, Jsphi_e):
        Jsphi_e = cast_to_2d_nparray(Jsphi_e, "Jsphi_e")
        self.__idxsphi_e = J_to_idx_slack(Jsphi_e)

    @property
    def x0(self):
        r"""
        :math:`x_0 \in \mathbb{R}^{n_x}` - initial state --
        Translated internally to :py:attr:`idxbx_0`, :py:attr:`lbx_0`, :py:attr:`ubx_0`, :py:attr:`idxbxe_0`
        """
        if self.has_x0:
            return self.lbx_0
        else:
            print("x0 is not set. You can set it or specify lbx_0, ubx_0, idxbx_0, idxbxe_0 to implement general bounds on x0.")
            print("")
            print("idxbx_0: ", self.__idxbx_0)
            print("lbx_0: ", self.__lbx_0)
            print("ubx_0: ", self.__ubx_0)
            print("idxbxe_0: ", self.__idxbxe_0)
            return None

    @x0.setter
    def x0(self, x0):
        if is_empty(x0):
            self.__has_x0 = False
            self.__lbx_0 = np.array([])
            self.__ubx_0 = np.array([])
            self.__idxbx_0 = np.array([])
            self.__idxbxe_0 = np.array([])
        else:
            x0 = cast_to_1d_nparray(x0, "x0")
            self.__lbx_0 = x0
            self.__ubx_0 = x0
            self.__idxbx_0 = np.arange(x0.size)
            self.__idxbxe_0 = np.arange(x0.size)
            self.__has_x0 = True

    @property
    def has_x0(self):
        """
        Internal variable to check if x0 is set.
        Cannot be set from outside.
        :bool: True if x0 is set, False otherwise.
        """
        return self.__has_x0

    def remove_x0_elimination(self):
        """Remove the elimination of x0 from the constraints, bounds on x0 are handled as general bounds on x."""
        self.__has_x0 = False
        self.idxbxe_0 = np.array([])

    def set(self, attr, value):
        setattr(self, attr, value)
