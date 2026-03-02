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
import casadi as ca

from .acados_ocp_qp import AcadosOcpQp


def _dm(arr) -> ca.DM:
    """Convert a numpy array or list to ca.DM."""
    return ca.DM(np.asarray(arr, dtype=float))


class AcadosCasadiOcpQp:
    """
    CasADi QP formulation of an OCP-structured QP given as an :class:`AcadosOcpQp`.

    Builds the standard NLP/QP form::

        min   f(w)
        s.t.  lbg <= g(w) <= ubg
              lbw <= w     <= ubw

    where the primal variable is ordered stage-by-stage::

        w = [x_0, u_0, sl_0, su_0,  x_1, u_1, sl_1, su_1,  ...,  x_N, sl_N, su_N]

    **Cost** (HPIPM convention: [u;x] ordering, cross-term S is (nu,nx))::

        f = sum_{i=0}^{N}  1/2 [u_i;x_i]' [R_i, S_i; S_i', Q_i] [u_i;x_i]
                         + [r_i;q_i]' [u_i;x_i]
                         + 1/2 sl_i' diag(Zl_i) sl_i  + zl_i' sl_i
                         + 1/2 su_i' diag(Zu_i) su_i  + zu_i' su_i

    **Constraints** assembled in g:

    * Dynamics equalities  ``x_{i+1} = A_i x_i + B_i u_i + b_i``  (lbg = ubg = 0)
    * Hard bound constraints on u/x                                 (via lbw / ubw)
    * Soft bound constraints with slack variables                   (in g)
    * Hard general linear constraints ``C_i x_i + D_i u_i``        (in g)
    * Soft general linear constraints with slack variables          (in g)
    * Slack lower bounds ``sl >= lls``, ``su >= lus``               (via lbw)
    * Equality bounds (idxe)                                        (lbw = ubw)
    * Constraint masks: inactive constraints use lbg=-inf/ubg=+inf  (or lbw=-inf/ubw=+inf)

    Access via :attr:`nlp`, :attr:`bounds`, :attr:`w0`, :attr:`index_map`.
    """

    def __init__(self, qp: AcadosOcpQp):

        qp.make_consistent()

        N = qp.N
        dims = qp.dims

        self._index_map = {
            'x_in_w': [],
            'u_in_w': [],
            'sl_in_w': [],
            'su_in_w': [],

            'lam_bx_in_lam_w': [],
            'lam_bu_in_lam_w': [],
            'lam_sl_in_lam_w': [],
            'lam_su_in_lam_w': [],
            'pi_in_lam_g': [],
            'lam_g_in_lam_g': [[] for _ in range(N+1)],
            'lam_g_sl_in_lam_g': [[] for _ in range(N+1)],
            'lam_g_su_in_lam_g': [[] for _ in range(N+1)],
            'lam_bx_sl_in_lam_g': [[] for _ in range(N+1)],
            'lam_bx_su_in_lam_g': [[] for _ in range(N+1)],
            'lam_bu_sl_in_lam_g': [[] for _ in range(N+1)],
            'lam_bu_su_in_lam_g': [[] for _ in range(N+1)],
        }
        # QP formulation
        w_sym_list = []
        lbw_list = []
        ubw_list = []
        w0_list = []

        # per-stage ca.SX symbol
        x_nodes_list = []
        u_nodes_list = []
        sl_nodes_list = []
        su_nodes_list = []
        # per-stage bounds list
        lbx_list = []
        ubx_list = []
        lbu_list = []
        ubu_list = []
        sl_lb_list = []
        sl_ub_list = []
        su_lb_list = []
        su_ub_list = []

        self.offset_w = 0
        self.offset_g = 0
        self.offset_lam = 0

        # 1. Create symbolic variables and record their w-indices
        for i in range(N+1):
            # First set all bounds for this stage
            self._setup_bounds('x', i, lbx_list, ubx_list, qp, dims)
            self._setup_bounds('u', i, lbu_list, ubu_list, qp, dims)
            self._setup_bounds('sl', i, sl_lb_list, sl_ub_list, qp, dims)
            self._setup_bounds('su', i, su_lb_list, su_ub_list, qp, dims)
            # Then create symbolics and append for this stage
            self._create_symbolics_and_append('x', x_nodes_list, w_sym_list, lbw_list, ubw_list, w0_list, lbx_list, ubx_list, i, qp.dims)
            self._create_symbolics_and_append('u', u_nodes_list, w_sym_list, lbw_list, ubw_list, w0_list, lbu_list, ubu_list, i, qp.dims)
            self._create_symbolics_and_append('sl', sl_nodes_list, w_sym_list, lbw_list, ubw_list, w0_list, sl_lb_list, sl_ub_list, i, qp.dims)
            self._create_symbolics_and_append('su', su_nodes_list, w_sym_list, lbw_list, ubw_list, w0_list, su_lb_list, su_ub_list, i, qp.dims)

        # Assemble
        w = ca.vertcat(*w_sym_list)
        lbw = np.concatenate(lbw_list).flatten()
        ubw = np.concatenate(ubw_list).flatten()
        w0 = np.concatenate(w0_list).flatten()

        # 2. Constraints: assemble into a single g with corresponding lbg/ubg, and record indices
        g_list = []
        lbg_list = []
        ubg_list = []

        # --- dynamics equalities: x_{i+1} = A_i x_i + B_i u_i + b_i ---
        for i in range(N):
            nx_next = int(dims.nx[i + 1])
            dyn = (x_nodes_list[i + 1]
                       - ca.mtimes(_dm(qp.A[i]), x_nodes_list[i])
                       - ca.mtimes(_dm(qp.B[i]), u_nodes_list[i])
                       - _dm(qp.b[i]))
            self._append_constraints(dyn, np.zeros(nx_next), np.zeros(nx_next), g_list, lbg_list, ubg_list, i, type='dyn')

        # --- bound constraints that are SOFT ---
        # and general linear constraints (hard and soft)
        for i in range(N + 1):
            nu = int(dims.nu[i])
            nb = int(dims.nb[i])
            ng = int(dims.ng[i])
            ns = int(dims.ns[i])
            nbu = int(dims.nbu[i])
            nbx = int(dims.nbx[i])

            idxs_rev = qp.idxs_rev[i]

            xi = x_nodes_list[i]
            ui = u_nodes_list[i]
            sli = sl_nodes_list[i]
            sui = su_nodes_list[i]

            lbu_mask = np.asarray(qp.lbu_mask[i]).reshape(-1) if qp.lbu_mask[i] is not None else np.ones(nbu)
            ubu_mask = np.asarray(qp.ubu_mask[i]).reshape(-1) if qp.ubu_mask[i] is not None else np.ones(nbu)
            lbx_mask = np.asarray(qp.lbx_mask[i]).reshape(-1) if qp.lbx_mask[i] is not None else np.ones(nbx)
            ubx_mask = np.asarray(qp.ubx_mask[i]).reshape(-1) if qp.ubx_mask[i] is not None else np.ones(nbx)
            lg_mask = np.asarray(qp.lg_mask[i]).reshape(-1) if qp.lg_mask[i] is not None else np.ones(ng)
            ug_mask = np.asarray(qp.ug_mask[i]).reshape(-1) if qp.ug_mask[i] is not None else np.ones(ng)

            # soft u bounds
            idxsbu_rev = idxs_rev[:nbu]
            valid_bool = idxsbu_rev >= 0
            soft_idx = np.where(valid_bool)[0]
            s_idx = idxs_rev[soft_idx] # value in idxsbu_rev is indice of slack variable in w
            if s_idx.size > 0:
                expr = ui[soft_idx]
                sl_k = sli[s_idx]
                su_k = sui[s_idx]
                # unmasked and soft constraints are handled in g
                lbu_full = -np.inf * np.ones((nbu, 1))
                ubu_full =  np.inf * np.ones((nbu, 1))
                lb_valid_bool = (lbu_mask > 0) & (idxsbu_rev >= 0)
                ub_valid_bool = (ubu_mask > 0) & (idxsbu_rev >= 0)
                lb_soft_idx = np.where(lb_valid_bool)[0]
                ub_soft_idx = np.where(ub_valid_bool)[0]
                lbu_full[lb_soft_idx] = qp.lbu[i][lb_soft_idx].reshape(-1, 1)
                ubu_full[ub_soft_idx] = qp.ubu[i][ub_soft_idx].reshape(-1, 1)
                self._append_constraints(expr+sl_k, lbu_full, np.full(lbu_full.shape, np.inf), g_list, lbg_list, ubg_list, i, type='bu_soft_lower')
                self._append_constraints(expr-su_k, np.full(ubu_full.shape, -np.inf), ubu_full, g_list, lbg_list, ubg_list, i, type='bu_soft_upper')

            # soft x bounds
            idxsxb_rev = idxs_rev[nbu:nbu+nbx]
            valid_bool = idxsxb_rev >= 0
            soft_idx = np.where(valid_bool)[0]
            s_idx = idxs_rev[nbu + soft_idx] # value in idxsxb_rev is indice of slack variable
            if s_idx.size > 0:
                expr = xi[soft_idx]
                sl_k = sli[s_idx]
                su_k = sui[s_idx]
                # unmasked and soft constraints are handled in g
                lbx_full = -np.inf * np.ones((nbx, 1))
                ubx_full =  np.inf * np.ones((nbx, 1))
                lb_valid_bool = (lbx_mask > 0) & (idxsxb_rev >= 0)
                ub_valid_bool = (ubx_mask > 0) & (idxsxb_rev >= 0)
                lb_soft_idx = np.where(lb_valid_bool)[0]
                ub_soft_idx = np.where(ub_valid_bool)[0]
                lbx_full[lb_soft_idx] = qp.lbx[i][lb_soft_idx].reshape(-1, 1)
                ubx_full[ub_soft_idx] = qp.ubx[i][ub_soft_idx].reshape(-1, 1)
                self._append_constraints(expr+sl_k, lbx_full, np.full(lbx_full.shape, np.inf), g_list, lbg_list, ubg_list, i, type='bx_soft_lower')
                self._append_constraints(expr-su_k, np.full(ubx_full.shape, -np.inf), ubx_full, g_list, lbg_list, ubg_list, i, type='bx_soft_upper')

            # general linear constraints: C x + D u
            if ng > 0:
                C_i = _dm(qp.C[i])
                D_i = _dm(qp.D[i])
                lg_i = np.asarray(qp.lg[i]).reshape(-1)
                ug_i = np.asarray(qp.ug[i]).reshape(-1)
                lin_expr = (ca.mtimes(C_i, xi) + ca.mtimes(D_i, ui)
                            if nu > 0 else ca.mtimes(C_i, xi))

                for j in range(ng):
                    s_idx = idxs_rev[nb + j]
                    lb_j = lg_i[j] if lg_mask[j] > 0 else -np.inf
                    ub_j = ug_i[j] if ug_mask[j] > 0 else np.inf
                    expr = lin_expr[j]

                    if s_idx >= 0:
                        # soft: split into two one-sided constraints with slack
                        sl_k = sli[s_idx]
                        su_k = sui[s_idx]
                        self._append_constraints(expr+sl_k, lb_j, np.inf, g_list, lbg_list, ubg_list, i, type='lg_soft_lower')
                        self._append_constraints(expr-su_k, -np.inf, ub_j, g_list, lbg_list, ubg_list, i, type='lg_soft_upper')
                    else:
                        # hard: single double-sided constraint
                        self._append_constraints(expr, lb_j, ub_j, g_list, lbg_list, ubg_list, i, type='hard')

        # 3. Cost
        f = ca.SX(0)

        for i in range(N + 1):
            xi = x_nodes_list[i]
            ui = u_nodes_list[i]
            sli = sl_nodes_list[i]
            sui = su_nodes_list[i]
            nu = int(dims.nu[i])
            ns = int(dims.ns[i])

            Q_i = _dm(qp.Q[i])
            q_i = _dm(qp.q[i])

            f += 0.5 * ca.bilin(Q_i, xi, xi) + ca.dot(q_i, xi)

            if nu > 0:
                R_i = _dm(qp.R[i])
                r_i = _dm(qp.r[i])
                S_i = _dm(qp.S[i])  # shape (nu, nx)

                # 1/2 [u;x]' [R,S;S',Q] [u;x]  →  1/2 u'Ru + u'Sx + 1/2 x'Qx
                f += 0.5 * ca.bilin(R_i, ui, ui)
                f += ca.mtimes(ui.T, ca.mtimes(S_i, xi))
                f += ca.dot(r_i, ui)

            if ns > 0:
                Zl_i = _dm(np.asarray(qp.Zl[i]).reshape(-1))
                Zu_i = _dm(np.asarray(qp.Zu[i]).reshape(-1))
                zl_i = _dm(np.asarray(qp.zl[i]).reshape(-1))
                zu_i = _dm(np.asarray(qp.zu[i]).reshape(-1))

                f += 0.5 * ca.bilin(ca.diag(Zl_i), sli, sli) + ca.dot(zl_i, sli)
                f += 0.5 * ca.bilin(ca.diag(Zu_i), sui, sui) + ca.dot(zu_i, sui)

        # 4. Assemble into NLP/QP form
        g = ca.vertcat(*g_list)
        lbg = np.concatenate(lbg_list)
        ubg = np.concatenate(ubg_list)

        self.__qp = {'x': w, 'f': f, 'g': g}
        self.__bounds = {'lbx': lbw, 'ubx': ubw, 'lbg': lbg, 'ubg': ubg}
        self.__w0 = w0

    def _setup_bounds(self, _field, i, lb_node_list, ub_node_list, qp: AcadosOcpQp, dims):

        nx  = dims.nx[i]
        nu  = dims.nu[i]
        nbu = dims.nbu[i]
        nbx = dims.nbx[i]
        ns  = dims.ns[i]
        idxs_rev = qp.idxs_rev[i]

        if _field == 'x':
            lbx_mask = np.asarray(qp.lbx_mask[i]).reshape(-1) if qp.lbx_mask[i] is not None else np.ones(nbx)
            ubx_mask = np.asarray(qp.ubx_mask[i]).reshape(-1) if qp.ubx_mask[i] is not None else np.ones(nbx)
            lbx_full = -np.inf * np.ones((nx, 1))
            ubx_full = np.inf * np.ones((nx, 1))
            idxsx_rev = idxs_rev[nbu:nbu+nbx]  # x bounds follow u bounds in idxs_rev
            lb_valid_bool = (lbx_mask > 0) & (idxsx_rev < 0)
            ub_valid_bool = (ubx_mask > 0) & (idxsx_rev < 0)
            lb_hard_idx = np.where(lb_valid_bool)[0]
            ub_hard_idx = np.where(ub_valid_bool)[0]
            lbx_full[lb_hard_idx] = qp.lbx[i][lb_hard_idx].reshape(-1, 1)
            ubx_full[ub_hard_idx] = qp.ubx[i][ub_hard_idx].reshape(-1, 1)
            lb_node_list.append(lbx_full)
            ub_node_list.append(ubx_full)
            self._index_map['lam_bx_in_lam_w'].append(list(range(self.offset_lam, self.offset_lam + nbx)))
            self.offset_lam += nx

        elif _field == 'u':
            lbu_mask = np.asarray(qp.lbu_mask[i]).reshape(-1)
            ubu_mask = np.asarray(qp.ubu_mask[i]).reshape(-1)
            idxsu_rev = idxs_rev[:nbu]  # x bounds follow u bounds in idxs_rev
            lbu_full = -np.inf * np.ones((nu, 1))
            ubu_full = np.inf * np.ones((nu, 1))
             # unmasked and hard, soft constraints are handled in g
            lb_valid_bool = (lbu_mask > 0) & (idxsu_rev < 0)
            ub_valid_bool = (ubu_mask > 0) & (idxsu_rev < 0)
            lb_hard_idx = np.where(lb_valid_bool)[0]
            ub_hard_idx = np.where(ub_valid_bool)[0]
            lbu_full[lb_hard_idx] = qp.lbu[i][lb_hard_idx].reshape(-1, 1)
            ubu_full[ub_hard_idx] = qp.ubu[i][ub_hard_idx].reshape(-1, 1)
            lb_node_list.append(lbu_full)
            ub_node_list.append(ubu_full)
            self._index_map['lam_bu_in_lam_w'].append(list(range(self.offset_lam, self.offset_lam + nbu)))
            self.offset_lam += nu

        elif _field == 'sl':
            lb_node_list.append(0 * np.ones((ns, 1)))
            ub_node_list.append(np.inf * np.ones((ns, 1)))
            self._index_map['lam_sl_in_lam_w'].append(list(range(self.offset_lam, self.offset_lam + ns)))
            self.offset_lam += ns

        elif _field == 'su':
            lb_node_list.append(0 * np.ones((ns, 1)))
            ub_node_list.append(np.inf * np.ones((ns, 1)))
            self._index_map['lam_su_in_lam_w'].append(list(range(self.offset_lam, self.offset_lam + ns)))
            self.offset_lam += ns

    def _create_symbolics_and_append(self, _field, node_list: list,
                                     w_sym_list, lbw_list, ubw_list, w0_list,
                                     lb_node_list, ub_node_list, 
                                     i, dims):
        nx = int(dims.nx[i])
        nu = int(dims.nu[i])
        ns = int(dims.ns[i])

        if _field == 'x':
            xi = ca.SX.sym(f'x_{i}', nx)
            node_list.append(xi)
            w_sym_list.append(xi)
            lbw_list.append(lb_node_list[i])
            ubw_list.append(ub_node_list[i])
            w0_list.append(np.zeros(nx))
            self._index_map['x_in_w'].append(list(range(self.offset_w, self.offset_w + nx)))
            self.offset_w += nx

        elif _field == 'u':
            if nu > 0:
                ui = ca.SX.sym(f'u_{i}', nu)
                node_list.append(ui)
                w_sym_list.append(ui)
                lbw_list.append(lb_node_list[i])
                ubw_list.append(ub_node_list[i])
                w0_list.append(np.zeros(nu))
                self._index_map['u_in_w'].append(list(range(self.offset_w, self.offset_w + nu)))
                self.offset_w += nu
            else:
                node_list.append(ca.SX.sym(f'u_{i}', 0))
                self._index_map['u_in_w'].append([])

        elif _field == 'sl':
            sli = ca.SX.sym(f'sl_{i}', ns)
            node_list.append(sli)
            w_sym_list.append(sli)
            lbw_list.append(lb_node_list[i])
            ubw_list.append(ub_node_list[i])
            w0_list.append(np.zeros(ns))
            self._index_map['sl_in_w'].append(list(range(self.offset_w, self.offset_w + ns)))
            self.offset_w += ns

        elif _field == 'su':
            sui = ca.SX.sym(f'su_{i}', ns)
            node_list.append(sui)
            w_sym_list.append(sui)
            lbw_list.append(lb_node_list[i])
            ubw_list.append(ub_node_list[i])
            w0_list.append(np.zeros(ns))
            self._index_map['su_in_w'].append(list(range(self.offset_w, self.offset_w + ns)))
            self.offset_w += ns

    def _append_constraints(self, expr, lb, ub, g_list, lbg_list, ubg_list, stage, type):
        if type == 'dyn':
            nx_next = expr.shape[0]
            g_list.append(expr)
            lbg_list.append(lb)
            ubg_list.append(ub)
            self._index_map['pi_in_lam_g'].append(list(range(self.offset_g, self.offset_g + nx_next)))
            self.offset_g += nx_next
        elif type == 'hard':
            g_list.append(expr)
            lbg_list.append(np.asarray(lb).reshape(-1))
            ubg_list.append(np.asarray(ub).reshape(-1))
            self._index_map['lam_g_in_lam_g'][stage].append(self.offset_g)
            self.offset_g += 1
        elif type == 'lg_soft_lower':
            g_list.append(expr)
            lbg_list.append(np.asarray(lb).reshape(-1))
            ubg_list.append(np.full((expr.shape[0],), np.inf))
            self._index_map['lam_g_sl_in_lam_g'][stage].append(self.offset_g)
            self.offset_g += 1
        elif type == 'lg_soft_upper':
            g_list.append(expr)
            lbg_list.append(np.full((expr.shape[0],), -np.inf))
            ubg_list.append(np.asarray(ub).reshape(-1))
            self._index_map['lam_g_su_in_lam_g'][stage].append(self.offset_g)
            self.offset_g += 1
        elif type == 'bx_soft_upper':
            g_list.append(expr)
            lbg_list.append(np.asarray(lb).reshape(-1))
            ubg_list.append(np.asarray(ub).reshape(-1))
            self._index_map['lam_bx_su_in_lam_g'][stage].append(self.offset_g)
            self.offset_g += 1
        elif type == 'bx_soft_lower':
            g_list.append(expr)
            lbg_list.append(np.asarray(lb).reshape(-1))
            ubg_list.append(np.asarray(ub).reshape(-1))
            self._index_map['lam_bx_sl_in_lam_g'][stage].append(self.offset_g)
            self.offset_g += 1
        elif type == 'bu_soft_upper':
            g_list.append(expr)
            lbg_list.append(np.asarray(lb).reshape(-1))
            ubg_list.append(np.asarray(ub).reshape(-1))
            self._index_map['lam_bu_su_in_lam_g'][stage].append(self.offset_g)
            self.offset_g += 1
        elif type == 'bu_soft_lower':
            g_list.append(expr)
            lbg_list.append(np.asarray(lb).reshape(-1))
            ubg_list.append(np.asarray(ub).reshape(-1))
            self._index_map['lam_bu_sl_in_lam_g'][stage].append(self.offset_g)
            self.offset_g += 1

    @property
    def qp(self) -> dict:
        """
        CasADi QP dictionary with keys ``'x'``, ``'f'``, ``'g'``.
        Pass directly to ``casadi.nlpsol`` or ``casadi.qpsol``.
        """
        return self.__qp

    @property
    def bounds(self) -> dict:
        """
        Bounds dictionary with keys ``'lbx'``, ``'ubx'``, ``'lbg'``, ``'ubg'``
        (numpy arrays).  Pass directly to the CasADi solver call.
        """
        return self.__bounds

    @property
    def w0(self) -> np.ndarray:
        """Default (zero) initial guess for the primal variable vector w."""
        return self.__w0

    @property
    def index_map(self) -> dict:
        """
        Dictionary mapping acados stage-wise quantities to indices in ``w`` / ``g``:
        """
        return self._index_map
