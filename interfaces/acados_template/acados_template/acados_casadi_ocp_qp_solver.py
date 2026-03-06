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

from typing import Union, List, Optional

import casadi as ca
import numpy as np

from .acados_ocp_qp import AcadosOcpQp
from .acados_ocp_iterate import AcadosOcpIterate, AcadosOcpFlattenedIterate
from .acados_casadi_ocp_qp import AcadosCasadiOcpQp


class AcadosCasadiOcpQpSolver:
    """
    CasADi-based reference NLP solver for OCP-structured QPs (:class:`AcadosOcpQp`).

    Wraps :class:`AcadosCasadiOcpQp` and a ``casadi.nlpsol`` backend (e.g. Ipopt).

    Typical usage::

        qp = AcadosOcpQp.from_json("my_qp.json")
        solver = AcadosCasadiOcpQpSolver(qp, solver="ipopt")
        u0 = solver.solve_for_x0(x0)

    :param qp:             The OCP-structured QP to solve.
    :param solver:         CasADi NLP solver name ('ipopt', 'fatrop'(?), …).
    :param verbose:        Whether to print solver output. (Not really implemented yet.)
    :param solver_opts:    Optional dictionary of solver options to pass to CasADi.
    """

    def __init__(
        self,
        qp: AcadosOcpQp,
        solver: str = "ipopt",
        verbose: bool = True,
        solver_opts: Optional[dict] = None,
    ):
        if not isinstance(qp, AcadosOcpQp):
            raise TypeError("qp must be an instance of AcadosOcpQp.")

        # build CasADi QP
        casadi_ocp_qp = AcadosCasadiOcpQp(qp)
        self.acados_ocp_qp = qp
        self.casadi_qp = casadi_ocp_qp.qp
        self.bounds = casadi_ocp_qp.bounds
        self.w0 = casadi_ocp_qp.w0
        self.index_map = casadi_ocp_qp.index_map

        if solver_opts is None:
            solver_opts = {}

        self._casadi_solver = ca.nlpsol('reference_solver', solver, self.casadi_qp, solver_opts)

        # dual cold-start
        nw = self.casadi_qp['x'].shape[0]
        ng = self.casadi_qp['g'].shape[0]
        self.lam_x0 = np.zeros(nw)
        self.lam_g0 = np.zeros(ng)

        # last solution
        self.qp_sol = None
        self.qp_sol_w = None
        self.qp_sol_lam_w = None
        self.qp_sol_lam_g = None
        self._status = None

    @property
    def status(self):
        """Status string of the last solve (e.g. ``"solved"``)."""
        return self._status

    def solve(self) -> str:
        """
        Solve the QP with the current bounds and warm-start.

        :return: solver status string.
        """
        sol = self._casadi_solver(
            x0=self.w0,
            lbx=self.bounds['lbx'],
            ubx=self.bounds['ubx'],
            lbg=self.bounds['lbg'],
            ubg=self.bounds['ubg'],
            lam_x0=self.lam_x0,
            lam_g0=self.lam_g0,
        )
        self.qp_sol = sol
        self.qp_sol_w = sol['x'].full().flatten()
        self.qp_sol_lam_w = sol['lam_x'].full().flatten()
        self.qp_sol_lam_g = sol['lam_g'].full().flatten() if self.casadi_qp['g'].shape[0] > 0 else np.empty(0)

        stats = self._casadi_solver.stats()
        self._status = stats.get('return_status', str(stats.get('success', 'unknown')))
        self._nlp_iter = stats.get('iter_count', None)
        self._time_total = stats.get('t_wall_total', None)
        return self._status

    def solve_for_x0(self, x0_bar: np.ndarray) -> np.ndarray:
        """
        Fix the initial state to *x0_bar*, solve the QP, and return *u_0*.

        :param x0_bar: initial state vector of length ``nx[0]``.
        :return: optimal control ``u_0``.
        """
        x0_bar = np.asarray(x0_bar).flatten()
        idx = self.index_map['x_in_w'][0]
        self.bounds['lbx'][idx] = x0_bar
        self.bounds['ubx'][idx] = x0_bar
        self.solve()
        return self.get(0, 'u')

    def get(self, stage: int, field: str) -> np.ndarray:
        """
        Retrieve a stage-wise quantity from the last solution.

        :return: numpy array.
        """
        if self.qp_sol is None:
            raise RuntimeError("No solution available; call solve() first.")
        if not isinstance(stage, int):
            raise TypeError("stage must be an integer.")

        w = self.qp_sol_w
        if field == 'x':
            return w[self.index_map['x_in_w'][stage]].copy()
        elif field == 'u':
            idx = self.index_map['u_in_w'][stage]
            return w[idx].copy() if idx else np.empty(0)
        elif field == 'sl':
            idx = self.index_map['sl_in_w'][stage]
            return w[idx].copy() if idx else np.empty(0)
        elif field == 'su':
            idx = self.index_map['su_in_w'][stage]
            return w[idx].copy() if idx else np.empty(0)
        elif field == 'pi':
            # dynamics equality multiplier: convention pi = -lam_g (acados sign)
            if stage >= self.acados_ocp_qp.N:
                return np.empty(0)
            return -self.qp_sol_lam_g[self.index_map['pi_in_lam_g'][stage]].copy()
        elif field == 'lam':
            return self._get_lam(stage)
        elif field == 'z':
            return np.empty((0,))
        else:
            raise NotImplementedError(f"get(): field '{field}' not implemented.")

    def _get_lam(self, stage: int) -> np.ndarray:
        """
        Return multipliers in acados convention:
        ``[lbu, lbx, lbg, ubu, ubx, ubg
           lsbu, lsbx, lsg, usbu, usbx, usg]``.
        """
        lam_w = self.qp_sol_lam_w
        lam_g = self.qp_sol_lam_g

        # --- hard u-bound multipliers ---
        bu_lam = lam_w[self.index_map['lam_bu_in_lam_w'][stage]]
        lbu_lam = np.maximum(0.0, -bu_lam)
        ubu_lam = np.maximum(0.0,  bu_lam)
        # --- hard x-bound multipliers ---
        bx_lam = lam_w[self.index_map['lam_bx_in_lam_w'][stage]]
        lbx_lam = np.maximum(0.0, -bx_lam)
        ubx_lam = np.maximum(0.0, bx_lam)
        # --- hard and soft constraint multipliers ---
        g_lam = lam_g[self.index_map['lam_g_in_lam_g'][stage]]
        if self.index_map['lam_g_su_in_lam_g'][stage] or self.index_map['lam_bx_su_in_lam_g'][stage] or self.index_map['lam_bu_su_in_lam_g'][stage]:
            lg_soft_lam = lam_g[self.index_map['lam_g_sl_in_lam_g'][stage] + 
                                self.index_map['lam_bx_sl_in_lam_g'][stage] + 
                                self.index_map['lam_bu_sl_in_lam_g'][stage]]
            ug_soft_lam = lam_g[self.index_map['lam_g_su_in_lam_g'][stage] + 
                                self.index_map['lam_bx_su_in_lam_g'][stage] + 
                                self.index_map['lam_bu_su_in_lam_g'][stage]]
            g_indices = np.array(self.index_map['lam_g_in_lam_g'][stage]+\
                                 self.index_map['lam_g_sl_in_lam_g'][stage]+\
                                 self.index_map['lam_bx_sl_in_lam_g'][stage]+\
                                 self.index_map['lam_bu_sl_in_lam_g'][stage])
            sorted_indices = np.argsort(g_indices) # get the right order
            g_lam_lower = np.concatenate((np.maximum(0, -g_lam), -lg_soft_lam))
            lbg_lam = g_lam_lower[sorted_indices]
            g_lam_upper = np.concatenate((np.maximum(0, g_lam), ug_soft_lam))
            ubg_lam = g_lam_upper[sorted_indices]
        else:
            lbg_lam = np.maximum(0, -g_lam)
            ubg_lam = np.maximum(0, g_lam)
        # --- slack variables multipliers ---
        # TODO: could be classified into bu, bx, g
        sl_lam = -lam_w[self.index_map['lam_sl_in_lam_w'][stage]]
        su_lam = -lam_w[self.index_map['lam_su_in_lam_w'][stage]]
        return np.concatenate((lbu_lam, lbx_lam, lbg_lam,
                                ubu_lam, ubx_lam, ubg_lam, 
                                sl_lam, su_lam)).flatten()

    def set(self, stage: int, field: str, value: np.ndarray):
        """
        Set a warm-start quantity.

        :param stage: shooting node index.
        :param field: ['x', 'u', 'pi', 'lam', 'sl', 'su'].
        :param value: numpy array.
        """
        if not isinstance(stage, int):
            raise TypeError("stage must be an integer.")
        v = np.asarray(value).flatten()

        if field == 'x':
            self.w0[self.index_map['x_in_w'][stage]] = v
        elif field == 'u':
            idx = self.index_map['u_in_w'][stage]
            if idx:
                self.w0[idx] = v
        elif field == 'sl':
            idx = self.index_map['sl_in_w'][stage]
            if idx:
                self.w0[idx] = v
        elif field == 'su':
            idx = self.index_map['su_in_w'][stage]
            if idx:
                self.w0[idx] = v
        elif field == 'pi':
            # dual warm-start for dynamics equality
            if stage < self.acados_ocp_qp.N:
                self.lam_g0[self.index_map['pi_in_lam_g'][stage]] = -v
        elif field == 'lam':
            self._set_lam(stage, v)
        else:
            raise NotImplementedError(f"set(): field '{field}' not implemented.")

    def _set_lam(self, stage: int, value: np.ndarray):
        """
        Warm-start dual variables from acados-convention lam vector
        ``[lbu, lbx, lbg, ubu, ubx, ubg
           lsbu, lsbx, lsg, usbu, usbx, usg ]``.
        """
        dims = self.acados_ocp_qp.dims
        nbu = dims.nbu[stage]
        nbx = dims.nbx[stage]
        ng_hard = len(self.index_map['lam_g_in_lam_g'][stage])
        nsg = len(self.index_map['lam_g_sl_in_lam_g'][stage])
        nsbx = len(self.index_map['lam_bx_sl_in_lam_g'][stage])
        nsbu = len(self.index_map['lam_bu_sl_in_lam_g'][stage])
        ng = ng_hard + nsg + nsbx + nsbu

        offset_upper = nbu + nbx + ng
        lbu_lam = value[:nbu]
        lbx_lam = value[nbu:nbu + nbx]
        lbg_lam = value[nbu + nbx:nbu + nbx + ng]
        ubu_lam = value[offset_upper:offset_upper + nbu]
        ubx_lam = value[offset_upper + nbu:offset_upper + nbu + nbx]
        ubg_lam = value[offset_upper + nbu + nbx:offset_upper + nbu + nbx + ng]
        # slacks
        offset_soft_l = 2 * (nbu + nbx + ng)
        offset_soft_u = offset_soft_l + nsg + nsbx + nsbu
        slbu_lam = value[offset_soft_l:offset_soft_l + nsbu]
        slbx_lam = value[offset_soft_l + nsbu:offset_soft_l + nsbu + nsbx]
        slg_lam = value[offset_soft_l + nsbu + nsbx:offset_soft_l + nsbu + nsbx + nsg]
        subu_lam = value[offset_soft_u:offset_soft_u + nsbu]
        subx_lam = value[offset_soft_u + nsbu:offset_soft_u + nsbu + nsbx]
        subg_lam = value[offset_soft_u + nsbu + nsbx:offset_soft_u + nsbu + nsbx + nsg]

        g_indices = np.array(
            self.index_map['lam_g_in_lam_g'][stage]
            + self.index_map['lam_g_sl_in_lam_g'][stage]
            + self.index_map['lam_bx_sl_in_lam_g'][stage]
            + self.index_map['lam_bu_sl_in_lam_g'][stage]
        )
        sorted_indices = np.argsort(g_indices)
        inverted_indices = np.argsort(sorted_indices)
        gnl_indices = inverted_indices[:len(self.index_map['lam_g_in_lam_g'][stage])]
        gnl_sl_indices = inverted_indices[len(self.index_map['lam_g_in_lam_g'][stage]):]
        lg_lam_hard = lbg_lam[gnl_indices]
        lg_lam_soft = lbg_lam[gnl_sl_indices]
        ug_lam_hard = ubg_lam[gnl_indices]
        ug_lam_soft = ubg_lam[gnl_sl_indices]

        if stage < self.acados_ocp_qp.N:
            self.lam_x0[self.index_map['lam_bu_in_lam_w'][stage]] = ubu_lam - lbu_lam
        self.lam_x0[self.index_map['lam_bx_in_lam_w'][stage]] = ubx_lam - lbx_lam
        self.lam_g0[self.index_map['lam_g_in_lam_g'][stage]] = lg_lam_hard - ug_lam_hard
        self.lam_g0[self.index_map['lam_g_sl_in_lam_g'][stage]] = -lg_lam_soft
        self.lam_g0[self.index_map['lam_g_su_in_lam_g'][stage]] = ug_lam_soft
        # TODO: separate soft contributions into bu, bx, g for more accurate warm-starting
        self.lam_x0[self.index_map['lam_sl_in_lam_w'][stage]] = -np.concatenate((slbu_lam, slbx_lam, slg_lam))
        self.lam_x0[self.index_map['lam_su_in_lam_w'][stage]] = -np.concatenate((subu_lam, subx_lam, subg_lam))

    def get_cost(self) -> float:
        """Optimal objective value of the last solve."""
        if self.qp_sol is None:
            raise RuntimeError("No solution available; call solve() first.")
        return float(self.qp_sol['f'])

    def get_stats(self, field: str) -> Union[int, float, None]:
        """
        Return a solver statistic.

        :param field: ``'nlp_iter'`` or ``'time_tot'``.
        """
        if field == 'nlp_iter':
            return self._nlp_iter
        elif field == 'time_tot':
            return self._time_total
        else:
            raise NotImplementedError(f"get_stats(): field '{field}' not implemented.")

    def get_iterate(self) -> AcadosOcpIterate:
        """Return the last solution as an :class:`AcadosOcpIterate`."""
        d = {}
        for field in ['x', 'u', 'z', 'sl', 'su', 'pi', 'lam']:
            traj = []
            for n in range(self.acados_ocp_qp.N + 1):
                if n < self.acados_ocp_qp.N or field not in ['u', 'pi', 'z']:
                    traj.append(self.get(n, field))
            d[f'{field}_traj'] = traj
        return AcadosOcpIterate(**d)

    def set_iterate(self, iterate: Union[AcadosOcpIterate, AcadosOcpFlattenedIterate]):
        """Load an iterate for warm-starting."""
        is_flat = isinstance(iterate, AcadosOcpFlattenedIterate)
        for key, traj in iterate.__dict__.items():
            field = key.replace('_traj', '')
            if field in ['x', 'u', 'pi', 'lam', 'sl', 'su']:
                if is_flat:
                    self.set_flat(field, getattr(iterate, field))
                else:
                    for n, val in enumerate(traj):
                        self.set(n, field, val)

    def get_flat(self, field: str) -> np.ndarray:
        """
        Return the concatenation of stage-wise values over all stages.

        :param field: ``'x'``, ``'u'``, ``'pi'``, ``'sl'``, ``'su'``, ``'lam'``,
                      or ``'z'``.
        """
        if self.qp_sol is None:
            raise RuntimeError("No solution available; call solve() first.")
        N = self.acados_ocp_qp.N

        if field in ['x', 'lam', 'sl', 'su']:
            return np.concatenate([self.get(i, field) for i in range(N + 1)])
        elif field in ['u', 'pi', 'z']:
            return np.concatenate([self.get(i, field) for i in range(N)])
        elif field == 'lam_x':
            return self.qp_sol_lam_w.copy()
        elif field == 'lam_g':
            return self.qp_sol_lam_g.copy()
        else:
            raise NotImplementedError(f"get_flat(): field '{field}' not implemented.")

    def set_flat(self, field: str, value: np.ndarray):
        """Set a concatenated stage-wise warm-start quantity."""
        v = np.asarray(value).flatten()
        N = self.acados_ocp_qp.N
        dims = self.acados_ocp_qp.dims

        if field == 'x':
            offset = 0
            for i in range(N + 1):
                nx_i = dims.nx[i]
                self.set(i, 'x', v[offset:offset + nx_i]); offset += nx_i
        elif field == 'u':
            offset = 0
            for i in range(N):
                nu_i = dims.nu[i]
                self.set(i, 'u', v[offset:offset + nu_i]); offset += nu_i
        elif field == 'pi':
            offset = 0
            for i in range(N):
                nx_next = dims.nx[i + 1]
                self.set(i, 'pi', v[offset:offset + nx_next]); offset += nx_next
        elif field in ['sl', 'su']:
            offset = 0
            for i in range(N + 1):
                ns_i = dims.ns[i]
                if ns_i > 0:
                    self.set(i, field, v[offset:offset + ns_i]); offset += ns_i
        elif field == 'lam':
            offset = 0
            for i in range(N + 1):
                nbu = dims.nbu[i]
                nbx = dims.nbx[i]
                ng_hard = len(self.index_map['lam_g_in_lam_g'][i])
                ns_i = dims.ns[i]
                n_lam_i = 2 * (nbu + nbx + ng_hard + ns_i)
                self.set(i, 'lam', v[offset:offset + n_lam_i]); offset += n_lam_i
        else:
            raise NotImplementedError(f"set_flat(): field '{field}' not implemented.")