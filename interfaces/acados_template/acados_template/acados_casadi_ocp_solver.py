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

import casadi as ca

from typing import Union, Optional, List

import numpy as np

from .utils import casadi_length
from .acados_ocp import AcadosOcp
from .acados_ocp_iterate import AcadosOcpIterate, AcadosOcpFlattenedIterate
from .acados_casadi_ocp import AcadosCasadiOcp

class AcadosCasadiOcpSolver:

    def __init__(self, ocp: AcadosOcp, solver: str = "ipopt", verbose=True,
                 casadi_solver_opts: Optional[dict] = None,
                 use_acados_hessian: bool = False,
                 use_single_shooting: bool = False):

        if not isinstance(ocp, AcadosOcp):
            raise TypeError('ocp should be of type AcadosOcp.')

        self.ocp = ocp
        self.multiple_shooting = not use_single_shooting
        # create casadi NLP formulation
        casadi_nlp_obj = AcadosCasadiOcp(ocp = ocp,
                                         with_hessian = use_acados_hessian,
                                         multiple_shooting= self.multiple_shooting
                                        )

        self.acados_casadi_ocp = casadi_nlp_obj

        self.casadi_nlp = casadi_nlp_obj.nlp
        self.bounds = casadi_nlp_obj.bounds
        self.w0 = casadi_nlp_obj.w0
        self.p = casadi_nlp_obj.p_nlp_values
        self.index_map = casadi_nlp_obj.index_map
        self.nlp_hess_l_custom = casadi_nlp_obj.nlp_hess_l_custom
        if use_single_shooting:
            self.x_traj_fun = casadi_nlp_obj._x_traj_fun

        # create NLP solver
        if casadi_solver_opts is None:
            casadi_solver_opts = {}

        if solver == "fatrop":
            pi_in_lam_g_flat = [idx for sublist in self.index_map['pi_in_lam_g'] for idx in sublist]
            is_equality_array = [True if i in pi_in_lam_g_flat else False for i in range(casadi_length(self.casadi_nlp['g']))]
            casadi_solver_opts['structure_detection'] = 'auto'
            casadi_solver_opts['equality'] = is_equality_array

        if use_acados_hessian:
            casadi_solver_opts["cache"] = {"nlp_hess_l": self.nlp_hess_l_custom}
        self.casadi_solver = ca.nlpsol("nlp_solver", solver, self.casadi_nlp, casadi_solver_opts)

        # create solution and initial guess
        self.lam_x0 = np.zeros(self.casadi_nlp['x'].shape).flatten()
        self.lam_g0 = np.zeros(self.casadi_nlp['g'].shape).flatten()
        self.nlp_sol = None

    def solve_for_x0(self, x0_bar):
        """
        Wrapper around `solve()` which sets initial state constraint, solves the OCP, and returns u0.
        """
        self.set(0, 'lbx', x0_bar)
        self.set(0, 'ubx', x0_bar)

        status = self.solve()

        u0 = self.get(0, "u")
        return u0

    def solve(self) -> int:
        """
        Solve the ocp with current input.

        :return: status of the solver
        """
        self.nlp_sol = self.casadi_solver(x0=self.w0, p=self.p,
                                          lam_g0=self.lam_g0, lam_x0=self.lam_x0,
                                          lbx=self.bounds['lbx'], ubx=self.bounds['ubx'],
                                          lbg=self.bounds['lbg'], ubg=self.bounds['ubg']
                                          )
        self.nlp_sol_w = self.nlp_sol['x'].full()
        self.nlp_sol_g = self.nlp_sol['g'].full()
        self.nlp_sol_lam_g = self.nlp_sol['lam_g'].full()
        self.nlp_sol_lam_w = self.nlp_sol['lam_x'].full()
        if not self.multiple_shooting:
            self.nlp_sol_x = [self.x_traj_fun[0]]
            for i in range(0, self.ocp.dims.N):
                x_current = self.nlp_sol_x[i]
                if self.ocp.solver_options.integrator_type == "DISCRETE":
                    x_next = self.x_traj_fun[i+1](x_current, self.nlp_sol_w[i], self.ocp.parameter_values, self.ocp.p_global_values)
                elif self.ocp.solver_options.integrator_type == "ERK":
                    param = np.concatenate([self.nlp_sol_w[i], self.ocp.parameter_values, self.ocp.p_global_values])
                    x_next = self.x_traj_fun[i+1](x_current, param, self.ocp.solver_options.time_steps[i])
                self.nlp_sol_x.append(x_next.full())

        # statistics
        solver_stats = self.casadi_solver.stats()
        # timing = solver_stats['t_proc_total']
        self.status = solver_stats['return_status'] if 'return_status' in solver_stats else solver_stats['success']
        self.nlp_iter = solver_stats['iter_count'] if 'iter_count' in solver_stats else None
        self.time_total = solver_stats['t_wall_total'] if 't_wall_total' in solver_stats else None
        self.solver_stats = solver_stats
        # nlp_res = ca.norm_inf(sol['g']).full()[0][0]
        # cost_val = ca.norm_inf(sol['f']).full()[0][0]
        return self.status

    def get_dim_flat(self, field: str):
        """
        Get dimension of flattened iterate.
        """
        if field not in ['x', 'u', 'z', 'pi', 'lam', 'sl', 'su', 'p']:
            raise ValueError(f'AcadosOcpSolver.get_dim_flat(field={field}): \'{field}\' is an invalid argument.')

        raise NotImplementedError()

    def get(self, stage: int, field: str):
        """
        Get the last solution of the solver.

        :param stage: integer corresponding to shooting node
        :param field: string in ['x', 'u', 'pi', 'p', 'lam', 'sl', 'su']

        .. note:: regarding lam: \n
        the inequalities are internally organized in the following order: \n
        [ lbu lbx lg lh lphi ubu ubx ug uh uphi; \n
        lsbu lsbx lsg lsh lsphi usbu usbx usg ush usphi]

        """
        if not isinstance(stage, int):
            raise TypeError('stage should be integer.')
        if self.nlp_sol is None:
            raise ValueError('No solution available. Please call solve() first.')
        dims = self.ocp.dims
        if field == 'x' and self.multiple_shooting:
            return self.nlp_sol_w[self.index_map['x_in_w'][stage]].flatten()
        elif field == 'x' and not self.multiple_shooting:
            return self.nlp_sol_x[stage].flatten()
        elif field == 'u':
            return self.nlp_sol_w[self.index_map['u_in_w'][stage]].flatten()
        elif field == 'pi' and self.multiple_shooting:
            return -self.nlp_sol_lam_g[self.index_map['pi_in_lam_g'][stage]].flatten()
        elif field == 'pi' and not self.multiple_shooting:
            return []
        elif field == 'p':
            return self.p[self.index_map['p_in_p_nlp'][stage]].flatten()
        elif field == 'sl':
            return self.nlp_sol_w[self.index_map['sl_in_w'][stage]].flatten()
        elif field == 'su':
            return self.nlp_sol_w[self.index_map['su_in_w'][stage]].flatten()
        elif field == 'lam':
            if stage == 0:
                bx_lam = self.nlp_sol_lam_w[self.index_map['lam_bx_in_lam_w'][stage]] if self.multiple_shooting else []
                bu_lam = self.nlp_sol_lam_w[self.index_map['lam_bu_in_lam_w'][stage]]
                g_lam = self.nlp_sol_lam_g[self.index_map['lam_gnl_in_lam_g'][stage]]
            elif stage < dims.N:
                bx_lam = self.nlp_sol_lam_w[self.index_map['lam_bx_in_lam_w'][stage]] if self.multiple_shooting else []
                bu_lam = self.nlp_sol_lam_w[self.index_map['lam_bu_in_lam_w'][stage]]
                g_lam = self.nlp_sol_lam_g[self.index_map['lam_gnl_in_lam_g'][stage]]
            elif stage == dims.N:
                bx_lam = self.nlp_sol_lam_w[self.index_map['lam_bx_in_lam_w'][stage]] if self.multiple_shooting else []
                bu_lam = np.empty((0, 1))
                g_lam = self.nlp_sol_lam_g[self.index_map['lam_gnl_in_lam_g'][stage]]

            lbx_lam = np.maximum(0, -bx_lam) if self.multiple_shooting else np.empty((0, 1))
            ubx_lam = np.maximum(0, bx_lam) if self.multiple_shooting else np.empty((0, 1))
            lbu_lam = np.maximum(0, -bu_lam)
            ubu_lam = np.maximum(0, bu_lam)
            if any([dims.ns_0, dims.ns, dims.ns_e]):
                lw_soft_lam = self.nlp_sol_lam_w[self.index_map['sl_in_w'][stage]]
                uw_soft_lam = self.nlp_sol_lam_w[self.index_map['su_in_w'][stage]]
                lg_soft_lam = self.nlp_sol_lam_g[self.index_map['lam_sl_in_lam_g'][stage]]
                ug_soft_lam = self.nlp_sol_lam_g[self.index_map['lam_su_in_lam_g'][stage]]
                if self.index_map['lam_su_in_lam_g'][stage]:
                    g_indices = np.array(self.index_map['lam_gnl_in_lam_g'][stage]+\
                                        self.index_map['lam_sl_in_lam_g'][stage])
                    sorted_indices = np.argsort(g_indices)
                    g_lam_lower = np.concatenate((np.maximum(0, -g_lam), -lg_soft_lam))
                    lbg_lam = g_lam_lower[sorted_indices]
                    g_lam_upper = np.concatenate((np.maximum(0, g_lam), ug_soft_lam))
                    ubg_lam = g_lam_upper[sorted_indices]
                else:
                    lbg_lam = np.abs(lg_soft_lam)
                    ubg_lam = np.abs(ug_soft_lam)
                lam_soft = np.concatenate((-lw_soft_lam, -uw_soft_lam))
            else:
                lbg_lam = np.maximum(0, -g_lam)
                ubg_lam = np.maximum(0, g_lam)
                lam_soft = np.empty((0, 1))
            lam = np.concatenate((lbu_lam, lbx_lam, lbg_lam, ubu_lam, ubx_lam, ubg_lam, lam_soft))
            return lam.flatten()
        elif field in ['z']:
            return np.empty((0,))  # Only empty is supported for now. TODO: extend.
        else:
            raise NotImplementedError(f"Field '{field}' is not implemented in AcadosCasadiOcpSolver")

    def get_flat(self, field_: str) -> np.ndarray:
        """
        Get concatenation of all stages of last solution of the solver.

        :param field: string in ['x', 'u', 'pi', 'lam', 'p', 'p_global', 'z', 'sl', 'su']

        .. note:: The parameter 'p_global' has no stage-wise structure and is processed in a memory saving manner by default. \n
                In order to read the 'p_global' parameter, the option 'save_p_global' must be set to 'True' upon instantiation. \n
        """
        if self.nlp_sol is None:
            raise ValueError('No solution available. Please call solve() first.')
        dims = self.ocp.dims
        result = []

        if field_ in ['x', 'lam', 'sl', 'su', 'p']:
            for i in range(dims.N+1):
                result.append(self.get(i, field_))
            return np.concatenate(result)
        elif field_ in ['u', 'pi', 'z']:
            for i in range(dims.N):
                result.append(self.get(i, field_))
            return np.concatenate(result)
        elif field_ == 'p_global':
            return self.p[self.index_map['p_global_in_p_nlp']].flatten()
        # casadi variables. TODO: maybe remove this.
        elif field_ == 'lam_x':
            return self.nlp_sol_lam_w.flatten()
        elif field_ == 'lam_g':
            return self.nlp_sol_lam_g.flatten()
        elif field_ == 'lam_p':
            return self.nlp_sol['lam_p'].full().flatten()
        else:
            raise NotImplementedError(f"Field '{field_}' is not implemented in get_flat().")

    def set_flat(self, field_: str, value_: np.ndarray) -> None:
        """
        Set concatenation solver initialization.

        :param field: string in ['x', 'u', 'lam', pi]
        """
        dims = self.ocp.dims
        if field_ == 'x':
            for i in range(dims.N+1):
                self.set(i, 'x', value_[i*dims.nx:(i+1)*dims.nx])
        elif field_ == 'u':
            for i in range(dims.N):
                self.set(i, 'u', value_[i*dims.nu:(i+1)*dims.nu])
        elif field_ == 'pi':
            for i in range(dims.N):
                self.set(i, 'pi', value_[i*dims.nx:(i+1)*dims.nx])
        elif field_ == 'lam':
            offset = 0
            for i in range(dims.N+1):
                if i == 0:
                    n_lam_i = 2 * (dims.nbx_0 + dims.nbu + dims.ng + dims.nh_0 + dims.nphi_0)
                elif i < dims.N:
                    n_lam_i = 2 * (dims.nbx + dims.nbu + dims.ng + dims.nh + dims.nphi)
                elif i == dims.N:
                    n_lam_i = 2 * (dims.nbx_e + dims.ng_e + dims.nh_e + dims.nphi_e)
                self.set(i, 'lam', value_[offset : offset + n_lam_i])
                offset += n_lam_i
        elif field_ in ['sl', 'su']:
            offset = 0
            for i in range(dims.N+1):
                if i == 0:
                    ns_i = dims.ns_0
                elif i < dims.N:
                    ns_i = dims.ns
                elif i == dims.N:
                    ns_i = dims.ns_e
                if ns_i > 0:
                    self.set(i, field_, value_[offset : offset + ns_i])
                    offset += ns_i
        else:
            raise NotImplementedError(f"Field '{field_}' is not yet implemented in set_flat().")

    def load_iterate(self, filename:str, verbose: bool = True):
        raise NotImplementedError()

    def store_iterate_to_obj(self) -> AcadosOcpIterate:
        """
        Returns the current iterate of the OCP solver as an AcadosOcpIterate.
        """
        d = {}
        for field in ["x", "u", "z", "sl", "su", "pi", "lam"]:
            traj = []
            for n in range(self.ocp.dims.N+1):
                if n < self.ocp.dims.N or not (field in ["u", "pi", "z"]):
                    traj.append(self.get(n, field))

            d[f"{field}_traj"] = traj

        return AcadosOcpIterate(**d)

    def load_iterate_from_obj(self, iterate: AcadosOcpIterate) -> None:
        """
        Loads the provided iterate into the OCP solver.
        Note: The iterate object does not contain the the parameters.
        """
        for key, traj in iterate.__dict__.items():
            field = key.replace('_traj', '')

            for n, val in enumerate(traj):
                if field in ['x', 'u', 'pi', 'lam', 'sl', 'su']:
                    self.set(n, field, val)

    def store_iterate_to_flat_obj(self) -> AcadosOcpFlattenedIterate:
        """
        Returns the current iterate of the OCP solver as an AcadosOcpFlattenedIterate.
        """
        return AcadosOcpFlattenedIterate(x = self.get_flat("x"),
                                         u = self.get_flat("u"),
                                         pi = self.get_flat("pi"),
                                         lam = self.get_flat("lam"),
                                         sl = self.get_flat("sl"),
                                         su = self.get_flat("su"),
                                         z = self.get_flat("z"))

    def load_iterate_from_flat_obj(self, iterate: AcadosOcpFlattenedIterate) -> None:
        """
        Loads the provided iterate into the OCP solver.
        Note: The iterate object does not contain the the parameters.
        """
        self.set_flat("x", iterate.x)
        self.set_flat("u", iterate.u)
        self.set_flat("pi", iterate.pi)
        self.set_flat("lam", iterate.lam)
        self.set_flat("sl", iterate.sl)
        self.set_flat("su", iterate.su)

    def get_stats(self, field_: str) -> Union[int, float, np.ndarray]:

        if field_ == "nlp_iter":
            return self.nlp_iter
        elif field_ == "time_tot":
            return self.time_total
        else:
            raise NotImplementedError()

    def get_cost(self) -> float:
        return self.nlp_sol['f'].full().item()

    def set(self, stage: int, field: str, value_: np.ndarray):
        """
        Set solver initialization to stages.

        :param stage: integer corresponding to shooting node
        :param field: string in ['x', 'u', 'pi', 'lam', 'p', 'sl', 'su']
        :value_:
        """
        dims = self.ocp.dims

        if field == 'x' and self.multiple_shooting:
            self.w0[self.index_map['x_in_w'][stage]] = value_.flatten()
        elif field == 'x' and not self.multiple_shooting:
            pass
        elif field == 'u':
            self.w0[self.index_map['u_in_w'][stage]] = value_.flatten()
        elif field == 'pi' and self.multiple_shooting:
            self.lam_g0[self.index_map['pi_in_lam_g'][stage]] = -value_.flatten()
        elif field == 'pi' and not self.multiple_shooting:
            pass
        elif field == 'p':
            self.p[self.index_map['p_in_p_nlp'][stage]] = value_.flatten()
        elif field == 'sl':
            self.w0[self.index_map['sl_in_w'][stage]] = value_.flatten()
        elif field == 'su':
            self.w0[self.index_map['su_in_w'][stage]] = value_.flatten()
        elif field == 'lam':
            if stage == 0:
                nbx = dims.nbx_0 if self.multiple_shooting else 0
                nbu = dims.nbu
                n_ghphi = dims.ng + dims.nh_0 + dims.nphi_0
                ns = dims.ns_0
            elif stage < dims.N:
                nbx = dims.nbx if self.multiple_shooting else 0
                nbu = dims.nbu
                n_ghphi = dims.ng + dims.nh + dims.nphi
                ns = dims.ns
            elif stage == dims.N:
                nbx = dims.nbx_e if self.multiple_shooting else 0
                nbu = 0
                n_ghphi = dims.ng_e + dims.nh_e + dims.nphi_e
                ns = dims.ns_e

            offset_u = (nbx+nbu+n_ghphi)
            lbu_lam = value_[:nbu]
            lbx_lam = value_[nbu:nbu+nbx]
            lg_lam = value_[nbu+nbx:nbu+nbx+n_ghphi]
            ubu_lam = value_[offset_u:offset_u+nbu]
            ubx_lam = value_[offset_u+nbu:offset_u+nbu+nbx]
            ug_lam = value_[offset_u+nbu+nbx:offset_u+nbu+nbx+n_ghphi]
            offset_soft = 2*offset_u
            soft_lam = value_[offset_soft:offset_soft + 2 * ns]

            g_indices = np.array(self.index_map['lam_gnl_in_lam_g'][stage]+\
                                self.index_map['lam_sl_in_lam_g'][stage])
            sorted = np.sort(g_indices)
            gnl_indices = [i for i, x in enumerate(sorted) if x in self.index_map['lam_gnl_in_lam_g'][stage]]
            sl_indices = [i for i, x in enumerate(sorted) if x in self.index_map['lam_sl_in_lam_g'][stage]]
            lg_lam_hard = lg_lam[gnl_indices]
            lg_lam_soft = lg_lam[sl_indices]
            ug_lam_hard = ug_lam[gnl_indices]
            ug_lam_soft = ug_lam[sl_indices]

            if stage != dims.N:
                if self.multiple_shooting:
                    self.lam_x0[self.index_map['lam_bx_in_lam_w'][stage]+self.index_map['lam_bu_in_lam_w'][stage]] = np.concatenate((ubx_lam-lbx_lam, ubu_lam-lbu_lam))
                else:
                    self.lam_x0[self.index_map['lam_bu_in_lam_w'][stage]] = ubu_lam-lbu_lam
                self.lam_g0[self.index_map['lam_gnl_in_lam_g'][stage]] =  ug_lam_hard-lg_lam_hard
                self.lam_g0[self.index_map['lam_sl_in_lam_g'][stage]] = -lg_lam_soft
                self.lam_g0[self.index_map['lam_su_in_lam_g'][stage]] = ug_lam_soft
                self.lam_x0[self.index_map['sl_in_w'][stage]+self.index_map['su_in_w'][stage]] = -soft_lam
            else:
                if self.multiple_shooting:
                    self.lam_x0[self.index_map['lam_bx_in_lam_w'][stage]] = ubx_lam-lbx_lam
                self.lam_g0[self.index_map['lam_gnl_in_lam_g'][stage]] = ug_lam_hard-lg_lam_hard
                self.lam_g0[self.index_map['lam_sl_in_lam_g'][stage]] = -lg_lam_soft
                self.lam_g0[self.index_map['lam_su_in_lam_g'][stage]] = ug_lam_soft
                self.lam_x0[self.index_map['sl_in_w'][stage]+self.index_map['su_in_w'][stage]] = -soft_lam
        elif field == 'lbx':
            self.bounds['lbx'][self.index_map['lam_bx_in_lam_w'][stage]] = value_.flatten()
        elif field == 'ubx':
            self.bounds['ubx'][self.index_map['lam_bx_in_lam_w'][stage]] = value_.flatten()
        elif field == 'yref':
            self.p[self.index_map['yref_in_p_nlp'][stage]] = value_.flatten()
        else:
            raise NotImplementedError(f"Field '{field}' is not yet implemented in set().")

    def set_params_sparse(self, stage_: int, idx_values_: np.ndarray, param_values_: np.ndarray):
        if not isinstance(stage_, int):
            raise TypeError('stage should be integer.')
        index = np.asarray(self.index_map['p_in_p_nlp'][stage_])[idx_values_]
        self.p[index] = param_values_.flatten()

    def cost_get(self, stage_: int, field_: str) -> np.ndarray:
        raise NotImplementedError()

    def cost_set(self, stage_: int, field_: str, value_):
        raise NotImplementedError()

    def get_constraints_value(self, stage: int):
        """
        Get the constraints values and lambda for a given stage.
        """
        if not isinstance(stage, int):
            raise TypeError('stage should be integer.')
        if self.nlp_sol is None:
            raise ValueError('No solution available. Please call solve() first.')

        # create constraints value and lambda in the order as [bx, bu, bg, bh, bphi]
        if stage < self.ocp.dims.N:
            constraints_value = np.concatenate((self.nlp_sol_w[self.index_map['lam_bx_in_lam_w'][stage]],
                                                self.nlp_sol_w[self.index_map['lam_bu_in_lam_w'][stage]],
                                                self.nlp_sol_g[self.index_map['pi_in_lam_g'][stage]],
                                                self.nlp_sol_g[self.index_map['lam_gnl_in_lam_g'][stage]])).flatten()
            lambda_values = np.concatenate((self.nlp_sol_lam_w[self.index_map['lam_bx_in_lam_w'][stage]],
                                            self.nlp_sol_lam_w[self.index_map['lam_bu_in_lam_w'][stage]],
                                            self.nlp_sol_lam_g[self.index_map['pi_in_lam_g'][stage]],
                                            self.nlp_sol_lam_g[self.index_map['lam_gnl_in_lam_g'][stage]])).flatten()
            lb = ca.vertcat(self.bounds['lbx'][self.index_map['lam_bx_in_lam_w'][stage]],
                            self.bounds['lbx'][self.index_map['lam_bu_in_lam_w'][stage]],
                            self.bounds['lbg'][self.index_map['pi_in_lam_g'][stage]],
                            self.bounds['lbg'][self.index_map['lam_gnl_in_lam_g'][stage]]).full().flatten()
            ub = ca.vertcat(self.bounds['ubx'][self.index_map['lam_bx_in_lam_w'][stage]],
                            self.bounds['ubx'][self.index_map['lam_bu_in_lam_w'][stage]],
                            self.bounds['ubg'][self.index_map['pi_in_lam_g'][stage]],
                            self.bounds['ubg'][self.index_map['lam_gnl_in_lam_g'][stage]]).full().flatten()
        elif stage == self.ocp.dims.N:
            constraints_value = np.concatenate((self.nlp_sol_w[self.index_map['lam_bx_in_lam_w'][stage]],
                                                self.nlp_sol_g[self.index_map['lam_gnl_in_lam_g'][stage]])).flatten()
            lambda_values = np.concatenate((self.nlp_sol_lam_w[self.index_map['lam_bx_in_lam_w'][stage]],
                                            self.nlp_sol_lam_g[self.index_map['lam_gnl_in_lam_g'][stage]])).flatten()
            lb = ca.vertcat(self.bounds['lbx'][self.index_map['lam_bx_in_lam_w'][stage]],
                            self.bounds['lbg'][self.index_map['lam_gnl_in_lam_g'][stage]]).full().flatten()
            ub = ca.vertcat(self.bounds['ubx'][self.index_map['lam_bx_in_lam_w'][stage]],
                            self.bounds['ubg'][self.index_map['lam_gnl_in_lam_g'][stage]]).full().flatten()
        return  constraints_value, lambda_values, lb, ub

    def get_constraints_indices(self, stage: int):
        """
        Get the indices of the constraints for a given stage.
        This function distinguishes between inequality and equality constraints
        returns indices of
        (inequality, equality for decision variables, equality for dynamic and gnl, lower active inequality, upper active inequality).
        """
        constraints_value, _, lb, ub = self.get_constraints_value(stage)
        tol = self.ocp.solver_options.nlp_solver_tol_ineq
        # distinguish between equality and inequality constraints
        if stage == 0:
            nbx = self.ocp.dims.nbx_0
            nbu = self.ocp.dims.nbu
        elif stage < self.ocp.dims.N:
            nbx = self.ocp.dims.nbx
            nbu = self.ocp.dims.nbu
        elif stage == self.ocp.dims.N:
            nbx = self.ocp.dims.nbx_e
            nbu = 0

        ineq_indices = []
        eq_indices_bounds = []
        eq_indices_ca_g = []

        for i in range(len(lb)):
            if lb[i] != ub[i]:
                ineq_indices.append(i)
            else:
                #distinguish between equality in decision variables and in constraints
                if i in range(nbx + nbu):
                    eq_indices_bounds.append(i)
                else:
                    eq_indices_ca_g.append(i)
        # get the inequality violations
        violations_ineq_lb = constraints_value[ineq_indices] - lb[ineq_indices]
        violations_ineq_ub = ub[ineq_indices] - constraints_value[ineq_indices]
        # any negative value in violations means infeasible constraint, raise an error
        if np.any(violations_ineq_lb < -tol) or np.any(violations_ineq_ub < -tol):
            raise ValueError('Constraints are violated. Please check the solution.')
        # get active inequality indices from inequality constraints
        active_ineq_lb_indices = np.take(ineq_indices, np.where(violations_ineq_lb < tol)[0])
        active_ineq_ub_indices = np.take(ineq_indices, np.where(violations_ineq_ub < tol)[0])
        return ineq_indices, eq_indices_bounds, eq_indices_ca_g, active_ineq_lb_indices, active_ineq_ub_indices

    def satisfies_strict_complementarity_stage_wise(self, stage: int, tol: float) -> bool:
        """
        Check if the solution satisfies strict complementarity conditions for a given stage.
        This checks that the Lagrange multipliers for active inequality constraints are strictly positive.
        Not tested yet.
        """
        tol = self.ocp.solver_options.nlp_solver_tol_ineq
        if self.nlp_sol is None:
            raise ValueError('No solution available. Please call solve() first.')

        _, lambda_value, _, _ = self.get_constraints_value(stage)
        _, _, _, active_ineq_lb_indices, active_ineq_ub_indices = self.get_constraints_indices(stage)

        for i in active_ineq_lb_indices:
            lam = np.maximum(0, -lambda_value[i])
            if lam < tol:
                return False
        for i in active_ineq_ub_indices:
            lam = np.maximum(0, lambda_value[i])
            if lam < tol:
                return False
        return True

    def satisfies_strict_complementarity_stages(self, tol: float) -> List[bool]:
        """
        Check if the solution satisfies strict complementarity conditions for all stages.
        Not tested yet.
        """
        tol = self.ocp.solver_options.nlp_solver_tol_ineq
        dims = self.ocp.dims
        complementarity = []
        for stage in range(dims.N + 1):
            complementarity.append(self.satisfies_strict_complementarity_stage_wise(stage, tol))
        return complementarity

    def satisfies_strict_complementarity(self) -> bool:
        """
        Check if the solution satisfies strict complementarity conditions for all stages.
        Not tested yet.
        """
        stage_wise_complementarity = self.satisfies_strict_complementarity_stages(self.ocp.solver_options.nlp_solver_tol_ineq)
        if all(stage_wise_complementarity):
            return True
        else:
            return False

    def satisfies_LICQ_stage_wise(self, stage) -> bool:
        """
        Check if the solution satisfies the Linear Independence Constraint Qualification (LICQ) for a given stage.
        """
        if self.nlp_sol is None:
            raise ValueError('No solution available. Please call solve() first.')

        _, eq_indices_bounds, eq_indices_ca_g, active_ineq_lb_indices, active_ineq_ub_indices = self.get_constraints_indices(stage)

        w, w_value, constraints_expr_stage, eq_indices = self._get_w_and_constraints_for_LICQ(stage, eq_indices_bounds, eq_indices_ca_g)

        eq_constraints = constraints_expr_stage[eq_indices] if len(eq_indices) != 0 else ca.vertcat()
        active_ineq_lb_constraints = constraints_expr_stage[active_ineq_lb_indices] if len(active_ineq_lb_indices) != 0 else ca.vertcat()
        active_ineq_ub_constraints = constraints_expr_stage[active_ineq_ub_indices] if len(active_ineq_ub_indices) != 0 else ca.vertcat()
        constraint_matrix = ca.vertcat(eq_constraints,
                                       -active_ineq_lb_constraints,
                                       active_ineq_ub_constraints)
        constraint_jac_expr = ca.Function('constraint_jac', [w], [ca.jacobian(constraint_matrix, w)])
        constraint_jac_value = constraint_jac_expr(w_value).full()

        # Check if the Jacobian of the constraints is full rank
        rank = np.linalg.matrix_rank(constraint_jac_value) if constraint_jac_value.any() else 0
        if rank == constraint_jac_value.shape[0]:
            return True
        else:
            return False

    def satisfies_LICQ_stages(self) -> List[bool]:
        """
        Check if the solution satisfies the Linear Independence Constraint Qualification (LICQ) for all stages.
        return a list of booleans, each indicating whether LICQ is satisfied for the corresponding stage.
        """
        dims = self.ocp.dims
        stage_wise_LICQ = []
        for stage in range(dims.N + 1):
            stage_wise_LICQ.append(self.satisfies_LICQ_stage_wise(stage))
        return stage_wise_LICQ

    def satisfies_LICQ(self) -> bool:
        """
        Check if the solution satisfies the Linear Independence Constraint Qualification (LICQ) for all stages.
        return True if LICQ is satisfied for all stages, otherwise False.
        """
        stage_wise_LICQ = self.satisfies_LICQ_stages()
        if all(stage_wise_LICQ):
            return True
        else:
            return False

    def _get_w_and_constraints_for_LICQ(self, stage: int, eq_indices_bounds, eq_indices_ca_g):
        """
        Helper function to get the w vector and constraints expression for a given stage.
        This is used to compute the Jacobian of the constraints for checking LICQ.
        """
        if stage == 0:
            w = ca.vertcat(self.casadi_nlp['x'][self.index_map['x_in_w'][stage]],
                           self.casadi_nlp['x'][self.index_map['u_in_w'][stage]])
            w_value = np.concatenate((self.nlp_sol_w[self.index_map['x_in_w'][stage]],
                                      self.nlp_sol_w[self.index_map['u_in_w'][stage]])).flatten()
            constraints_expr_stage = ca.vertcat(self.casadi_nlp['x'][self.index_map['lam_bx_in_lam_w'][stage]],
                                                self.casadi_nlp['x'][self.index_map['lam_bu_in_lam_w'][stage]],
                                                self.casadi_nlp['g'][self.index_map['pi_in_lam_g'][stage]],
                                                self.casadi_nlp['g'][self.index_map['lam_gnl_in_lam_g'][stage]])
            eq_indices = eq_indices_ca_g
        elif stage < self.ocp.dims.N:
            w = ca.vertcat(self.casadi_nlp['x'][self.index_map['x_in_w'][stage]],
                        self.casadi_nlp['x'][self.index_map['u_in_w'][stage]])
            w_value = np.concatenate((self.nlp_sol_w[self.index_map['x_in_w'][stage]],
                                    self.nlp_sol_w[self.index_map['u_in_w'][stage]])).flatten()
            constraints_expr_stage = ca.vertcat(self.casadi_nlp['x'][self.index_map['lam_bx_in_lam_w'][stage]],
                                                self.casadi_nlp['x'][self.index_map['lam_bu_in_lam_w'][stage]],
                                                self.casadi_nlp['g'][self.index_map['pi_in_lam_g'][stage]],
                                                self.casadi_nlp['g'][self.index_map['lam_gnl_in_lam_g'][stage]])
            eq_indices = eq_indices_bounds + eq_indices_ca_g
        elif stage == self.ocp.dims.N:
            w = self.casadi_nlp['x'][self.index_map['x_in_w'][stage]]
            w_value = self.nlp_sol_w[self.index_map['x_in_w'][stage]].flatten()
            constraints_expr_stage = ca.vertcat(self.casadi_nlp['x'][self.index_map['lam_bx_in_lam_w'][stage]],
                                                self.casadi_nlp['g'][self.index_map['lam_gnl_in_lam_g'][stage]])
            eq_indices = eq_indices_bounds + eq_indices_ca_g
        return w, w_value, constraints_expr_stage, eq_indices