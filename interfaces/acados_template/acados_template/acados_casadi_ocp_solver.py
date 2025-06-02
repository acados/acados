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

import casadi as ca

from typing import Union, Tuple

import numpy as np
from .acados_ocp import AcadosOcp
from .acados_ocp_iterate import AcadosOcpIterate, AcadosOcpIterates, AcadosOcpFlattenedIterate

class AcadosCasadiOcpSolver:

    @classmethod
    def create_casadi_nlp_formulation(cls, ocp: AcadosOcp) -> Tuple[dict, dict, np.ndarray]:
        """
        Creates an equivalent CasADi NLP formulation of the OCP.
        Experimental, not fully implemented yet.

        :return: nlp_dict, bounds_dict, w0 (initial guess)
        """
        ocp.make_consistent()

        # unpack
        model = ocp.model
        dims = ocp.dims
        constraints = ocp.constraints
        solver_options = ocp.solver_options
        index_map = {'x': [], 'u': [], 
                    'lam_g_dynamic': [], 'lam_g_constraint': []}

        if any([dims.ns_0, dims.ns, dims.ns_e]):
            raise NotImplementedError("CasADi NLP formulation not implemented for formulations with soft constraints yet.")

        # create variables indexed by shooting nodes, vectorize in the end
        ca_symbol = model.get_casadi_symbol()
        xtraj_node = [ca_symbol(f'x{i}', dims.nx, 1) for i in range(solver_options.N_horizon+1)]
        utraj_node = [ca_symbol(f'u{i}', dims.nu, 1) for i in range(solver_options.N_horizon)]
        if dims.nz > 0:
            raise NotImplementedError("CasADi NLP formulation not implemented for models with algebraic variables (z).")
        # parameters
        ptraj_node = [ca_symbol(f'p{i}', dims.np, 1) for i in range(solver_options.N_horizon)]

        ### Constraints: bounds
        # setup state bounds
        lb_xtraj_node = [-np.inf * ca.DM.ones((dims.nx, 1)) for _ in range(solver_options.N_horizon+1)]
        ub_xtraj_node = [np.inf * ca.DM.ones((dims.nx, 1)) for _ in range(solver_options.N_horizon+1)]
        lb_xtraj_node[0][constraints.idxbx_0] = constraints.lbx_0
        ub_xtraj_node[0][constraints.idxbx_0] = constraints.ubx_0
        offset = 0
        for i in range(1, solver_options.N_horizon):
            lb_xtraj_node[i][constraints.idxbx] = constraints.lbx
            ub_xtraj_node[i][constraints.idxbx] = constraints.ubx
        lb_xtraj_node[-1][constraints.idxbx_e] = constraints.lbx_e
        ub_xtraj_node[-1][constraints.idxbx_e] = constraints.ubx_e

        # setup control bounds
        lb_utraj_node = [-np.inf * ca.DM.ones((dims.nu, 1)) for _ in range(solver_options.N_horizon)]
        ub_utraj_node = [np.inf * ca.DM.ones((dims.nu, 1)) for _ in range(solver_options.N_horizon)]
        for i in range(solver_options.N_horizon):
            lb_utraj_node[i][constraints.idxbu] = constraints.lbu
            ub_utraj_node[i][constraints.idxbu] = constraints.ubu

        ### Nonlinear constraints
        g = []
        lbg = []
        ubg = []
        offset = 0
        # dynamics
        if solver_options.integrator_type == "DISCRETE":
            f_discr_fun = ca.Function('f_discr_fun', [model.x, model.u, model.p, model.p_global], [model.disc_dyn_expr])
        elif solver_options.integrator_type == "ERK":
            ca_expl_ode = ca.Function('ca_expl_ode', [model.x, model.u], [model.f_expl_expr])
                                    #   , model.p, model.p_global]
            f_discr_fun = ca.simpleRK(ca_expl_ode, solver_options.sim_method_num_steps[0], solver_options.sim_method_num_stages[0])
        else:
            raise NotImplementedError(f"Integrator type {solver_options.integrator_type} not supported.")

        for i in range(solver_options.N_horizon):
            # add dynamics constraints
            if solver_options.integrator_type == "DISCRETE":
                g.append(f_discr_fun(xtraj_node[i], utraj_node[i], ptraj_node[i], model.p_global) - xtraj_node[i+1])
            elif solver_options.integrator_type == "ERK":
                g.append(f_discr_fun(xtraj_node[i], utraj_node[i], solver_options.time_steps[i]) - xtraj_node[i+1])
            lbg.append(np.zeros((dims.nx, 1)))
            ubg.append(np.zeros((dims.nx, 1)))
            index_map['lam_g_dynamic'].append(list(range(offset, offset+dims.nx)))
            offset += dims.nx

        # nonlinear constraints -- initial stage
        h0_fun = ca.Function('h0_fun', [model.x, model.u, model.p, model.p_global], [model.con_h_expr_0])
        g.append(h0_fun(xtraj_node[0], utraj_node[0], ptraj_node[0], model.p_global))
        lbg.append(constraints.lh_0)
        ubg.append(constraints.uh_0)
        index_map['lam_g_constraint'].append(list(range(offset, offset + dims.nh_0)))
        offset += dims.nh_0

        if dims.nphi_0 > 0:
            conl_constr_expr_0 = ca.substitute(model.con_phi_expr_0, model.con_r_in_phi_0, model.con_r_expr_0)
            conl_constr_0_fun = ca.Function('conl_constr_0_fun', [model.x, model.u, model.p, model.p_global], [conl_constr_expr_0])
            g.append(conl_constr_0_fun(xtraj_node[0], utraj_node[0], ptraj_node[0], model.p_global))
            lbg.append(constraints.lphi_0)
            ubg.append(constraints.uphi_0)
            index_map['lam_g_constraint'].append(list(range(offset, offset + dims.nphi_0)))
            offset += dims.nphi_0

        # nonlinear constraints -- intermediate stages
        h_fun = ca.Function('h_fun', [model.x, model.u, model.p, model.p_global], [model.con_h_expr])

        if dims.nphi > 0:
            conl_constr_expr = ca.substitute(model.con_phi_expr, model.con_r_in_phi, model.con_r_expr)
            conl_constr_fun = ca.Function('conl_constr_fun', [model.x, model.u, model.p, model.p_global], [conl_constr_expr])

        for i in range(1, solver_options.N_horizon):
            g.append(h_fun(xtraj_node[i], utraj_node[i], ptraj_node[i], model.p_global))
            lbg.append(constraints.lh)
            ubg.append(constraints.uh)
            index_map['lam_g_constraint'].append(list(range(offset, offset + dims.nh)))
            offset += dims.nh

            if dims.nphi > 0:
                g.append(conl_constr_fun(xtraj_node[i], utraj_node[i], ptraj_node[i], model.p_global))
                lbg.append(constraints.lphi)
                ubg.append(constraints.uphi)
                index_map['lam_g_constraint'].append(list(range(offset, offset + dims.nphi)))
                offset += dims.nphi

        # nonlinear constraints -- terminal stage
        h_e_fun = ca.Function('h_e_fun', [model.x, model.p, model.p_global], [model.con_h_expr_e])

        g.append(h_e_fun(xtraj_node[-1], ptraj_node[-1], model.p_global))
        lbg.append(constraints.lh_e)
        ubg.append(constraints.uh_e)
        index_map['lam_g_constraint'].append(list(range(offset, offset + dims.nh_e)))
        offset += dims.nh_e

        if dims.nphi_e > 0:
            conl_constr_expr_e = ca.substitute(model.con_phi_expr_e, model.con_r_in_phi_e, model.con_r_expr_e)
            conl_constr_e_fun = ca.Function('conl_constr_e_fun', [model.x, model.p, model.p_global], [conl_constr_expr_e])
            g.append(conl_constr_e_fun(xtraj_node[-1], ptraj_node[-1], model.p_global))
            lbg.append(constraints.lphi_e)
            ubg.append(constraints.uphi_e)
            index_map['lam_g_constraint'].append(list(range(offset, offset + dims.nphi_e)))
            offset += dims.nphi_e

        ### Cost
        # initial cost term
        nlp_cost = 0
        cost_expr_0 = ocp.get_initial_cost_expression()
        cost_fun_0 = ca.Function('cost_fun_0', [model.x, model.u, model.p, model.p_global], [cost_expr_0])
        nlp_cost += solver_options.cost_scaling[0] * cost_fun_0(xtraj_node[0], utraj_node[0], ptraj_node[0], model.p_global)

        # intermediate cost term
        cost_expr = ocp.get_path_cost_expression()
        cost_fun = ca.Function('cost_fun', [model.x, model.u, model.p, model.p_global], [cost_expr])
        for i in range(1, solver_options.N_horizon):
            nlp_cost += solver_options.cost_scaling[i] * cost_fun(xtraj_node[i], utraj_node[i], ptraj_node[i], model.p_global)

        # terminal cost term
        cost_expr_e = ocp.get_terminal_cost_expression()
        cost_fun_e = ca.Function('cost_fun_e', [model.x, model.p, model.p_global], [cost_expr_e])
        nlp_cost += solver_options.cost_scaling[-1] * cost_fun_e(xtraj_node[-1], ptraj_node[-1], model.p_global)

        ### Formulation
        # interleave primary variables w and bounds
        # w = [x0, u0, x1, u1, ...]
        w_sym_list = []
        lbw_list = []
        ubw_list = []
        w0_list = []
        offset = 0
        x_guess = ocp.constraints.x0 if ocp.constraints.has_x0 else np.zeros((dims.nx,))
        for i in range(solver_options.N_horizon):
            # add x
            w_sym_list.append(xtraj_node[i])
            lbw_list.append(lb_xtraj_node[i])
            ubw_list.append(ub_xtraj_node[i])
            w0_list.append(x_guess)
            index_map['x'].append(list(range(offset, offset + dims.nx)))
            offset += dims.nx
            # add u
            w_sym_list.append(utraj_node[i])
            lbw_list.append(lb_utraj_node[i])
            ubw_list.append(ub_utraj_node[i])
            w0_list.append(np.zeros((dims.nu,)))
            index_map['u'].append(list(range(offset, offset + dims.nu)))
            offset += dims.nu
        ## terminal stage
        # add x
        w_sym_list.append(xtraj_node[-1])
        lbw_list.append(lb_xtraj_node[-1])
        ubw_list.append(ub_xtraj_node[-1])
        w0_list.append(x_guess)
        index_map['x'].append(list(range(offset, offset + dims.nx)))
        offset += dims.nx

        # vectorize
        w = ca.vertcat(*w_sym_list)
        lbw = ca.vertcat(*lbw_list)
        ubw = ca.vertcat(*ubw_list)
        p_nlp = ca.vertcat(*ptraj_node, model.p_global)

        # create NLP
        nlp = {"x": w, "p": p_nlp, "g": ca.vertcat(*g), "f": nlp_cost}
        bounds = {"lbx": lbw, "ubx": ubw, "lbg": ca.vertcat(*lbg), "ubg": ca.vertcat(*ubg)}
        w0 = np.concatenate(w0_list)

        return nlp, bounds, w0, index_map


    def __init__(self, acados_ocp: AcadosOcp, solver: str = "ipopt", verbose=True):

        if not isinstance(acados_ocp, AcadosOcp):
            raise TypeError('acados_ocp should be of type AcadosOcp.')

        self.acados_ocp = acados_ocp
        self.casadi_nlp, self.bounds, self.w0, self.index_map = self.create_casadi_nlp_formulation(acados_ocp)
        self.casadi_solver = ca.nlpsol("nlp_solver", solver, self.casadi_nlp)
        self.lam_x0 = None
        self.lam_g0 = None
        self.lam_g0_dynamic = None
        self.lam_g0_constraint = None
        self.nlp_sol = None


    def solve_for_x0(self, x0_bar, fail_on_nonzero_status=True, print_stats_on_failure=True):
         raise NotImplementedError()


    def solve(self) -> int:
        """
        Solve the ocp with current input.

        :return: status of the solver
        """
        if  any(x is not None for x in [self.lam_x0, self.lam_g0_dynamic, self.lam_g0_constraint]):
            self.lam_g0 = np.concatenate((self.lam_g0_dynamic, self.lam_g0_constraint)).flatten()
            self.nlp_sol = self.casadi_solver(
                                            x0=self.w0,
                                            lam_x0=self.lam_x0,
                                            lam_g0=self.lam_g0,
                                            lbx=self.bounds['lbx'], ubx=self.bounds['ubx'],
                                            lbg=self.bounds['lbg'], ubg=self.bounds['ubg']
                                            )
        else:
            self.nlp_sol = self.casadi_solver(
                                            x0=self.w0,
                                            lbx=self.bounds['lbx'], ubx=self.bounds['ubx'],
                                            lbg=self.bounds['lbg'], ubg=self.bounds['ubg']
                                            )
        self.nlp_sol_w = self.nlp_sol['x'].full()
        self.nlp_sol_lam_g = self.nlp_sol['lam_g'].full()
        self.nlp_sol_lam_x = self.nlp_sol['lam_x'].full()
        # TODO: return correct status
        return 0

    # def get_dim_flat(self, field: str):
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
        :param field: string in ['x', 'u', 'pi', 'lam']

        """
        if not isinstance(stage, int):
            raise TypeError('stage should be integer.')
        if self.nlp_sol is None:
            raise ValueError('No solution available. Please call solve() first.')
        dims = self.acados_ocp.dims
        if field == 'x':
            return self.nlp_sol_w[self.index_map['x'][stage]].flatten()
        elif field == 'u':
            return self.nlp_sol_w[self.index_map['u'][stage]].flatten()
        elif field == 'pi':
            return self.nlp_sol_lam_g[self.index_map['lam_g_dynamic'][stage]].flatten()
        elif field == 'lam':
            # get the lamda of (ubx-lbx), (ubu-lbu) and (ug-lg), (uphi-lphi)
            if stage == 0:
                bx_length = self.acados_ocp.constraints.lbx_0.size
                bu_length = self.acados_ocp.constraints.lbu.size
                bx_lam = self.nlp_sol_lam_x[self.index_map['x'][stage]] if bx_length else np.empty((0, 1))
                bu_lam = self.nlp_sol_lam_x[self.index_map['u'][stage]] if bu_length else np.empty((0, 1))
                g_lam = self.nlp_sol_lam_g[self.index_map['lam_g_constraint'][stage]]
            elif stage < dims.N:
                bx_length = self.acados_ocp.constraints.lbx.size
                bu_length = self.acados_ocp.constraints.lbu.size
                bx_lam = self.nlp_sol_lam_x[self.index_map['x'][stage]] if bx_length else np.empty((0, 1))
                bu_lam = self.nlp_sol_lam_x[self.index_map['u'][stage]] if bu_length else np.empty((0, 1))
                g_lam = self.nlp_sol_lam_g[self.index_map['lam_g_constraint'][stage]]
            elif stage == dims.N:
                bx_length = self.acados_ocp.constraints.lbx_e.size
                bu_length = self.acados_ocp.constraints.lbu.size
                bx_lam = self.nlp_sol_lam_x[self.index_map['x'][stage]] if bx_length else np.empty((0, 1))
                bu_lam = self.nlp_sol_lam_x[self.index_map['u'][stage]] if bu_length else np.empty((0, 1))
                g_lam = self.nlp_sol_lam_g[self.index_map['lam_g_constraint'][stage]]

            lbx_lam = np.maximum(0,-bx_lam)
            lbu_lam = np.maximum(0,-bu_lam)
            lg_lam = np.maximum(0,-g_lam)
            ubx_lam = np.maximum(0,bx_lam)
            ubu_lam = np.maximum(0,bu_lam)
            ug_lam = np.maximum(0,g_lam)
            lam = np.concatenate((lbu_lam, lbx_lam, lg_lam, ubu_lam, ubx_lam, ug_lam))
            return lam.flatten()
        else:
            raise NotImplementedError(f"Field '{field}' is not implemented in AcadosCasadiOcpSolver")

    def get_flat(self, field_: str) -> np.ndarray:
        """
        Get concatenation of all stages of last solution of the solver.

        :param field: string in ['x', 'u', 'pi','lam']

        .. note:: The parameter 'p_global' has no stage-wise structure and is processed in a memory saving manner by default. \n
                In order to read the 'p_global' parameter, the option 'save_p_global' must be set to 'True' upon instantiation. \n
        """
        if self.nlp_sol is None:
            raise ValueError('No solution available. Please call solve() first.')
        dims = self.acados_ocp.dims
        result = []

        if field_ == 'x':
            for i in range(dims.N+1):
                result.append(self.get(i, 'x'))
            return np.concatenate(result)
        elif field_ == 'u':
            for i in range(dims.N):
                result.append(self.get(i, 'u'))
            return np.concatenate(result)
        elif field_ == 'lam':
            for i in range(dims.N):
                result.append(self.get(i, 'lam'))
            return np.concatenate(result)
        elif field_ == 'pi':
            for i in range(dims.N):
                result.append(self.get(i, 'pi'))
            return np.concatenate(result)
        elif field_ == 'lam_x':
            return self.nlp_sol_lam_x.flatten()
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
        :value_:
        """
        dims = self.acados_ocp.dims

        if field_ == 'x':
            nx = dims.nx
            for i in range(dims.N+1):
                self.set(i, 'x', value_[i*nx:(i+1)*nx])
        elif field_ == 'u':
            nu = dims.nu
            for i in range(dims.N):
                self.set(i, 'u', value_[i*nu:(i+1)*nu])
        elif field_ == 'pi':
            nx = dims.nx
            for i in range(dims.N):
                self.set(i, 'pi', value_[i*nx:(i+1)*nx])
        elif field_ == 'lam':
            for i in range(dims.N):
                if i == 0:
                    self.set(i, 'lam', value_[i*(dims.nx+dims.nu+dims.nh_0+dims.nphi_0)
                                              :(i+1)*(dims.nx+dims.nu+dims.nh_0+dims.nphi_0)])
                elif i < dims.N:
                    self.set(i, 'lam', value_[i*(dims.nx+dims.nu+dims.nh+dims.nphi)
                                              :(i+1)*(dims.nx+dims.nu+dims.nh+dims.nphi)])
                elif i == dims.N:
                    self.set(i, 'lam', value_[i*(dims.nx+dims.nu+dims.nh_e+dims.nphi_e)
                                              :(i+1)*(dims.nx+dims.nu+dims.nh_e+dims.nphi_e)])
        else:
            raise NotImplementedError(f"Field '{field_}' is not implemented in set_flat().")

    # def load_iterate(self, filename:str, verbose: bool = True):
        raise NotImplementedError()

    def store_iterate_to_obj(self) -> dict:
        """
        Returns the current iterate of the OCP solver as an dictionary.
        """
        self.last_iterate = dict()
        for field in ['x', 'u', 'lam_x', 'lam_g', 'lam_p']:
            traj = []
            for n in range(self.acados_ocp.dims.N+1):
                if n < self.acados_ocp.dims.N or not (field in ["u"]):
                    traj.append(self.get(n, field))

            self.last_iterate[f"{field}_traj"] = traj

        return self.last_iterate

    def load_iterate_from_obj(self, iterate: dict) -> None:
        """
        Loads the provided iterate into the OCP solver.
        Note: The iterate object does not contain the the parameters.
        """
        for key, traj in iterate.items():
            field = key.replace('_traj', '')
            for n, val in enumerate(traj):
                self.set(n, field, val)

    def store_iterate_to_flat_obj(self) -> AcadosOcpFlattenedIterate:
        """
        Returns the current iterate of the OCP solver as an AcadosOcpFlattenedIterate.
        """
        self.acados_ocp.constraints.x0 = np.concatenate(self.get_flat('x'),self.get_flat('u')).flatten()
        self.laxm_x0 = self.get_flat('lam_x')
        self.lam_g0 = self.get_flat('lam_g')
        return self.last_iterate

    def load_iterate_from_flat_obj(self, iterate: dict) -> None:
        """
        Loads the provided iterate into the OCP solver.
        Note: The iterate object does not contain the the parameters.
        """
        self.set_flat('x')
        self.set_flat('u')
        self.set_flat('lam')
        self.set_flat('pi')

    # def get_stats(self, field_: str) -> Union[int, float, np.ndarray]:
        raise NotImplementedError()

    # def get_cost(self) -> float:
        raise NotImplementedError()

    def set(self, stage: int, field: str, value_: np.ndarray):
        """
        Set solver initialization to stages.
        
        :param stage: integer corresponding to shooting node
        :param field: string in ['x', 'u', 'pi', 'lam']
        :value_: 
        """
        dims = self.acados_ocp.dims
        offset_w = stage*(dims.nx+dims.nu)
        offset_lam = stage*dims.nx

        if field == 'x':
            if value_.size == dims.nx:
                self.w0[offset_w:offset_w+dims.nx] = value_.flatten()
            else:
                raise ValueError(f"Expected {dims.nx} values for 'x', got {value_.size}")
        elif field == 'u':
            if value_.size == dims.nu:
                self.w0[offset_w+dims.nx:offset_w+dims.nx+dims.nu] = value_.flatten()
            else:
                raise ValueError(f"Expected {dims.nu} values for 'u', got {value_.size}")
        elif field == 'pi':
            if value_.size == dims.nx:
                self.lam_g0_dynamic[offset_lam:offset_lam+dims.nx] = -value_.flatten()
            else:
                raise ValueError(f"Expected {dims.nx} values for 'pi', got {value_.size}")
        elif field == 'lam':
            if stage == 0:
                offset_l = stage*(dims.nx+dims.nu+dims.nh_0+dims.nphi_0)
                offset_u = (stage+1)*(dims.nx+dims.nu+dims.nh_0+dims.nphi_0)

                bx_length = self.acados_ocp.constraints.lbx_0.size
                bu_length = self.acados_ocp.constraints.lbu.size

                lbu_lam = value_[offset_l:offset_l+bu_length]
                lbx_lam = value_[offset_l+dims.nu:offset_l+dims.nu+bx_length]
                lg_lam = value_[offset_l+dims.nu+dims.nx:offset_l+dims.nx+dims.nu+dims.nh_0+dims.nphi_0]

                ubu_lam = value_[offset_u:offset_u+bu_length]
                ubx_lam = value_[offset_u+dims.nu:offset_u+dims.nu+bx_length]
                ug_lam = value_[offset_u+dims.nu+dims.nx:offset_u+dims.nx+dims.nu+dims.nh_0+dims.nphi_0]

                self.lam_x0[offset_l:offset_l+bu_length+bx_length] = np.concatenate((ubu_lam-lbu_lam, ubx_lam-lbx_lam))
                self.lam_g0_constraint[stage:stage+dims.nh_0+dims.nphi_0] =  ug_lam-lg_lam
            
            elif stage < dims.N:
                offset_l = stage*(dims.nx+dims.nu+dims.nh+dims.nphi)
                offset_u = (stage+1)*(dims.nx+dims.nu+dims.nh+dims.nphi)

                bx_length = self.acados_ocp.constraints.lbx.size
                bu_length = self.acados_ocp.constraints.lbu.size

                lbu_lam = value_[offset_l:offset_l+bu_length]
                lbx_lam = value_[offset_l+dims.nu:offset_l+dims.nu+bx_length]
                lg_lam = value_[offset_l+dims.nu+dims.nx:offset_l+dims.nx+dims.nu+dims.nh+dims.nphi]

                ubu_lam = value_[offset_u:offset_u+bu_length]
                ubx_lam = value_[offset_u+dims.nu:offset_u+dims.nu+bx_length]
                ug_lam = value_[offset_u+dims.nu+dims.nx:offset_u+dims.nx+dims.nu+dims.nh+dims.nphi]

                self.lam_x0[offset_l:offset_l+bu_length+bx_length] = np.concatenate((ubu_lam-lbu_lam, ubx_lam-lbx_lam))
                self.lam_g0_dynamic[offset_l:offset_l+dims.nh+dims.nphi] = ug_lam-lg_lam

            elif stage == dims.N:
                offset_l = stage*(dims.nx+dims.nu+dims.nh+dims.nphi_e)
                offset_u = (stage+1)*(dims.nx+dims.nu+dims.nh+dims.nphi_e)

                bx_length = self.acados_ocp.constraints.lbx_e.size
                bu_length = self.acados_ocp.constraints.lbu.size

                lbu_lam = value_[offset_l:offset_l+bu_length]
                lbx_lam = value_[offset_l+dims.nu:offset_l+dims.nu+bx_length]
                lg_lam = value_[offset_l+dims.nu+dims.nx:offset_l+dims.nx+dims.nu+dims.nh+dims.nphi_e]

                ubu_lam = value_[offset_u:offset_u+bu_length]
                ubx_lam = value_[offset_u+dims.nu:offset_u+dims.nu+bx_length]
                ug_lam = value_[offset_u+dims.nu+dims.nx:offset_u+dims.nx+dims.nu+dims.nh+dims.nphi_e]

                self.lam_x0[offset_l:offset_l+bu_length+bx_length] = np.concatenate((ubu_lam-lbu_lam, ubx_lam-lbx_lam))
                self.lam_g0_dynamic[offset_l:offset_l+dims.nh_e+dims.nphi_e] = ug_lam-lg_lam

    # def cost_get(self, stage_: int, field_: str) -> np.ndarray:
        raise NotImplementedError()

    # def cost_set(self, stage_: int, field_: str, value_):
        raise NotImplementedError()