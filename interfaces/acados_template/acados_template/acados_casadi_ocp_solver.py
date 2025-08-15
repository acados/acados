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

from typing import Union, Optional, List

import numpy as np

from .utils import casadi_length, is_casadi_SX
from .acados_ocp import AcadosOcp
from .acados_ocp_iterate import AcadosOcpIterate, AcadosOcpFlattenedIterate

class AcadosCasadiOcp:

    def __init__(self, ocp: AcadosOcp, with_hessian=False, multiple_shooting=True):
        """
        Creates an equivalent CasADi NLP formulation of the OCP.
        Experimental, not fully implemented yet.

        :return: nlp_dict, bounds_dict, w0 (initial guess)
        """
        ocp.make_consistent()

        # create index map for variables
        self._index_map = {
            # indices of variables within w
            'x_in_w': [],
            'u_in_w': [],
            # indices of slack variables within w
            'sl_in_w': [],
            'su_in_w': [],
            # indices of parameters within p_nlp
            'p_in_p_nlp': [],
            'p_global_in_p_nlp': [],
            # indices of state bounds and control bounds within lam_x(lam_w) in casadi formulation
            'lam_bx_in_lam_w':[],
            'lam_bu_in_lam_w': [],
            # indices of dynamic constraints within g in casadi formulation
            'pi_in_lam_g': [],
            # indicies to [g, h, phi] in acados formulation within lam_g in casadi formulation
            'lam_gnl_in_lam_g': [],
            # indices of slack variables within lam_g in casadi formulation
            'lam_sl_in_lam_g': [],
            'lam_su_in_lam_g': [],
        }
        self.offset_w = 0  # offset for the indices in index_map
        self.offset_gnl = 0
        self.offset_lam = 0

        # unpack
        model = ocp.model
        dims = ocp.dims
        constraints = ocp.constraints
        cost = ocp.cost
        solver_options = ocp.solver_options
        N_horizon = solver_options.N_horizon

        # check what is not supported yet
        if any([dims.nsbx, dims.nsbx_e, dims.nsbu]):
            raise NotImplementedError("AcadosCasadiOcpSolver does not support slack variables (s) for variables (x and u) yet.")
        if any([dims.nsg, dims.nsg_e, dims.nsphi, dims.nsphi_e]):
            raise NotImplementedError("AcadosCasadiOcpSolver does not support slack variables (s)  for general linear and convex-over-nonlinear constraints (g, phi).")
        if dims.nz > 0:
            raise NotImplementedError("AcadosCasadiOcpSolver does not support algebraic variables (z) yet.")
        if ocp.solver_options.integrator_type not in ["DISCRETE", "ERK"]:
            raise NotImplementedError(f"AcadosCasadiOcpSolver does not support integrator type {ocp.solver_options.integrator_type} yet.")

        # create primal variables and slack variables
        ca_symbol = model.get_casadi_symbol()
        xtraj_nodes = []
        utraj_nodes = []
        sl_nodes = []
        su_nodes = []
        if multiple_shooting:
            for i in range(N_horizon+1):
                self._append_node('x', ca_symbol, xtraj_nodes, i, dims)
                self._append_node('u', ca_symbol, utraj_nodes, i, dims)
                self._append_node('sl', ca_symbol, sl_nodes, i, dims)
                self._append_node('su', ca_symbol, su_nodes, i, dims)
        else: # single_shooting
            self._x_traj_fun = []
            for i in range(N_horizon):
                self._append_node('u', ca_symbol, utraj_nodes, i, dims)
                self._append_node('sl', ca_symbol, sl_nodes, i, dims)
                self._append_node('su', ca_symbol, su_nodes, i, dims)

        # parameters
        ptraj_nodes = [ca_symbol(f'p{i}', dims.np, 1) for i in range(N_horizon+1)]

        # setup state and control bounds
        lb_xtraj_nodes = [-np.inf * ca.DM.ones((dims.nx, 1)) for _ in range(N_horizon+1)]
        ub_xtraj_nodes = [np.inf * ca.DM.ones((dims.nx, 1)) for _ in range(N_horizon+1)]
        lb_utraj_nodes = [-np.inf * ca.DM.ones((dims.nu, 1)) for _ in range(N_horizon)]
        ub_utraj_nodes = [np.inf * ca.DM.ones((dims.nu, 1)) for _ in range(N_horizon)]
        # setup slack variables
        # TODO: speicify different bounds for lsbu, lsbx, lsg, lsh ,lsphi
        lb_slack_nodes = ([0 * ca.DM.ones((dims.ns_0, 1))] if dims.ns_0 else []) + \
                        ([0* ca.DM.ones((dims.ns, 1)) for _ in range(N_horizon-1)] if dims.ns else []) + \
                        ([0 * ca.DM.ones((dims.ns_e, 1))] if dims.ns_e else [])
        ub_slack_nodes = ([np.inf * ca.DM.ones((dims.ns, 1))] if dims.ns_0 else []) + \
                        ([np.inf * ca.DM.ones((dims.ns, 1)) for _ in range(N_horizon-1)] if dims.ns else []) + \
                        ([np.inf * ca.DM.ones((dims.ns_e, 1))] if dims.ns_e else [])
        if multiple_shooting:
            for i in range(0, N_horizon+1):
                self._set_bounds_indices('x', i, lb_xtraj_nodes, ub_xtraj_nodes, constraints, dims)
                self._set_bounds_indices('u', i, lb_utraj_nodes, ub_utraj_nodes, constraints, dims)
        else: # single_shooting
            for i in range(0, N_horizon):
                self._set_bounds_indices('u', i, lb_utraj_nodes, ub_utraj_nodes, constraints, dims)

        ### Concatenate primal variables and bounds
        # w = [x0, u0, sl0, su0, x1, u1, ...]
        w_sym_list = []
        lbw_list = []
        ubw_list = []
        w0_list = []
        p_list = []
        offset_p = 0
        x_guess = ocp.constraints.x0 if ocp.constraints.has_x0 else np.zeros((dims.nx,))
        if multiple_shooting:
            for i in range(N_horizon+1):
                self._append_variables_and_bounds('x', w_sym_list, lbw_list, ubw_list, w0_list, xtraj_nodes, lb_xtraj_nodes, ub_xtraj_nodes, i, dims, x_guess)
                if i < N_horizon:
                    self._append_variables_and_bounds('u',w_sym_list, lbw_list, ubw_list, w0_list, utraj_nodes, lb_utraj_nodes, ub_utraj_nodes, i, dims, x_guess)
                self._append_variables_and_bounds('slack', w_sym_list, lbw_list, ubw_list, w0_list, [sl_nodes, su_nodes], lb_slack_nodes, ub_slack_nodes, i, dims, x_guess)
                p_list.append(ocp.parameter_values)
                self._index_map['p_in_p_nlp'].append(list(range(offset_p, offset_p+dims.np)))
                offset_p += dims.np
        else: # single_shooting
            xtraj_nodes.append(x_guess)
            self._x_traj_fun.append(x_guess)
            for i in range(N_horizon+1):
                if i < N_horizon:
                    self._append_variables_and_bounds('u', w_sym_list, lbw_list, ubw_list, w0_list, utraj_nodes, lb_utraj_nodes, ub_utraj_nodes, i, dims, x_guess)
                self._append_variables_and_bounds('slack', w_sym_list, lbw_list, ubw_list, w0_list, [sl_nodes, su_nodes], lb_slack_nodes, ub_slack_nodes, i, dims, x_guess)
                p_list.append(ocp.parameter_values)
                self._index_map['p_in_p_nlp'].append(list(range(offset_p, offset_p+dims.np)))
                offset_p += dims.np
        p_list.append(ocp.p_global_values)
        self._index_map['p_global_in_p_nlp'].append(list(range(offset_p, offset_p+dims.np_global)))
        offset_p += dims.np_global

        nw = self.offset_w  # number of primal variables

        # vectorize
        w = ca.vertcat(*w_sym_list)
        lbw = ca.vertcat(*lbw_list)
        ubw = ca.vertcat(*ubw_list)
        p_nlp = ca.vertcat(*ptraj_nodes, model.p_global)

        ### Create Constraints
        g = []
        lbg = []
        ubg = []
        if with_hessian:
            lam_g = []
            hess_l = ca.DM.zeros((nw, nw))
        # dynamics constraints
        if solver_options.integrator_type == "DISCRETE":
            f_discr_fun = ca.Function('f_discr_fun', [model.x, model.u, model.p, model.p_global], [model.disc_dyn_expr])
        elif solver_options.integrator_type == "ERK":
            param = ca.vertcat(model.u, model.p, model.p_global)
            ca_expl_ode = ca.Function('ca_expl_ode', [model.x, param], [model.f_expl_expr])
            f_discr_fun = ca.simpleRK(ca_expl_ode, solver_options.sim_method_num_steps[0], solver_options.sim_method_num_stages[0])
        else:
            raise NotImplementedError(f"Integrator type {solver_options.integrator_type} not supported.")

        for i in range(N_horizon+1):
            # add dynamics constraints
            if multiple_shooting:
                if i < N_horizon:
                    if solver_options.integrator_type == "DISCRETE":
                        dyn_equality = xtraj_nodes[i+1] - f_discr_fun(xtraj_nodes[i], utraj_nodes[i], ptraj_nodes[i], model.p_global)
                    elif solver_options.integrator_type == "ERK":
                        param = ca.vertcat(utraj_nodes[i], ptraj_nodes[i], model.p_global)
                        dyn_equality = xtraj_nodes[i+1] - f_discr_fun(xtraj_nodes[i], param, solver_options.time_steps[i])
                    self._append_constraints(i, 'dyn', g, lbg, ubg,
                                            g_expr = dyn_equality,
                                            lbg_expr = np.zeros((dims.nx, 1)),
                                            ubg_expr = np.zeros((dims.nx, 1)),
                                            cons_dim=dims.nx)
                    if with_hessian:
                        # add hessian of dynamics constraints
                        lam_g_dyn = ca_symbol(f'lam_g_dyn{i}', dims.nx, 1)
                        lam_g.append(lam_g_dyn)
                        if ocp.solver_options.hessian_approx == 'EXACT' and ocp.solver_options.exact_hess_dyn:
                            adj = ca.jtimes(dyn_equality, w, lam_g_dyn, True)
                            hess_l += ca.jacobian(adj, w, {"symmetric": is_casadi_SX(model.x)})
            else: # single_shooting
                if i < N_horizon:
                    x_current = xtraj_nodes[i]
                    if solver_options.integrator_type == "DISCRETE":
                        x_next = f_discr_fun(x_current, utraj_nodes[i], ptraj_nodes[i], model.p_global)
                    elif solver_options.integrator_type == "ERK":
                        param = ca.vertcat(utraj_nodes[i], ptraj_nodes[i], model.p_global)
                        x_next = f_discr_fun(x_current, param, solver_options.time_steps[i])
                    xtraj_nodes.append(x_next)
                    self._x_traj_fun.append(f_discr_fun)
            # Nonlinear Constraints
            # initial stage
            lg, ug, lh, uh, lphi, uphi, ng, nh, nphi, nsg, nsh, nsphi, idxsh, linear_constr_expr, h_i_nlp_expr, conl_constr_fun =\
            self._get_constraint_node(i, N_horizon, xtraj_nodes, utraj_nodes, ptraj_nodes, model, constraints, dims)

            # add linear constraints
            if ng > 0:
                self._append_constraints(i, 'gnl', g, lbg, ubg,
                                         g_expr = linear_constr_expr,
                                         lbg_expr = lg,
                                         ubg_expr = ug,
                                         cons_dim=ng)

            # add nonlinear constraints
            if nh > 0:
                if nsh > 0:
                    # h_fun with slack variables
                    soft_h_indices = idxsh
                    hard_h_indices = np.array([h for h in range(len(lh)) if h not in idxsh])
                    for index_in_nh in range(nh):
                        if index_in_nh in soft_h_indices:
                            index_in_soft = soft_h_indices.tolist().index(index_in_nh)
                            self._append_constraints(i, 'gnl', g, lbg, ubg,
                                                     g_expr = h_i_nlp_expr[index_in_nh] + sl_nodes[i][index_in_soft],
                                                     lbg_expr = lh[index_in_nh],
                                                     ubg_expr = np.inf * ca.DM.ones((1, 1)),
                                                     cons_dim=1,
                                                     sl=True)
                            self._append_constraints(i, 'gnl', g, lbg, ubg,
                                                     g_expr = h_i_nlp_expr[index_in_nh] - su_nodes[i][index_in_soft],
                                                     lbg_expr = -np.inf * ca.DM.ones((1, 1)),
                                                     ubg_expr = uh[index_in_nh],
                                                     cons_dim=1,
                                                     su=True)
                        elif index_in_nh in hard_h_indices:
                            self._append_constraints(i, 'gnl', g, lbg, ubg,
                                                     g_expr = h_i_nlp_expr[index_in_nh],
                                                     lbg_expr = lh[index_in_nh],
                                                     ubg_expr = uh[index_in_nh],
                                                     cons_dim=1)
                else:
                    self._append_constraints(i, 'gnl', g, lbg, ubg,
                                             g_expr = h_i_nlp_expr,
                                             lbg_expr = lh,
                                             ubg_expr = uh,
                                             cons_dim=nh)
                if with_hessian:
                    # add hessian contribution
                    lam_h = ca_symbol(f'lam_h_{i}', dims.nh, 1)
                    lam_g.append(lam_h)
                    if ocp.solver_options.hessian_approx == 'EXACT' and ocp.solver_options.exact_hess_constr:
                        adj = ca.jtimes(h_i_nlp_expr, w, lam_h, True)
                        hess_l += ca.jacobian(adj, w, {"symmetric": is_casadi_SX(model.x)})

            # add convex-over-nonlinear constraints
            if nphi > 0:
                self._append_constraints(i, 'gnl', g, lbg, ubg,
                                         g_expr = conl_constr_fun(xtraj_nodes[i], utraj_nodes[i], ptraj_nodes[i], model.p_global),
                                         lbg_expr = lphi,
                                         ubg_expr = uphi,
                                         cons_dim=nphi)
                if with_hessian:
                    lam_phi = ca_symbol(f'lam_phi', nphi, 1)
                    lam_g.append(lam_phi)
                    # always use CONL Hessian approximation here, disregarding inner second derivative
                    outer_hess_r = ca.vertcat(*[ca.hessian(model.con_phi_expr[i], model.con_r_in_phi)[0] for i in range(dims.nphi)])
                    outer_hess_r = ca.substitute(outer_hess_r, model.con_r_in_phi, model.con_r_expr)
                    r_in_nlp = ca.substitute(model.con_r_expr, model.x, xtraj_nodes[-1])
                    dr_dw = ca.jacobian(r_in_nlp, w)
                    hess_l += dr_dw.T @ outer_hess_r @ dr_dw

        ### Cost
        nlp_cost = 0
        for i in range(N_horizon+1):
            xtraj_node_i, utraj_node_i, ptraj_node_i, sl_node_i, su_node_i, cost_expr_i, ns, zl, Zl, zu, Zu = \
            self._get_cost_node(i, N_horizon, xtraj_nodes, utraj_nodes, ptraj_nodes, sl_nodes, su_nodes, ocp, dims, cost)

            cost_fun_i = ca.Function(f'cost_fun_{i}', [model.x, model.u, model.p, model.p_global], [cost_expr_i])
            nlp_cost += solver_options.cost_scaling[i] * cost_fun_i(xtraj_node_i, utraj_node_i, ptraj_node_i, model.p_global)
            if ns:
                penalty_expr_i = 0.5 * ca.mtimes(sl_node_i.T, ca.mtimes(np.diag(Zl), sl_node_i)) + \
                    ca.mtimes(zl.reshape(-1, 1).T, sl_node_i) + \
                    0.5 * ca.mtimes(su_node_i.T, ca.mtimes(np.diag(Zu), su_node_i)) + \
                    ca.mtimes(zu.reshape(-1, 1).T, su_node_i)
                nlp_cost += solver_options.cost_scaling[i] * penalty_expr_i

        if with_hessian:
            lam_f = ca_symbol('lam_f', 1, 1)
            if ocp.solver_options.hessian_approx == 'EXACT' or \
                (ocp.cost.cost_type == "LINEAR_LS" and ocp.cost.cost_type_0 == "LINEAR_LS" and ocp.cost.cost_type_e == "LINEAR_LS"):
                hess_l += lam_f * ca.hessian(nlp_cost, w)[0]
            else:
                raise NotImplementedError("Hessian approximation not implemented for this cost type.")
            lam_g_vec = ca.vertcat(*lam_g)
            nlp_hess_l_custom = ca.Function('nlp_hess_l', [w, p_nlp, lam_f, lam_g_vec], [ca.triu(hess_l)])
            assert casadi_length(lam_g_vec) == casadi_length(ca.vertcat(*g)), f"Number of nonlinear constraints does not match the expected number, got {casadi_length(lam_g_vec)} != {casadi_length(ca.vertcat(*g))}."
        else:
            nlp_hess_l_custom = None
            hess_l = None

        # create NLP
        nlp = {"x": w, "p": p_nlp, "g": ca.vertcat(*g), "f": nlp_cost}
        bounds = {"lbx": lbw, "ubx": ubw, "lbg": ca.vertcat(*lbg), "ubg": ca.vertcat(*ubg)}
        w0 = np.concatenate(w0_list)
        p = np.concatenate(p_list)

        self.__nlp = nlp
        self.__bounds = bounds
        self.__w0 = w0
        self.__p = p
        self.__index_map = self._index_map
        self.__nlp_hess_l_custom = nlp_hess_l_custom
        self.__hess_approx_expr = hess_l

    def _append_node(self, _field, ca_symbol, node:list, i, dims):
        """
        Helper function to append a node to the NLP formulation.
        """
        if i == 0:
            ns = dims.ns_0
        elif i < dims.N:
            ns = dims.ns
        else:
            ns = dims.ns_e
        if _field == 'x':
            node.append(ca_symbol(f'x{i}', dims.nx, 1))
        elif _field == 'u':
            node.append(ca_symbol(f'u{i}', dims.nu, 1))
        elif _field == 'sl':
            if ns > 0:
                node.append(ca_symbol(f'sl{i}', ns, 1))
            else:
                node.append([])
        elif _field == 'su':
            if ns > 0:
                node.append(ca_symbol(f'su{i}', ns, 1))
            else:
                node.append([])

    def _set_bounds_indices(self, _field, i, lb_node, ub_node, constraints, dims):
        """
        Helper function to set bounds and indices for the primal variables.
        """
        if i == 0:
            idxbx = constraints.idxbx_0
            idxbu = constraints.idxbu
            lbx = constraints.lbx_0
            ubx = constraints.ubx_0
            lbu = constraints.lbu
            ubu = constraints.ubu
        elif i < dims.N:
            idxbx = constraints.idxbx
            idxbu = constraints.idxbu
            lbx = constraints.lbx
            ubx = constraints.ubx
            lbu = constraints.lbu
            ubu = constraints.ubu
        elif i == dims.N:
            idxbx = constraints.idxbx_e
            lbx = constraints.lbx_e
            ubx = constraints.ubx_e
        if i < dims.N:
            if _field == 'x':
                lb_node[i][idxbx] = lbx
                ub_node[i][idxbx] = ubx
                self._index_map['lam_bx_in_lam_w'].append(list(self.offset_lam + idxbx))
                self.offset_lam += dims.nx
            elif _field == 'u':
                lb_node[i][idxbu] = lbu
                ub_node[i][idxbu] = ubu
                self._index_map['lam_bu_in_lam_w'].append(list(self.offset_lam + idxbu))
                self.offset_lam += dims.nu
        elif i == dims.N:
            if _field == 'x':
                lb_node[i][idxbx] = lbx
                ub_node[i][idxbx] = ubx
                self._index_map['lam_bx_in_lam_w'].append(list(self.offset_lam + idxbx))
                self.offset_lam += dims.nx

    def _append_variables_and_bounds(self, _field, w_sym_list, lbw_list, ubw_list, w0_list,
                                     node_list, lb_node_list, ub_node_list, i, dims, x_guess):
        """
        Unified helper function to add a primal or slack variable to the NLP formulation.
        """
        if _field == "x":
            # Add state variable
            w_sym_list.append(node_list[i])
            lbw_list.append(lb_node_list[i])
            ubw_list.append(ub_node_list[i])
            w0_list.append(x_guess)
            self._index_map['x_in_w'].append(list(range(self.offset_w, self.offset_w + dims.nx)))
            self.offset_w += dims.nx

        elif _field == "u":
            # Add control variable
            w_sym_list.append(node_list[i])
            lbw_list.append(lb_node_list[i])
            ubw_list.append(ub_node_list[i])
            w0_list.append(np.zeros((dims.nu,)))
            self._index_map['u_in_w'].append(list(range(self.offset_w, self.offset_w + dims.nu)))
            self.offset_w += dims.nu

        elif _field == "slack":
            # Add slack variables (sl and su)
            if i == 0 and dims.ns_0:
                ns = dims.ns_0
            elif i < dims.N and dims.ns:
                ns = dims.ns
            elif i == dims.N and dims.ns_e:
                ns = dims.ns_e
            else:
                self._index_map['sl_in_w'].append([])
                self._index_map['su_in_w'].append([])
                return

            # Add sl
            w_sym_list.append(node_list[0][i])
            lbw_list.append(lb_node_list[i])
            ubw_list.append(ub_node_list[i])
            w0_list.append(np.zeros((ns,)))
            # Add su
            w_sym_list.append(node_list[1][i])
            lbw_list.append(lb_node_list[i])
            ubw_list.append(ub_node_list[i])
            w0_list.append(np.zeros((ns,)))

            self._index_map['sl_in_w'].append(list(range(self.offset_w, self.offset_w + ns)))
            self._index_map['su_in_w'].append(list(range(self.offset_w + ns, self.offset_w + 2 * ns)))
            self.offset_w += 2 * ns

        else:
            raise ValueError(f"Unsupported for: {_field}")

    def _append_constraints(self, i, _field, g, lbg, ubg, g_expr, lbg_expr, ubg_expr, cons_dim, sl=False, su=False):
        """
        Helper function to append constraints to the NLP formulation.
        """
        g.append(g_expr)
        lbg.append(lbg_expr)
        ubg.append(ubg_expr)
        if _field == 'dyn':
            self._index_map['pi_in_lam_g'].append(list(range(self.offset_gnl, self.offset_gnl + cons_dim)))
            self.offset_gnl += cons_dim
        elif _field == 'gnl':
            if not sl and not su:
                self._index_map['lam_gnl_in_lam_g'][i].extend(list(range(self.offset_gnl, self.offset_gnl + cons_dim)))
                self.offset_gnl += cons_dim
            elif sl:
                self._index_map['lam_sl_in_lam_g'][i].append(self.offset_gnl)
                self.offset_gnl += 1
            elif su:
                self._index_map['lam_su_in_lam_g'][i].append(self.offset_gnl)
                self.offset_gnl += 1

    def _get_cost_node(self, i, N_horizon, xtraj_node, utraj_node, ptraj_node, sl_node, su_node, ocp, dims, cost):
        """
        Helper function to get the cost node for a given stage.
        """
        if i == 0:
            return (xtraj_node[0],
                    utraj_node[0],
                    ptraj_node[0],
                    sl_node[0],
                    su_node[0],
                    ocp.get_initial_cost_expression(),
                    dims.ns_0, cost.zl_0, cost.Zl_0, cost.zu_0, cost.Zu_0)
        elif i < N_horizon:
            return (xtraj_node[i],
                    utraj_node[i],
                    ptraj_node[i],
                    sl_node[i],
                    su_node[i],
                    ocp.get_path_cost_expression(),
                    dims.ns, cost.zl, cost.Zl, cost.zu, cost.Zu)
        else:
            return (xtraj_node[-1],
                    [],
                    ptraj_node[-1],
                    sl_node[-1],
                    su_node[-1],
                    ocp.get_terminal_cost_expression(),
                    dims.ns_e, cost.zl_e, cost.Zl_e, cost.zu_e, cost.Zu_e)

    def _get_constraint_node(self, i, N_horizon, xtraj_node, utraj_node, ptraj_node, model, constraints, dims):
        """
        Helper function to get the constraint node for a given stage.
        """
        if i == 0 and N_horizon > 0:
            lg, ug = constraints.lg, constraints.ug
            lh, uh = constraints.lh_0, constraints.uh_0
            lphi, uphi = constraints.lphi_0, constraints.uphi_0
            ng, nh, nphi = dims.ng, dims.nh_0, dims.nphi_0
            nsg, nsh, nsphi, idxsh = dims.nsg, dims.nsh_0, dims.nsphi_0, constraints.idxsh_0

            # linear function
            linear_constr_expr = None
            if dims.ng > 0:
                C = constraints.C
                D = constraints.D
                linear_constr_expr = ca.mtimes(C, xtraj_node[i]) + ca.mtimes(D, utraj_node[i])
            # nonlinear function
            h_fun = ca.Function('h_0_fun', [model.x, model.u, model.p, model.p_global], [model.con_h_expr_0])
            h_i_nlp_expr = h_fun(xtraj_node[i], utraj_node[i], ptraj_node[i], model.p_global)
            # compound nonlinear constraint
            conl_constr_fun = None
            if dims.nphi_0 > 0:
                conl_expr = ca.substitute(model.con_phi_expr_0, model.con_r_in_phi_0, model.con_r_expr_0)
                conl_constr_fun = ca.Function('conl_constr_0_fun', [model.x, model.u, model.p, model.p_global], [conl_expr])

        elif i < N_horizon:
            lg, ug = constraints.lg, constraints.ug
            lh, uh = constraints.lh, constraints.uh
            lphi, uphi = constraints.lphi, constraints.uphi
            ng, nh, nphi = dims.ng, dims.nh, dims.nphi
            nsg, nsh, nsphi, idxsh = dims.nsg, dims.nsh, dims.nsphi, constraints.idxsh

            linear_constr_expr = None
            if dims.ng > 0:
                C = constraints.C
                D = constraints.D
                linear_constr_expr = ca.mtimes(C, xtraj_node[i]) + ca.mtimes(D, utraj_node[i])
            h_fun = ca.Function('h_fun', [model.x, model.u, model.p, model.p_global], [model.con_h_expr])
            h_i_nlp_expr = h_fun(xtraj_node[i], utraj_node[i], ptraj_node[i], model.p_global)
            conl_constr_fun = None
            if dims.nphi > 0:
                conl_expr = ca.substitute(model.con_phi_expr, model.con_r_in_phi, model.con_r_expr)
                conl_constr_fun = ca.Function('conl_constr_fun', [model.x, model.u, model.p, model.p_global], [conl_expr])

        else:
            lg, ug = constraints.lg_e, constraints.ug_e
            lh, uh = constraints.lh_e, constraints.uh_e
            lphi, uphi = constraints.lphi_e, constraints.uphi_e
            ng, nh, nphi = dims.ng_e, dims.nh_e, dims.nphi_e
            nsg, nsh, nsphi, idxsh = dims.nsg_e, dims.nsh_e, dims.nsphi_e, constraints.idxsh_e

            linear_constr_expr = None
            if dims.ng_e > 0:
                C = constraints.C_e
                linear_constr_expr = ca.mtimes(C, xtraj_node[i])
            h_fun = ca.Function('h_e_fun', [model.x, model.p, model.p_global], [model.con_h_expr_e])
            h_i_nlp_expr = h_fun(xtraj_node[i], ptraj_node[i], model.p_global)
            conl_constr_fun = None
            if dims.nphi_e > 0:
                conl_expr = ca.substitute(model.con_phi_expr_e, model.con_r_in_phi_e, model.con_r_expr_e)
                conl_constr_fun = ca.Function('conl_constr_e_fun', [model.x, model.p, model.p_global], [conl_expr])

        self._index_map['lam_gnl_in_lam_g'].append([])
        self._index_map['lam_sl_in_lam_g'].append([])
        self._index_map['lam_su_in_lam_g'].append([])

        return (
            lg, ug,
            lh, uh,
            lphi, uphi,
            ng, nh, nphi,
            nsg, nsh, nsphi,
            idxsh,
            linear_constr_expr,
            h_i_nlp_expr,
            conl_constr_fun
        )

    @property
    def nlp(self):
        """
        Dict containing all symbolics needed to create a `casadi.nlpsol` solver, namely entries 'x', 'p', 'g', 'f'.
        """
        return self.__nlp

    @property
    def w0(self):
        """
        Default initial guess for primal variable vector w for given NLP.
        """
        return self.__w0

    @property
    def p_nlp_values(self):
        """
        Default parameter vector p_nlp in the form of [p_0,..., p_N, p_global] for given NLP.
        """
        return self.__p

    @property
    def bounds(self):
        """
        Dict containing all bounds needed to call a `casadi.nlpsol` solver.
        """
        return self.__bounds

    @property
    def index_map(self):
        """
        Dict containing indices corresponding to stage-wise values of the original OCP, specifically:
        - 'x_in_w': indices of x variables within primal variable vector w
        - 'u_in_w': indices of u variables within primal variable vector w
        - 'pi_in_lam_g': indices of dynamic constraints within g in casadi formulation
        - 'lam_gnl_in_lam_g' indicies to [g, h, phi] in acados formulation within lam_g in casadi formulation
        """
        return self.__index_map

    @property
    def nlp_hess_l_custom(self):
        """
        CasADi Function that computes the Hessian approximation of the Lagrangian in the format required by `casadi.nlpsol`, i.e. as upper triangular matrix.
        The Hessian is set up to match the Hessian that would be used in acados and depends on the solver options.
        """
        return self.__nlp_hess_l_custom

    @property
    def hess_approx_expr(self):
        """
        CasADi expression corresponding to the Hessian approximation of the Lagrangian.
        Expression corresponding to what is output by the `nlp_hess_l_custom` function.
        """
        return self.__hess_approx_expr

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


    def solve_for_x0(self, x0_bar, fail_on_nonzero_status=True, print_stats_on_failure=True):
        raise NotImplementedError()


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
            return np.empty((0, 1))  # Only empty is supported for now. TODO: extend.
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
        else:
            raise NotImplementedError(f"Field '{field}' is not yet implemented in set().")

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