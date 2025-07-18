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

from typing import Union, Optional

import numpy as np

from .utils import casadi_length, is_casadi_SX
from .acados_ocp import AcadosOcp
from .acados_ocp_iterate import AcadosOcpIterate, AcadosOcpFlattenedIterate

class AcadosCasadiOcp:

    def __init__(self, ocp: AcadosOcp, with_hessian=False):
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
        xtraj_node = []
        utraj_node = []
        sl_node = []
        su_node = []
        for i in range(N_horizon+1):
            self._append_node(ca_symbol, xtraj_node, utraj_node, sl_node, su_node, i, dims)

        # parameters
        ptraj_node = [ca_symbol(f'p{i}', dims.np, 1) for i in range(N_horizon+1)]

        # setup state and control bounds
        lb_xtraj_node = [-np.inf * ca.DM.ones((dims.nx, 1)) for _ in range(N_horizon+1)]
        ub_xtraj_node = [np.inf * ca.DM.ones((dims.nx, 1)) for _ in range(N_horizon+1)]
        lb_utraj_node = [-np.inf * ca.DM.ones((dims.nu, 1)) for _ in range(N_horizon)]
        ub_utraj_node = [np.inf * ca.DM.ones((dims.nu, 1)) for _ in range(N_horizon)]
        # setup slack variables
        # TODO: speicify different bounds for lsbu, lsbx, lsg, lsh ,lsphi
        lb_slack_node = ([0 * ca.DM.ones((dims.ns_0, 1))] if dims.ns_0 else []) + \
                        ([0* ca.DM.ones((dims.ns, 1)) for _ in range(N_horizon-1)] if dims.ns else []) + \
                        ([0 * ca.DM.ones((dims.ns_e, 1))] if dims.ns_e else [])
        ub_slack_node = ([np.inf * ca.DM.ones((dims.ns, 1))] if dims.ns_0 else []) + \
                        ([np.inf * ca.DM.ones((dims.ns, 1)) for _ in range(N_horizon-1)] if dims.ns else []) + \
                        ([np.inf * ca.DM.ones((dims.ns_e, 1))] if dims.ns_e else [])
        for i in range(0, N_horizon+1):
            self._set_bounds_indices(i, lb_xtraj_node, ub_xtraj_node, lb_utraj_node, ub_utraj_node, constraints, dims)

        ### Concatenate primal variables and bounds
        # w = [x0, u0, sl0, su0, x1, u1, ...]
        w_sym_list = []
        lbw_list = []
        ubw_list = []
        w0_list = []
        p_list = []
        offset_p = 0
        x_guess = ocp.constraints.x0 if ocp.constraints.has_x0 else np.zeros((dims.nx,))
        for i in range(N_horizon+1):
            if i < N_horizon:
                # add x
                self._append_variables_and_bounds('x', w_sym_list, lbw_list, ubw_list, w0_list, xtraj_node, lb_xtraj_node, ub_xtraj_node, i, dims, x_guess)
                # add u
                self._append_variables_and_bounds('u',w_sym_list, lbw_list, ubw_list, w0_list, utraj_node, lb_utraj_node, ub_utraj_node, i, dims, x_guess)
                # add slack variables
                self._append_variables_and_bounds('slack', w_sym_list, lbw_list, ubw_list, w0_list, [sl_node, su_node], lb_slack_node, ub_slack_node, i, dims, x_guess)
                # add parameters
                p_list.append(ocp.parameter_values)
                self._index_map['p_in_p_nlp'].append(list(range(offset_p, offset_p+dims.np)))
                offset_p += dims.np
            else:
                ## terminal stage
                # add x
                self._append_variables_and_bounds('x', w_sym_list, lbw_list, ubw_list, w0_list, xtraj_node, lb_xtraj_node, ub_xtraj_node, i, dims, x_guess)
                # add slack variables
                self._append_variables_and_bounds('slack', w_sym_list, lbw_list, ubw_list, w0_list, [sl_node, su_node], lb_slack_node, ub_slack_node, i, dims, x_guess)
                # add parameters
                p_list.append(ocp.parameter_values)
                self._index_map['p_in_p_nlp'].append(list(range(offset_p, offset_p+dims.np)))
                offset_p += dims.np
                # add global parameters
                p_list.append(ocp.p_global_values)
                self._index_map['p_global_in_p_nlp'].append(list(range(offset_p, offset_p+dims.np_global)))
                offset_p += dims.np_global

        nw = self.offset_w  # number of primal variables

        # vectorize
        w = ca.vertcat(*w_sym_list)
        lbw = ca.vertcat(*lbw_list)
        ubw = ca.vertcat(*ubw_list)
        p_nlp = ca.vertcat(*ptraj_node, model.p_global)

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
            para = ca.vertcat(model.u, model.p, model.p_global)
            ca_expl_ode = ca.Function('ca_expl_ode', [model.x, para], [model.f_expl_expr])
            f_discr_fun = ca.simpleRK(ca_expl_ode, solver_options.sim_method_num_steps[0], solver_options.sim_method_num_stages[0])
        else:
            raise NotImplementedError(f"Integrator type {solver_options.integrator_type} not supported.")

        for i in range(N_horizon+1):
            # add dynamics constraints
            if i < N_horizon:
                if solver_options.integrator_type == "DISCRETE":
                    dyn_equality = xtraj_node[i+1] - f_discr_fun(xtraj_node[i], utraj_node[i], ptraj_node[i], model.p_global)
                elif solver_options.integrator_type == "ERK":
                    para = ca.vertcat(utraj_node[i], ptraj_node[i], model.p_global)
                    dyn_equality = xtraj_node[i+1] - f_discr_fun(xtraj_node[i], para, solver_options.time_steps[i])
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

            # Nonlinear Constraints
            # initial stage
            lg, ug, lh, uh, lphi, uphi, ng, nh, nphi, nsg, nsh, nsphi, idxsh, linear_constr_expr, h_i_nlp_expr, conl_constr_fun =\
            self._get_constraint_node(i, N_horizon, xtraj_node, utraj_node, ptraj_node, model, constraints, dims)

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
                                                     g_expr = h_i_nlp_expr[index_in_nh] + sl_node[i][index_in_soft],
                                                     lbg_expr = lh[index_in_nh],
                                                     ubg_expr = np.inf * ca.DM.ones((1, 1)),
                                                     cons_dim=1,
                                                     sl=True)
                            self._append_constraints(i, 'gnl', g, lbg, ubg,
                                                     g_expr = h_i_nlp_expr[index_in_nh] - su_node[i][index_in_soft],
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

            # add compound nonlinear constraints
            if nphi > 0:
                self._append_constraints(i, 'gnl', g, lbg, ubg,
                                         g_expr = conl_constr_fun(xtraj_node[i], utraj_node[i], ptraj_node[i], model.p_global),
                                         lbg_expr = lphi,
                                         ubg_expr = uphi,
                                         cons_dim=nphi)
                if with_hessian:
                    lam_phi = ca_symbol(f'lam_phi', nphi, 1)
                    lam_g.append(lam_phi)
                    # always use CONL Hessian approximation here, disregarding inner second derivative
                    outer_hess_r = ca.vertcat(*[ca.hessian(model.con_phi_expr[i], model.con_r_in_phi)[0] for i in range(dims.nphi)])
                    outer_hess_r = ca.substitute(outer_hess_r, model.con_r_in_phi, model.con_r_expr)
                    r_in_nlp = ca.substitute(model.con_r_expr, model.x, xtraj_node[-1])
                    dr_dw = ca.jacobian(r_in_nlp, w)
                    hess_l += dr_dw.T @ outer_hess_r @ dr_dw

        ### Cost
        nlp_cost = 0
        for i in range(N_horizon+1):
            xtraj_node_i, utraj_node_i, ptraj_node_i, sl_node_i, su_node_i, cost_expr_i, ns, zl, Zl, zu, Zu = \
            self._get_cost_node(i, N_horizon, xtraj_node, utraj_node, ptraj_node, sl_node, su_node, ocp, dims, cost)

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

    def _append_node(self, ca_symbol, xtraj_node:list, utraj_node:list, sl_node:list, su_node:list, i, dims):
        """
        Helper function to append a node to the NLP formulation.
        """
        if i == 0:
            ns = dims.ns_0
        elif i < dims.N:
            ns = dims.ns
        else:
            ns = dims.ns_e
        xtraj_node.append(ca_symbol(f'x{i}', dims.nx, 1))
        utraj_node.append(ca_symbol(f'u{i}', dims.nu, 1))
        if ns > 0:
            sl_node.append(ca_symbol(f'sl_0', ns, 1))
            su_node.append(ca_symbol(f'su_0', ns, 1))
        else:
            sl_node.append([])
            su_node.append([])

    def _set_bounds_indices(self, i, lb_xtraj_node, ub_xtraj_node, lb_utraj_node, ub_utraj_node, constraints, dims):
        """
        Helper function to set bounds and indices for the primal variables.
        """
        if i == 0:
            lb_xtraj_node[i][constraints.idxbx_0] = constraints.lbx_0
            ub_xtraj_node[i][constraints.idxbx_0] = constraints.ubx_0
            self._index_map['lam_bx_in_lam_w'].append(list(self.offset_lam + constraints.idxbx_0))
            self.offset_lam += dims.nx
        elif i < dims.N:
            lb_xtraj_node[i][constraints.idxbx] = constraints.lbx
            ub_xtraj_node[i][constraints.idxbx] = constraints.ubx
            self._index_map['lam_bx_in_lam_w'].append(list(self.offset_lam + constraints.idxbx))
            self.offset_lam += dims.nx
        elif i == dims.N:
            lb_xtraj_node[-1][constraints.idxbx_e] = constraints.lbx_e
            ub_xtraj_node[-1][constraints.idxbx_e] = constraints.ubx_e
            self._index_map['lam_bx_in_lam_w'].append(list(self.offset_lam + constraints.idxbx_e))
            self.offset_lam += dims.nx
        if i < dims.N:
            lb_utraj_node[i][constraints.idxbu] = constraints.lbu
            ub_utraj_node[i][constraints.idxbu] = constraints.ubu
            self._index_map['lam_bu_in_lam_w'].append(list(self.offset_lam + constraints.idxbu))
            self.offset_lam += dims.nu
            self.offset_lam += 2*dims.ns_0 if i == 0 else 2*dims.ns

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
                 use_acados_hessian: bool = False):

        if not isinstance(ocp, AcadosOcp):
            raise TypeError('ocp should be of type AcadosOcp.')

        self.ocp = ocp

        # create casadi NLP formulation
        casadi_nlp_obj = AcadosCasadiOcp(ocp, with_hessian=use_acados_hessian)

        self.acados_casadi_ocp = casadi_nlp_obj

        self.casadi_nlp = casadi_nlp_obj.nlp
        self.bounds = casadi_nlp_obj.bounds
        self.w0 = casadi_nlp_obj.w0
        self.p = casadi_nlp_obj.p_nlp_values
        self.index_map = casadi_nlp_obj.index_map
        self.nlp_hess_l_custom = casadi_nlp_obj.nlp_hess_l_custom

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
        self.lam_x0 = np.empty(self.casadi_nlp['x'].shape).flatten()
        self.lam_g0 = np.empty(self.casadi_nlp['g'].shape).flatten()
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
        self.nlp_sol_lam_g = self.nlp_sol['lam_g'].full()
        self.nlp_sol_lam_x = self.nlp_sol['lam_x'].full()

        # statistics
        solver_stats = self.casadi_solver.stats()
        # timing = solver_stats['t_proc_total']
        self.status = solver_stats['return_status']
        self.nlp_iter = solver_stats['iter_count']
        self.time_total = solver_stats['t_wall_total']
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
        if field == 'x':
            return self.nlp_sol_w[self.index_map['x_in_w'][stage]].flatten()
        elif field == 'u':
            return self.nlp_sol_w[self.index_map['u_in_w'][stage]].flatten()
        elif field == 'pi':
            return -self.nlp_sol_lam_g[self.index_map['pi_in_lam_g'][stage]].flatten()
        elif field == 'p':
            return self.p[self.index_map['p_in_p_nlp'][stage]].flatten()
        elif field == 'sl':
            return self.nlp_sol_w[self.index_map['sl_in_w'][stage]].flatten()
        elif field == 'su':
            return self.nlp_sol_w[self.index_map['su_in_w'][stage]].flatten()
        elif field == 'lam':
            if stage == 0:
                bx_lam = self.nlp_sol_lam_x[self.index_map['lam_bx_in_lam_w'][stage]]
                bu_lam = self.nlp_sol_lam_x[self.index_map['lam_bu_in_lam_w'][stage]]
                g_lam = self.nlp_sol_lam_g[self.index_map['lam_gnl_in_lam_g'][stage]]
            elif stage < dims.N:
                bx_lam = self.nlp_sol_lam_x[self.index_map['lam_bx_in_lam_w'][stage]]
                bu_lam = self.nlp_sol_lam_x[self.index_map['lam_bu_in_lam_w'][stage]]
                g_lam = self.nlp_sol_lam_g[self.index_map['lam_gnl_in_lam_g'][stage]]
            elif stage == dims.N:
                bx_lam = self.nlp_sol_lam_x[self.index_map['lam_bx_in_lam_w'][stage]]
                bu_lam = np.empty((0, 1))
                g_lam = self.nlp_sol_lam_g[self.index_map['lam_gnl_in_lam_g'][stage]]

            lbx_lam = np.maximum(0, -bx_lam)
            ubx_lam = np.maximum(0, bx_lam)
            lbu_lam = np.maximum(0, -bu_lam)
            ubu_lam = np.maximum(0, bu_lam)
            if any([dims.ns_0, dims.ns, dims.ns_e]):
                lw_soft_lam = self.nlp_sol_lam_x[self.index_map['sl_in_w'][stage]]
                uw_soft_lam = self.nlp_sol_lam_x[self.index_map['su_in_w'][stage]]
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

        if field == 'x':
            self.w0[self.index_map['x_in_w'][stage]] = value_.flatten()
        elif field == 'u':
            self.w0[self.index_map['u_in_w'][stage]] = value_.flatten()
        elif field == 'pi':
            self.lam_g0[self.index_map['pi_in_lam_g'][stage]] = -value_.flatten()
        elif field == 'p':
            self.p[self.index_map['p_in_p_nlp'][stage]] = value_.flatten()
        elif field == 'sl':
            self.w0[self.index_map['sl_in_w'][stage]] = value_.flatten()
        elif field == 'su':
            self.w0[self.index_map['su_in_w'][stage]] = value_.flatten()
        elif field == 'lam':
            if stage == 0:
                nbx = dims.nbx_0
                nbu = dims.nbu
                n_ghphi = dims.ng + dims.nh_0 + dims.nphi_0
                ns = dims.ns_0
            elif stage < dims.N:
                nbx = dims.nbx
                nbu = dims.nbu
                n_ghphi = dims.ng + dims.nh + dims.nphi
                ns = dims.ns
            elif stage == dims.N:
                nbx = dims.nbx_e
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
                self.lam_x0[self.index_map['lam_bx_in_lam_w'][stage]+self.index_map['lam_bu_in_lam_w'][stage]] = np.concatenate((ubx_lam-lbx_lam, ubu_lam-lbu_lam))
                self.lam_g0[self.index_map['lam_gnl_in_lam_g'][stage]] =  ug_lam_hard-lg_lam_hard
                self.lam_g0[self.index_map['lam_sl_in_lam_g'][stage]] = -lg_lam_soft
                self.lam_g0[self.index_map['lam_su_in_lam_g'][stage]] = ug_lam_soft
                self.lam_x0[self.index_map['sl_in_w'][stage]+self.index_map['su_in_w'][stage]] = -soft_lam
            else:
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
