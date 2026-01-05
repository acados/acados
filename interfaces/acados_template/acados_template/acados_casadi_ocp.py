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

import numpy as np

from .utils import casadi_length, is_casadi_SX, is_empty
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
            'yref_in_p_nlp': [],
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
        self.offset_p = 0

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
                self._append_node_for_variables('x', ca_symbol, xtraj_nodes, i, dims)
                self._append_node_for_variables('u', ca_symbol, utraj_nodes, i, dims)
                self._append_node_for_variables('sl', ca_symbol, sl_nodes, i, dims)
                self._append_node_for_variables('su', ca_symbol, su_nodes, i, dims)
        else: # single_shooting
            self._x_traj_fun = []
            for i in range(N_horizon):
                self._append_node_for_variables('u', ca_symbol, utraj_nodes, i, dims)
                self._append_node_for_variables('sl', ca_symbol, sl_nodes, i, dims)
                self._append_node_for_variables('su', ca_symbol, su_nodes, i, dims)

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
        x_guess = ocp.constraints.x0 if ocp.constraints.has_x0 else np.zeros((dims.nx,))
        if multiple_shooting:
            for i in range(N_horizon+1):
                self._append_variables_and_bounds('x', w_sym_list, lbw_list, ubw_list, w0_list, xtraj_nodes, lb_xtraj_nodes, ub_xtraj_nodes, i, dims, x_guess)
                if i < N_horizon:
                    self._append_variables_and_bounds('u',w_sym_list, lbw_list, ubw_list, w0_list, utraj_nodes, lb_utraj_nodes, ub_utraj_nodes, i, dims, x_guess)
                self._append_variables_and_bounds('slack', w_sym_list, lbw_list, ubw_list, w0_list, [sl_nodes, su_nodes], lb_slack_nodes, ub_slack_nodes, i, dims, x_guess)
        else: # single_shooting
            xtraj_nodes.append(x_guess)
            self._x_traj_fun.append(x_guess)
            for i in range(N_horizon+1):
                if i < N_horizon:
                    self._append_variables_and_bounds('u', w_sym_list, lbw_list, ubw_list, w0_list, utraj_nodes, lb_utraj_nodes, ub_utraj_nodes, i, dims, x_guess)
                self._append_variables_and_bounds('slack', w_sym_list, lbw_list, ubw_list, w0_list, [sl_nodes, su_nodes], lb_slack_nodes, ub_slack_nodes, i, dims, x_guess)

        nw = self.offset_w  # number of primal variables

        # setup parameter nodes and value
        ptraj_nodes = []
        p_list = []
        for i in range(N_horizon+1):
            self._append_node_and_value_for_params(ca_symbol, p_list, ptraj_nodes, i, ocp)
        p_list.append(ocp.p_global_values)
        self._index_map['p_global_in_p_nlp'].append(list(range(self.offset_p, self.offset_p+dims.np_global)))
        self.offset_p += dims.np_global

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
                    utraj_node = utraj_nodes[i] if dims.nu > 0 else ca_symbol('dummy_u', 0, 1)
                    ptraj_node = ptraj_nodes[i][:dims.np] if dims.np > 0 else ca_symbol('dummy_p', 0, 1)
                    if solver_options.integrator_type == "DISCRETE":
                        dyn_equality = xtraj_nodes[i+1] - f_discr_fun(xtraj_nodes[i], utraj_node, ptraj_node, model.p_global)
                    elif solver_options.integrator_type == "ERK":
                        param = ca.vertcat(utraj_node, ptraj_node, model.p_global)
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
                    utraj_node = utraj_nodes[i] if dims.nu > 0 else ca_symbol('dummy_u', 0, 1)
                    ptraj_node = ptraj_nodes[i][:dims.np] if dims.np > 0 else ca_symbol('dummy_p', 0, 1)
                    x_current = xtraj_nodes[i]
                    if solver_options.integrator_type == "DISCRETE":
                        x_next = f_discr_fun(x_current, utraj_node, ptraj_node, model.p_global)
                    elif solver_options.integrator_type == "ERK":
                        param = ca.vertcat(utraj_node, ptraj_node, model.p_global)
                        x_next = f_discr_fun(x_current, param, solver_options.time_steps[i])
                    xtraj_nodes.append(x_next)
                    self._x_traj_fun.append(f_discr_fun)

            # Nonlinear Constraints
            constraint_dict = self._get_constraint_node(i, N_horizon, xtraj_nodes, utraj_nodes, ptraj_nodes, model, constraints, dims)

            # add linear constraints
            if constraint_dict['ng'] > 0:
                self._append_constraints(i, 'gnl', g, lbg, ubg,
                                         g_expr = constraint_dict['linear_constr_expr'],
                                         lbg_expr = constraint_dict['lg'],
                                         ubg_expr = constraint_dict['ug'],
                                         cons_dim=constraint_dict['ng'])

            # add nonlinear constraints using constraint_dict directly (no locals)
            if constraint_dict['nh'] > 0:
                if constraint_dict['nsh'] > 0:
                    # h_fun with slack variables
                    soft_h_indices = constraint_dict['idxsh']
                    hard_h_indices = np.array([h for h in range(len(constraint_dict['lh'])) if h not in constraint_dict['idxsh']])
                    for index_in_nh in range(constraint_dict['nh']):
                        if index_in_nh in soft_h_indices:
                            index_in_soft = soft_h_indices.tolist().index(index_in_nh)
                            self._append_constraints(i, 'gnl', g, lbg, ubg,
                                                     g_expr = constraint_dict['h_i_nlp_expr'][index_in_nh] + sl_nodes[i][index_in_soft],
                                                     lbg_expr = constraint_dict['lh'][index_in_nh],
                                                     ubg_expr = np.inf * ca.DM.ones((1, 1)),
                                                     cons_dim=1,
                                                     sl=True)
                            self._append_constraints(i, 'gnl', g, lbg, ubg,
                                                     g_expr = constraint_dict['h_i_nlp_expr'][index_in_nh] - su_nodes[i][index_in_soft],
                                                     lbg_expr = -np.inf * ca.DM.ones((1, 1)),
                                                     ubg_expr = constraint_dict['uh'][index_in_nh],
                                                     cons_dim=1,
                                                     su=True)
                        elif index_in_nh in hard_h_indices:
                            self._append_constraints(i, 'gnl', g, lbg, ubg,
                                                     g_expr = constraint_dict['h_i_nlp_expr'][index_in_nh],
                                                     lbg_expr = constraint_dict['lh'][index_in_nh],
                                                     ubg_expr = constraint_dict['uh'][index_in_nh],
                                                     cons_dim=1)
                else:
                    self._append_constraints(i, 'gnl', g, lbg, ubg,
                                             g_expr = constraint_dict['h_i_nlp_expr'],
                                             lbg_expr = constraint_dict['lh'],
                                             ubg_expr = constraint_dict['uh'],
                                             cons_dim=constraint_dict['nh'])
                if with_hessian:
                    # add hessian contribution
                    lam_h = ca_symbol(f'lam_h_{i}', constraint_dict['nh'], 1)
                    lam_g.append(lam_h)
                    if ocp.solver_options.hessian_approx == 'EXACT' and ocp.solver_options.exact_hess_constr:
                        adj = ca.jtimes(constraint_dict['h_i_nlp_expr'], w, lam_h, True)
                        hess_l += ca.jacobian(adj, w, {"symmetric": is_casadi_SX(model.x)})

            # add convex-over-nonlinear constraints
            if constraint_dict['nphi'] > 0:
                conl_constr_fun = constraint_dict['conl_constr_fun']
                utraj_node = utraj_nodes[i] if dims.nu > 0 else ca_symbol('dummy_u', 0, 1)
                ptraj_node = ptraj_nodes[i][:dims.np] if dims.np > 0 else ca_symbol('dummy_p', 0, 1)
                self._append_constraints(i, 'gnl', g, lbg, ubg,
                                         g_expr = conl_constr_fun(xtraj_nodes[i], utraj_node, ptraj_node, model.p_global),
                                         lbg_expr = constraint_dict['lphi'],
                                         ubg_expr = constraint_dict['uphi'],
                                         cons_dim=constraint_dict['nphi'])
                if with_hessian:
                    lam_phi = ca_symbol(f'lam_phi_{i}', constraint_dict['nphi'], 1)
                    lam_g.append(lam_phi)
                    # always use CONL Hessian approximation here, disregarding inner second derivative
                    outer_hess_r = ca.vertcat(*[ca.hessian(constraint_dict['con_phi_expr'][j], constraint_dict['con_r_in_phi'])[0] for j in range(constraint_dict['nphi'])])
                    outer_hess_r = ca.substitute(outer_hess_r, constraint_dict['con_r_in_phi'], constraint_dict['con_r_expr'])
                    r_in_nlp = ca.substitute(constraint_dict['con_r_expr'], model.x, xtraj_nodes[-1])
                    dr_dw = ca.jacobian(r_in_nlp, w)
                    hess_l += dr_dw.T @ outer_hess_r @ dr_dw

        ### Cost and residual
        nlp_cost = 0
        residual_list = []
        for i in range(N_horizon+1):
            cost_dict = self._get_cost_node(i, N_horizon, xtraj_nodes, utraj_nodes, p_nlp, sl_nodes, su_nodes, ocp)

            cost_fun_i = ca.Function(f'cost_fun_{i}', [model.x, model.u, cost_dict['p_for_model'], model.p_global], [cost_dict['cost_expr']])
            cost_i = cost_fun_i(cost_dict['xtraj_node'], cost_dict['utraj_node'], cost_dict['ptraj_node'], model.p_global)
            nlp_cost += solver_options.cost_scaling[i] * cost_i

            if cost_dict['residual_expr'] is not None:
                residual_fun_i = ca.Function(f'residual_fun_{i}', [model.x, model.u, cost_dict['p_for_model'], model.p_global], [cost_dict['residual_expr']])
                residual_i = ca.sqrt(solver_options.cost_scaling[i]) * ca.sqrt(cost_dict['W_mat']) @ residual_fun_i(cost_dict['xtraj_node'], cost_dict['utraj_node'], cost_dict['ptraj_node'], model.p_global)
                residual_list.append(residual_i)
            if cost_dict['ns'] > 0:
                penalty_expr_i = 0.5 * ca.mtimes(cost_dict['sl_node'].T, ca.mtimes(np.diag(cost_dict['Zl']), cost_dict['sl_node'])) + \
                    ca.mtimes(cost_dict['zl'].reshape(-1, 1).T, cost_dict['sl_node']) + \
                    0.5 * ca.mtimes(cost_dict['su_node'].T, ca.mtimes(np.diag(cost_dict['Zu']), cost_dict['su_node'])) + \
                    ca.mtimes(cost_dict['zu'].reshape(-1, 1).T, cost_dict['su_node'])
                nlp_cost += solver_options.cost_scaling[i] * penalty_expr_i

        if with_hessian:
            lam_f = ca_symbol('lam_f', 1, 1)
            if (ocp.cost.cost_type == "LINEAR_LS" and ocp.cost.cost_type_0 == "LINEAR_LS" and ocp.cost.cost_type_e == "LINEAR_LS"):
                hess_l += lam_f * ca.hessian(nlp_cost, w)[0]
            elif (ocp.cost.cost_type == "NONLINEAR_LS" and ocp.cost.cost_type_0 == "NONLINEAR_LS" and ocp.cost.cost_type_e == "NONLINEAR_LS"):
                if  ocp.solver_options.hessian_approx == 'EXACT':
                    hess_l += lam_f * ca.hessian(nlp_cost, w)[0]
                elif ocp.solver_options.hessian_approx == 'GAUSS_NEWTON':
                    Gauss_Newton = 0
                    for i in range(len(residual_list)):
                        dr_dw = ca.jacobian(residual_list[i], w)
                        Gauss_Newton += dr_dw.T @ dr_dw
                    hess_l += lam_f * Gauss_Newton
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

    def _append_node_for_variables(self, _field, ca_symbol, node:list, i, dims):
        """
        Helper function to append a node to the NLP formulation.
        """
        if i == 0:
            ns = dims.ns_0
            nx = dims.nx
            nu = dims.nu
        elif i < dims.N:
            ns = dims.ns
            nx = dims.nx
            nu = dims.nu
        else:
            ns = dims.ns_e
            nx = dims.nx
            nu = 0

        if _field == 'x':
            node.append(ca_symbol(f'x{i}', nx, 1))
        elif _field == 'u':
            node.append(ca_symbol(f'u{i}', nu, 1))
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

    def _append_node_and_value_for_params(self, ca_symbol, p_list, node, stage, ocp:AcadosOcp):
        """
        Helper function to append parameter values to the NLP formulation.
        """
        if stage == 0:
            ny = ocp.dims.ny_0
            yref = ocp.cost.yref_0
        elif stage < ocp.solver_options.N_horizon:
            ny = ocp.dims.ny
            yref = ocp.cost.yref
        else:
            ny = ocp.dims.ny_e
            yref = ocp.cost.yref_e

        node.append(ca.vertcat(ca_symbol(f'p{stage}', ocp.dims.np, 1), ca_symbol(f'yref{stage}', ny, 1)))
        p_list.append(ocp.parameter_values)
        self._index_map['p_in_p_nlp'].append(list(range(self.offset_p, self.offset_p + ocp.dims.np)))
        self.offset_p += ocp.dims.np
        p_list.append(yref)
        self._index_map['yref_in_p_nlp'].append(list(range(self.offset_p, self.offset_p + ny)))
        self.offset_p += ny

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

    def _get_cost_node(self, i, N_horizon, xtraj_node, utraj_node, p_nlp, sl_node, su_node, ocp: AcadosOcp):
        """
        Helper function to get the cost node for a given stage.
        """
        cost_dict = {}
        dims = ocp.dims
        model = ocp.model
        cost = ocp.cost
        p_index = self._index_map['p_in_p_nlp'][i]
        yref_index = self._index_map['yref_in_p_nlp'][i]
        yref = p_nlp[yref_index]
        if i == 0:
            if cost.cost_type_0 == "NONLINEAR_LS" and not is_empty(ocp.model.cost_y_expr_0):
                y = ocp.model.cost_y_expr_0
                residual_expr = y - yref
            elif cost.cost_type_0 == "LINEAR_LS" and not is_empty(cost.Vx_0) and not is_empty(cost.Vu_0):
                y = cost.Vx_0 @ model.x + cost.Vu_0 @ model.u
                residual_expr = y - yref
            else:
                residual_expr = None
            cost_dict['xtraj_node'] = xtraj_node[0]
            cost_dict['utraj_node'] = utraj_node[0]
            cost_dict['ptraj_node'] = p_nlp[p_index + yref_index]
            cost_dict['sl_node'] = sl_node[0]
            cost_dict['su_node'] = su_node[0]
            cost_dict['cost_expr'] = ocp.get_initial_cost_expression(yref)
            cost_dict['residual_expr'] = residual_expr
            cost_dict['p_for_model'] = ca.vertcat(ocp.model.p, yref)
            cost_dict['ns'] = dims.ns_0
            cost_dict['W_mat'] = cost.W_0
            cost_dict['zl'] = cost.zl_0
            cost_dict['Zl'] = cost.Zl_0
            cost_dict['zu'] = cost.zu_0
            cost_dict['Zu'] = cost.Zu_0
        elif i < N_horizon:
            if cost.cost_type == "NONLINEAR_LS" and not is_empty(ocp.model.cost_y_expr):
                y = ocp.model.cost_y_expr
                residual_expr = y - yref
            elif cost.cost_type == "LINEAR_LS" and not is_empty(cost.Vx) and not is_empty(cost.Vu):
                y = cost.Vx @ model.x + cost.Vu @ model.u
                residual_expr = y - yref
            else:
                residual_expr = None
            cost_dict['xtraj_node'] = xtraj_node[i]
            cost_dict['utraj_node'] = utraj_node[i]
            cost_dict['ptraj_node'] = p_nlp[p_index + yref_index]
            cost_dict['sl_node'] = sl_node[i]
            cost_dict['su_node'] = su_node[i]
            cost_dict['cost_expr'] = ocp.get_path_cost_expression(yref)
            cost_dict['residual_expr'] = residual_expr
            cost_dict['p_for_model'] = ca.vertcat(ocp.model.p, yref)
            cost_dict['ns'] = dims.ns
            cost_dict['W_mat'] = cost.W
            cost_dict['zl'] = cost.zl
            cost_dict['Zl'] = cost.Zl
            cost_dict['zu'] = cost.zu
            cost_dict['Zu'] = cost.Zu
        else:
            if cost.cost_type_e == "NONLINEAR_LS" and not is_empty(ocp.model.cost_y_expr_e):
                y = ocp.model.cost_y_expr_e
                residual_expr = y - yref
            elif cost.cost_type_e == "LINEAR_LS" and not is_empty(cost.Vx_e):
                y = cost.Vx_e @ model.x
                residual_expr = y - yref
            else:
                residual_expr = None
            cost_dict['xtraj_node'] = xtraj_node[-1]
            cost_dict['utraj_node'] = utraj_node[-1]
            cost_dict['ptraj_node'] = p_nlp[p_index + yref_index]
            cost_dict['sl_node'] = sl_node[-1]
            cost_dict['su_node'] = su_node[-1]
            cost_dict['cost_expr'] = ocp.get_terminal_cost_expression(yref)
            cost_dict['residual_expr'] = residual_expr
            cost_dict['p_for_model'] = ca.vertcat(ocp.model.p, yref)
            cost_dict['ns'] = dims.ns_e
            cost_dict['W_mat'] = cost.W_e
            cost_dict['zl'] = cost.zl_e
            cost_dict['Zl'] = cost.Zl_e
            cost_dict['zu'] = cost.zu_e
            cost_dict['Zu'] = cost.Zu_e

        return cost_dict

    def _get_constraint_node(self, i, N_horizon, xtraj_node, utraj_node, ptraj_node, model, constraints, dims):
        """
        Helper function to get the constraint node for a given stage.
        """
        # dict to store constraint function for each node
        cons_dict = {}
        if i == 0 and N_horizon > 0:
            cons_dict['lg'], cons_dict['ug'] = constraints.lg, constraints.ug
            cons_dict['lh'], cons_dict['uh'] = constraints.lh_0, constraints.uh_0
            cons_dict['lphi'], cons_dict['uphi'] = constraints.lphi_0, constraints.uphi_0
            cons_dict['ng'], cons_dict['nh'], cons_dict['nphi'] = dims.ng, dims.nh_0, dims.nphi_0
            cons_dict['nsg'], cons_dict['nsh'], cons_dict['nsphi'], cons_dict['idxsh'] = dims.nsg, dims.nsh_0, dims.nsphi_0, constraints.idxsh_0

            # linear function
            if dims.ng > 0:
                cons_dict['C'] = constraints.C
                cons_dict['D'] = constraints.D
                cons_dict['linear_constr_expr'] = ca.mtimes(cons_dict['C'], xtraj_node[i]) + ca.mtimes(cons_dict['D'], utraj_node[i])
            # nonlinear function
            cons_dict['h_fun'] = ca.Function('h_0_fun', [model.x, model.u, model.p, model.p_global], [model.con_h_expr_0])
            cons_dict['h_i_nlp_expr'] = cons_dict['h_fun'](xtraj_node[i], utraj_node[i], ptraj_node[i][:dims.np], model.p_global)
            # compound nonlinear constraint
            cons_dict['conl_constr_fun'] = None
            if dims.nphi_0 > 0:
                cons_dict['con_phi_expr'], cons_dict['con_r_in_phi'], cons_dict['con_r_expr'] = model.con_phi_expr_0, model.con_r_in_phi_0, model.con_r_expr_0
                conl_expr = ca.substitute(model.con_phi_expr_0, model.con_r_in_phi_0, model.con_r_expr_0)
                cons_dict['conl_constr_fun'] = ca.Function('conl_constr_0_fun', [model.x, model.u, model.p, model.p_global], [conl_expr])

        elif i < N_horizon:
            # populate cons_dict for intermediate stage (mirror initial-stage style)
            cons_dict['lg'], cons_dict['ug'] = constraints.lg, constraints.ug
            cons_dict['lh'], cons_dict['uh'] = constraints.lh, constraints.uh
            cons_dict['lphi'], cons_dict['uphi'] = constraints.lphi, constraints.uphi
            cons_dict['ng'], cons_dict['nh'], cons_dict['nphi'] = dims.ng, dims.nh, dims.nphi
            cons_dict['nsg'], cons_dict['nsh'], cons_dict['nsphi'], cons_dict['idxsh'] = dims.nsg, dims.nsh, dims.nsphi, constraints.idxsh

            # linear function
            if dims.ng > 0:
                cons_dict['C'] = constraints.C
                cons_dict['D'] = constraints.D
                cons_dict['linear_constr_expr'] = ca.mtimes(cons_dict['C'], xtraj_node[i]) + ca.mtimes(cons_dict['D'], utraj_node[i])
            # nonlinear function
            cons_dict['h_fun'] = ca.Function('h_fun', [model.x, model.u, model.p, model.p_global], [model.con_h_expr])
            cons_dict['h_i_nlp_expr'] = cons_dict['h_fun'](xtraj_node[i], utraj_node[i], ptraj_node[i][:dims.np], model.p_global)
            # compound nonlinear constraint
            cons_dict['conl_constr_fun'] = None
            if dims.nphi > 0:
                cons_dict['con_phi_expr'], cons_dict['con_r_in_phi'], cons_dict['con_r_expr'] = model.con_phi_expr, model.con_r_in_phi, model.con_r_expr
                conl_expr = ca.substitute(model.con_phi_expr, model.con_r_in_phi, model.con_r_expr)
                cons_dict['conl_constr_fun'] = ca.Function('conl_constr_fun', [model.x, model.u, model.p, model.p_global], [conl_expr])

        else:
            # populate cons_dict for terminal stage
            cons_dict['lg'], cons_dict['ug'] = constraints.lg_e, constraints.ug_e
            cons_dict['lh'], cons_dict['uh'] = constraints.lh_e, constraints.uh_e
            cons_dict['lphi'], cons_dict['uphi'] = constraints.lphi_e, constraints.uphi_e
            cons_dict['ng'], cons_dict['nh'], cons_dict['nphi'] = dims.ng_e, dims.nh_e, dims.nphi_e
            cons_dict['nsg'], cons_dict['nsh'], cons_dict['nsphi'], cons_dict['idxsh'] = dims.nsg_e, dims.nsh_e, dims.nsphi_e, constraints.idxsh_e

            # linear function
            if dims.ng_e > 0:
                cons_dict['C'] = constraints.C_e
                cons_dict['linear_constr_expr'] = ca.mtimes(cons_dict['C'], xtraj_node[i])
            # nonlinear function
            cons_dict['h_fun'] = ca.Function('h_e_fun', [model.x, model.p, model.p_global], [model.con_h_expr_e])
            cons_dict['h_i_nlp_expr'] = cons_dict['h_fun'](xtraj_node[i], ptraj_node[i][:dims.np], model.p_global)
            # compound nonlinear constraint
            cons_dict['conl_constr_fun'] = None
            if dims.nphi_e > 0:
                cons_dict['con_phi_expr'], cons_dict['con_r_in_phi'], cons_dict['con_r_expr'] = model.con_phi_expr_e, model.con_r_in_phi_e, model.con_r_expr_e
                conl_expr = ca.substitute(model.con_phi_expr_e, model.con_r_in_phi_e, model.con_r_expr_e)
                cons_dict['conl_constr_fun'] = ca.Function('conl_constr_e_fun', [model.x, model.u, model.p, model.p_global], [conl_expr])

        self._index_map['lam_gnl_in_lam_g'].append([])
        self._index_map['lam_sl_in_lam_g'].append([])
        self._index_map['lam_su_in_lam_g'].append([])

        return cons_dict

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