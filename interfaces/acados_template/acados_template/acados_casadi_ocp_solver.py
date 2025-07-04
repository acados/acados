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
        index_map = {
            # indices of variables within w
            'x_in_w': [],
            'u_in_w': [],
            # indices of parameters within p_nlp
            'p_in_p_nlp': [],
            'p_global_in_p_nlp': [],
            # indices of state bounds within lam_x(lam_w) in casadi formulation
            'lam_bx_in_lam_w':[],
            # indices of control bounds within lam_x(lam_w) in casadi formulation
            'lam_bu_in_lam_w': [],
            # indices of dynamic constraints within g in casadi formulation
            'pi_in_lam_g': [],
            # indicies to [g, h, phi] in acados formulation within lam_g in casadi formulation
            'lam_gnl_in_lam_g': []
        }

        # unpack
        model = ocp.model
        dims = ocp.dims
        constraints = ocp.constraints
        solver_options = ocp.solver_options
        N_horizon = solver_options.N_horizon

        # check what is not supported yet
        if any([dims.ns_0, dims.ns, dims.ns_e]):
            raise NotImplementedError("AcadosCasadiOcpSolver does not support soft constraints yet.")
        if dims.nz > 0:
            raise NotImplementedError("AcadosCasadiOcpSolver does not support algebraic variables (z) yet.")
        if ocp.solver_options.integrator_type not in ["DISCRETE", "ERK"]:
            raise NotImplementedError(f"AcadosCasadiOcpSolver does not support integrator type {ocp.solver_options.integrator_type} yet.")

        # create primal variables indexed by shooting nodes
        ca_symbol = model.get_casadi_symbol()
        xtraj_node = [ca_symbol(f'x{i}', dims.nx, 1) for i in range(N_horizon+1)]
        utraj_node = [ca_symbol(f'u{i}', dims.nu, 1) for i in range(N_horizon)]
        if dims.nz > 0:
            raise NotImplementedError("CasADi NLP formulation not implemented for models with algebraic variables (z).")

        # parameters
        ptraj_node = [ca_symbol(f'p{i}', dims.np, 1) for i in range(N_horizon+1)]

        # setup state and control bounds
        lb_xtraj_node = [-np.inf * ca.DM.ones((dims.nx, 1)) for _ in range(N_horizon+1)]
        ub_xtraj_node = [np.inf * ca.DM.ones((dims.nx, 1)) for _ in range(N_horizon+1)]
        lb_utraj_node = [-np.inf * ca.DM.ones((dims.nu, 1)) for _ in range(N_horizon)]
        ub_utraj_node = [np.inf * ca.DM.ones((dims.nu, 1)) for _ in range(N_horizon)]
        offset = 0
        for i in range(0, N_horizon+1):
            if i == 0:
                lb_xtraj_node[i][constraints.idxbx_0] = constraints.lbx_0
                ub_xtraj_node[i][constraints.idxbx_0] = constraints.ubx_0
                index_map['lam_bx_in_lam_w'].append(list(offset + constraints.idxbx_0))
                offset += dims.nx
            elif i < N_horizon:
                lb_xtraj_node[i][constraints.idxbx] = constraints.lbx
                ub_xtraj_node[i][constraints.idxbx] = constraints.ubx
                index_map['lam_bx_in_lam_w'].append(list(offset + constraints.idxbx))
                offset += dims.nx
            elif i == N_horizon:
                lb_xtraj_node[-1][constraints.idxbx_e] = constraints.lbx_e
                ub_xtraj_node[-1][constraints.idxbx_e] = constraints.ubx_e
                index_map['lam_bx_in_lam_w'].append(list(offset + constraints.idxbx_e))
                offset += dims.nx
            if i < N_horizon:
                lb_utraj_node[i][constraints.idxbu] = constraints.lbu
                ub_utraj_node[i][constraints.idxbu] = constraints.ubu
                index_map['lam_bu_in_lam_w'].append(list(offset + constraints.idxbu))
                offset += dims.nu

        ### Concatenate primal variables and bounds
        # w = [x0, u0, x1, u1, ...]
        w_sym_list = []
        lbw_list = []
        ubw_list = []
        w0_list = []
        p_list = []
        offset = 0
        offset_p = 0
        x_guess = ocp.constraints.x0 if ocp.constraints.has_x0 else np.zeros((dims.nx,))
        for i in range(N_horizon):
            # add x
            w_sym_list.append(xtraj_node[i])
            lbw_list.append(lb_xtraj_node[i])
            ubw_list.append(ub_xtraj_node[i])
            w0_list.append(x_guess)
            index_map['x_in_w'].append(list(range(offset, offset + dims.nx)))
            offset += dims.nx
            # add u
            w_sym_list.append(utraj_node[i])
            lbw_list.append(lb_utraj_node[i])
            ubw_list.append(ub_utraj_node[i])
            w0_list.append(np.zeros((dims.nu,)))
            index_map['u_in_w'].append(list(range(offset, offset + dims.nu)))
            offset += dims.nu
            # add parameters
            p_list.append(ocp.parameter_values)
            index_map['p_in_p_nlp'].append(list(range(offset_p, offset_p+dims.np)))
            offset_p += dims.np
        ## terminal stage
        # add x
        w_sym_list.append(xtraj_node[-1])
        lbw_list.append(lb_xtraj_node[-1])
        ubw_list.append(ub_xtraj_node[-1])
        w0_list.append(x_guess)
        index_map['x_in_w'].append(list(range(offset, offset + dims.nx)))
        offset += dims.nx
        # add parameters
        p_list.append(ocp.parameter_values)
        index_map['p_in_p_nlp'].append(list(range(offset_p, offset_p+dims.np)))
        offset_p += dims.np
        p_list.append(ocp.p_global_values)
        index_map['p_global_in_p_nlp'].append(list(range(offset_p, offset_p+dims.np_global)))
        offset_p += dims.np_global

        nw = offset  # number of primal variables

        # vectorize
        w = ca.vertcat(*w_sym_list)
        lbw = ca.vertcat(*lbw_list)
        ubw = ca.vertcat(*ubw_list)
        p_nlp = ca.vertcat(*ptraj_node, model.p_global)

        ### Nonlinear constraints
        # dynamics
        if solver_options.integrator_type == "DISCRETE":
            f_discr_fun = ca.Function('f_discr_fun', [model.x, model.u, model.p, model.p_global], [model.disc_dyn_expr])
        elif solver_options.integrator_type == "ERK":
            para = ca.vertcat(model.u, model.p, model.p_global)
            ca_expl_ode = ca.Function('ca_expl_ode', [model.x, para], [model.f_expl_expr])
            f_discr_fun = ca.simpleRK(ca_expl_ode, solver_options.sim_method_num_steps[0], solver_options.sim_method_num_stages[0])
        else:
            raise NotImplementedError(f"Integrator type {solver_options.integrator_type} not supported.")

        # initial
        h_0_fun = ca.Function('h_0_fun', [model.x, model.u, model.p, model.p_global], [model.con_h_expr_0])

        # intermediate
        h_fun = ca.Function('h_fun', [model.x, model.u, model.p, model.p_global], [model.con_h_expr])
        if dims.nphi > 0:
            conl_constr_expr = ca.substitute(model.con_phi_expr, model.con_r_in_phi, model.con_r_expr)
            conl_constr_fun = ca.Function('conl_constr_fun', [model.x, model.u, model.p, model.p_global], [conl_constr_expr])

        # terminal
        h_e_fun = ca.Function('h_e_fun', [model.x, model.p, model.p_global], [model.con_h_expr_e])
        if dims.nphi_e > 0:
            conl_constr_expr_e = ca.substitute(model.con_phi_expr_e, model.con_r_in_phi_e, model.con_r_expr_e)
            conl_constr_e_fun = ca.Function('conl_constr_e_fun', [model.x, model.p, model.p_global], [conl_constr_expr_e])

        # create nonlinear constraints
        g = []
        lbg = []
        ubg = []
        offset = 0
        if with_hessian:
            lam_g = []
            hess_l = ca.DM.zeros((nw, nw))

        for i in range(N_horizon+1):
            # add dynamics constraints
            if i < N_horizon:
                if solver_options.integrator_type == "DISCRETE":
                    dyn_equality = xtraj_node[i+1] - f_discr_fun(xtraj_node[i], utraj_node[i], ptraj_node[i], model.p_global)
                elif solver_options.integrator_type == "ERK":
                    para = ca.vertcat(utraj_node[i], ptraj_node[i], model.p_global)
                    dyn_equality = xtraj_node[i+1] - f_discr_fun(xtraj_node[i], para, solver_options.time_steps[i])
                g.append(dyn_equality)
                lbg.append(np.zeros((dims.nx, 1)))
                ubg.append(np.zeros((dims.nx, 1)))
                index_map['pi_in_lam_g'].append(list(range(offset, offset+dims.nx)))
                offset += dims.nx

                if with_hessian:
                    # add hessian of dynamics constraints
                    lam_g_dyn = ca_symbol(f'lam_g_dyn{i}', dims.nx, 1)
                    lam_g.append(lam_g_dyn)
                    if ocp.solver_options.hessian_approx == 'EXACT' and ocp.solver_options.exact_hess_dyn:
                        adj = ca.jtimes(dyn_equality, w, lam_g_dyn, True)
                        hess_l += ca.jacobian(adj, w, {"symmetric": is_casadi_SX(model.x)})

            # nonlinear constraints
            # initial stage
            if i == 0 and N_horizon > 0:
                if dims.ng > 0:
                    C = constraints.C
                    D = constraints.D
                    linear_constr_expr = ca.mtimes(C, xtraj_node[0]) + ca.mtimes(D, utraj_node[0])
                    g.append(linear_constr_expr)
                    lbg.append(constraints.lg)
                    ubg.append(constraints.ug)

                if dims.nh_0 > 0:
                    # h_0
                    h_0_nlp_expr = h_0_fun(xtraj_node[0], utraj_node[0], ptraj_node[0], model.p_global)
                    g.append(h_0_nlp_expr)
                    lbg.append(constraints.lh_0)
                    ubg.append(constraints.uh_0)
                    if with_hessian:
                        lam_h_0 = ca_symbol(f'lam_h_0', dims.nh_0, 1)
                        lam_g.append(lam_h_0)
                        # add hessian contribution
                        if ocp.solver_options.hessian_approx == 'EXACT' and ocp.solver_options.exact_hess_constr:
                            adj = ca.jtimes(h_0_nlp_expr, w, lam_h_0, True)
                            hess_l += ca.jacobian(adj, w, {"symmetric": is_casadi_SX(model.x)})

                if dims.nphi_0 > 0:
                    conl_constr_expr_0 = ca.substitute(model.con_phi_expr_0, model.con_r_in_phi_0, model.con_r_expr_0)
                    conl_constr_0_fun = ca.Function('conl_constr_0_fun', [model.x, model.u, model.p, model.p_global], [conl_constr_expr_0])
                    g.append(conl_constr_0_fun(xtraj_node[0], utraj_node[0], ptraj_node[0], model.p_global))
                    lbg.append(constraints.lphi_0)
                    ubg.append(constraints.uphi_0)
                    if with_hessian:
                        lam_phi_0 = ca_symbol(f'lam_phi_0', dims.nphi_0, 1)
                        lam_g.append(lam_phi_0)
                        # always use CONL Hessian approximation here, disregarding inner second derivative
                        outer_hess_r = ca.vertcat(*[ca.hessian(model.con_phi_expr_0[i], model.con_r_in_phi_0)[0] for i in range(dims.nphi_0)])
                        outer_hess_r = ca.substitute(outer_hess_r, model.con_r_in_phi_0, model.con_r_expr_0)
                        r_in_nlp = ca.substitute(model.con_r_expr_0, model.x, xtraj_node[-1])
                        dr_dw = ca.jacobian(r_in_nlp, w)
                        hess_l += dr_dw.T @ outer_hess_r @ dr_dw

                index_map['lam_gnl_in_lam_g'].append(list(range(offset, offset + dims.ng + dims.nh_0 + dims.nphi_0)))
                offset += dims.ng + dims.nh_0 + dims.nphi_0
            
            # intermediate stages
            elif i < N_horizon:
                if dims.ng > 0:
                    C = constraints.C
                    D = constraints.D
                    linear_constr_expr = ca.mtimes(C, xtraj_node[i]) + ca.mtimes(D, utraj_node[i])
                    g.append(linear_constr_expr)
                    lbg.append(constraints.lg)
                    ubg.append(constraints.ug)

                if dims.nh > 0:
                    h_i_nlp_expr = h_fun(xtraj_node[i], utraj_node[i], ptraj_node[i], model.p_global)
                    g.append(h_i_nlp_expr)
                    lbg.append(constraints.lh)
                    ubg.append(constraints.uh)
                    if with_hessian and dims.nh > 0:
                        # add hessian contribution
                        lam_h = ca_symbol(f'lam_h_{i}', dims.nh, 1)
                        lam_g.append(lam_h)
                        if ocp.solver_options.hessian_approx == 'EXACT' and ocp.solver_options.exact_hess_constr:
                            adj = ca.jtimes(h_i_nlp_expr, w, lam_h, True)
                            hess_l += ca.jacobian(adj, w, {"symmetric": is_casadi_SX(model.x)})

                if dims.nphi > 0:
                    g.append(conl_constr_fun(xtraj_node[i], utraj_node[i], ptraj_node[i], model.p_global))
                    lbg.append(constraints.lphi)
                    ubg.append(constraints.uphi)
                    if with_hessian:
                        lam_phi = ca_symbol(f'lam_phi', dims.nphi, 1)
                        lam_g.append(lam_phi)
                        # always use CONL Hessian approximation here, disregarding inner second derivative
                        outer_hess_r = ca.vertcat(*[ca.hessian(model.con_phi_expr[i], model.con_r_in_phi)[0] for i in range(dims.nphi)])
                        outer_hess_r = ca.substitute(outer_hess_r, model.con_r_in_phi, model.con_r_expr)
                        r_in_nlp = ca.substitute(model.con_r_expr, model.x, xtraj_node[-1])
                        dr_dw = ca.jacobian(r_in_nlp, w)
                        hess_l += dr_dw.T @ outer_hess_r @ dr_dw

                index_map['lam_gnl_in_lam_g'].append(list(range(offset, offset + dims.ng + dims.nh + dims.nphi)))
                offset += dims.ng + dims.nphi + dims.nh
            
            # terminal stage
            else:
                if dims.ng_e > 0:
                    C_e = constraints.C_e
                    linear_constr_expr_e = ca.mtimes(C_e, xtraj_node[-1])
                    g.append(linear_constr_expr_e)
                    lbg.append(constraints.lg_e)
                    ubg.append(constraints.ug_e)

                if dims.nh_e > 0:
                    h_e_nlp_expr = h_e_fun(xtraj_node[-1], ptraj_node[-1], model.p_global)
                    g.append(h_e_nlp_expr)
                    lbg.append(constraints.lh_e)
                    ubg.append(constraints.uh_e)
                    if with_hessian and dims.nh_e > 0:
                        # add hessian contribution
                        lam_h_e = ca_symbol(f'lam_h_e', dims.nh_e, 1)
                        lam_g.append(lam_h_e)
                        if ocp.solver_options.hessian_approx == 'EXACT' and ocp.solver_options.exact_hess_constr:
                            adj = ca.jtimes(h_e_nlp_expr, w, lam_h_e, True)
                            hess_l += ca.jacobian(adj, w, {"symmetric": is_casadi_SX(model.x)})

                if dims.nphi_e > 0:
                    g.append(conl_constr_e_fun(xtraj_node[-1], ptraj_node[-1], model.p_global))
                    lbg.append(constraints.lphi_e)
                    ubg.append(constraints.uphi_e)
                    if with_hessian:
                        lam_phi_e = ca_symbol(f'lam_phi_e', dims.nphi_e, 1)
                        lam_g.append(lam_phi_e)
                        # always use CONL Hessian approximation here, disregarding inner second derivative
                        outer_hess_r = ca.vertcat(*[ca.hessian(model.con_phi_expr_e[i], model.con_r_in_phi_e)[0] for i in range(dims.nphi_e)])
                        outer_hess_r = ca.substitute(outer_hess_r, model.con_r_in_phi_e, model.con_r_expr_e)
                        r_in_nlp = ca.substitute(model.con_r_expr_e, model.x, xtraj_node[-1])
                        dr_dw = ca.jacobian(r_in_nlp, w)
                        hess_l += dr_dw.T @ outer_hess_r @ dr_dw

                index_map['lam_gnl_in_lam_g'].append(list(range(offset, offset + dims.ng_e + dims.nh_e + dims.nphi_e)))
                offset += dims.ng_e + dims.nh_e + dims.nphi_e

        ### Cost
        # initial cost term
        nlp_cost = 0
        cost_expr_0 = ocp.get_initial_cost_expression()
        cost_fun_0 = ca.Function('cost_fun_0', [model.x, model.u, model.p, model.p_global], [cost_expr_0])
        nlp_cost += solver_options.cost_scaling[0] * cost_fun_0(xtraj_node[0], utraj_node[0], ptraj_node[0], model.p_global)

        # intermediate cost term
        cost_expr = ocp.get_path_cost_expression()
        cost_fun = ca.Function('cost_fun', [model.x, model.u, model.p, model.p_global], [cost_expr])
        for i in range(1, N_horizon):
            nlp_cost += solver_options.cost_scaling[i] * cost_fun(xtraj_node[i], utraj_node[i], ptraj_node[i], model.p_global)

        # terminal cost term
        cost_expr_e = ocp.get_terminal_cost_expression()
        cost_fun_e = ca.Function('cost_fun_e', [model.x, model.p, model.p_global], [cost_expr_e])
        nlp_cost += solver_options.cost_scaling[-1] * cost_fun_e(xtraj_node[-1], ptraj_node[-1], model.p_global)

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

        # sanity check

        # create NLP
        nlp = {"x": w, "p": p_nlp, "g": ca.vertcat(*g), "f": nlp_cost}
        bounds = {"lbx": lbw, "ubx": ubw, "lbg": ca.vertcat(*lbg), "ubg": ca.vertcat(*ubg)}
        w0 = np.concatenate(w0_list)
        p = np.concatenate(p_list)

        self.__nlp = nlp
        self.__bounds = bounds
        self.__w0 = w0
        self.__p = p
        self.__index_map = index_map
        self.__nlp_hess_l_custom = nlp_hess_l_custom
        self.__hess_approx_expr = hess_l

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
        :param field: string in ['x', 'u', 'pi', 'p', 'lam']

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
            lbu_lam = np.maximum(0, -bu_lam)
            lg_lam = np.maximum(0, -g_lam)
            ubx_lam = np.maximum(0, bx_lam)
            ubu_lam = np.maximum(0, bu_lam)
            ug_lam = np.maximum(0, g_lam)
            lam = np.concatenate((lbu_lam, lbx_lam, lg_lam, ubu_lam, ubx_lam, ug_lam))
            return lam.flatten()
        elif field in ['sl', 'su', 'z']:
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
        # TODO: add slacks
        for key, traj in iterate.__dict__.items():
            field = key.replace('_traj', '')

            for n, val in enumerate(traj):
                if field in ['x', 'u', 'pi', 'lam']:
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

    def get_stats(self, field_: str) -> Union[int, float, np.ndarray]:

        if field_ == "nlp_iter":
            return self.nlp_iter
        else:
            raise NotImplementedError()

    def get_cost(self) -> float:
        return self.nlp_sol['f'].full().item()

    def set(self, stage: int, field: str, value_: np.ndarray):
        """
        Set solver initialization to stages.

        :param stage: integer corresponding to shooting node
        :param field: string in ['x', 'u', 'pi', 'lam']
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
        elif field == 'lam':
            if stage == 0:
                nbx = dims.nbx_0
                nbu = dims.nbu
                n_ghphi = dims.ng + dims.nh_0 + dims.nphi_0
            elif stage < dims.N:
                nbx = dims.nbx
                nbu = dims.nbu
                n_ghphi = dims.ng + dims.nh + dims.nphi
            elif stage == dims.N:
                nbx = dims.nbx_e
                nbu = 0
                n_ghphi = dims.ng_e + dims.nh_e + dims.nphi_e

            offset_u = (nbx+nbu+n_ghphi)
            lbu_lam = value_[:nbu]
            lbx_lam = value_[nbu:nbu+nbx]
            lg_lam = value_[nbu+nbx:nbu+nbx+n_ghphi]
            ubu_lam = value_[offset_u:offset_u+nbu]
            ubx_lam = value_[offset_u+nbu:offset_u+nbu+nbx]
            ug_lam = value_[offset_u+nbu+nbx:offset_u+nbu+nbx+n_ghphi]
            if stage != dims.N:
                self.lam_x0[self.index_map['lam_bx_in_lam_w'][stage]+self.index_map['lam_bu_in_lam_w'][stage]] = np.concatenate((ubx_lam-lbx_lam, ubu_lam-lbu_lam))
                self.lam_g0[self.index_map['lam_gnl_in_lam_g'][stage]] =  ug_lam-lg_lam
            else:
                self.lam_x0[self.index_map['lam_bx_in_lam_w'][stage]] = ubx_lam-lbx_lam
                self.lam_g0[self.index_map['lam_gnl_in_lam_g'][stage]] = ug_lam-lg_lam
        elif field in ['sl', 'su']:
            # do nothing for now, only empty is supported
            pass
        else:
            raise NotImplementedError(f"Field '{field}' is not yet implemented in set().")

    def cost_get(self, stage_: int, field_: str) -> np.ndarray:
        raise NotImplementedError()

    def cost_set(self, stage_: int, field_: str, value_):
        raise NotImplementedError()
