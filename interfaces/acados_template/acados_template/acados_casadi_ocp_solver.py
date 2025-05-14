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
    def create_casadi_nlp_formulation(cls, ocp: AcadosOcp) -> Tuple[dict, dict]:
        """
        Creates an equivalent CasADi NLP formulation of the OCP.
        Experimental, not fully implemented yet.

        :return: nlp_dict, bounds_dict
        """
        ocp.make_consistent()

        # unpack
        model = ocp.model
        dims = ocp.dims
        constraints = ocp.constraints
        solver_options = ocp.solver_options

        if any([dims.ns_0, dims.ns, dims.ns_e]):
            raise NotImplementedError("CasADi NLP formulation not implemented for formulations with soft constraints yet.")

        # create variables
        ca_symbol = model.get_casadi_symbol()
        xtraj = ca_symbol('x', dims.nx, solver_options.N_horizon+1)
        utraj = ca_symbol('u', dims.nu, solver_options.N_horizon)
        if dims.nz > 0:
            raise NotImplementedError("CasADi NLP formulation not implemented for models with algebraic variables (z).")
        # parameters
        ptraj = ca_symbol('p', dims.np, solver_options.N_horizon+1)

        ### Constraints: bounds
        # setup state bounds
        lb_xtraj = -np.inf * ca.DM.ones((dims.nx, solver_options.N_horizon+1))
        ub_xtraj = np.inf * ca.DM.ones((dims.nx, solver_options.N_horizon+1))
        lb_xtraj[constraints.idxbx_0, 0] = constraints.lbx_0
        ub_xtraj[constraints.idxbx_0, 0] = constraints.ubx_0
        for i in range(1, solver_options.N_horizon):
            lb_xtraj[constraints.idxbx, i] = constraints.lbx
            ub_xtraj[constraints.idxbx, i] = constraints.ubx
        lb_xtraj[constraints.idxbx_e, -1] = constraints.lbx_e
        ub_xtraj[constraints.idxbx_e, -1] = constraints.ubx_e
        # setup control bounds
        lb_utraj = -np.inf * ca.DM.ones((dims.nu, solver_options.N_horizon))
        ub_utraj = np.inf * ca.DM.ones((dims.nu, solver_options.N_horizon))
        for i in range(solver_options.N_horizon):
            lb_utraj[constraints.idxbu, i] = constraints.lbu
            ub_utraj[constraints.idxbu, i] = constraints.ubu

        ### Nonlinear constraints
        g = []
        lbg = []
        ubg = []
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
                g.append(xtraj[:, i+1] - f_discr_fun(xtraj[:, i], utraj[:, i], ptraj[:, i], model.p_global))
            elif solver_options.integrator_type == "ERK":
                g.append(xtraj[:, i+1] - f_discr_fun(xtraj[:, i], utraj[:, i], solver_options.time_steps[i]))
            lbg.append(np.zeros((dims.nx, 1)))
            ubg.append(np.zeros((dims.nx, 1)))

        # nonlinear constraints -- initial stage
        h0_fun = ca.Function('h0_fun', [model.x, model.u, model.p, model.p_global], [model.con_h_expr_0])
        g.append(h0_fun(xtraj[:, 0], utraj[:, 0], ptraj[:, 0], model.p_global))
        lbg.append(constraints.lh_0)
        ubg.append(constraints.uh_0)

        if dims.nphi_0 > 0:
            conl_constr_expr_0 = ca.substitute(model.con_phi_expr_0, model.con_r_in_phi_0, model.con_r_expr_0)
            conl_constr_0_fun = ca.Function('conl_constr_0_fun', [model.x, model.u, model.p, model.p_global], [conl_constr_expr_0])
            g.append(conl_constr_0_fun(xtraj[:, 0], utraj[:, 0], ptraj[:, 0], model.p_global))
            lbg.append(constraints.lphi_0)
            ubg.append(constraints.uphi_0)

        # nonlinear constraints -- intermediate stages
        h_fun = ca.Function('h_fun', [model.x, model.u, model.p, model.p_global], [model.con_h_expr])

        if dims.nphi > 0:
            conl_constr_expr = ca.substitute(model.con_phi_expr, model.con_r_in_phi, model.con_r_expr)
            conl_constr_fun = ca.Function('conl_constr_fun', [model.x, model.u, model.p, model.p_global], [conl_constr_expr])

        for i in range(1, solver_options.N_horizon):
            g.append(h_fun(xtraj[:, i], utraj[:, i], ptraj[:, i], model.p_global))
            lbg.append(constraints.lh)
            ubg.append(constraints.uh)

            if dims.nphi > 0:
                g.append(conl_constr_fun(xtraj[:, i], utraj[:, i], ptraj[:, i], model.p_global))
                lbg.append(constraints.lphi)
                ubg.append(constraints.uphi)

        # nonlinear constraints -- terminal stage
        h_e_fun = ca.Function('h_e_fun', [model.x, model.p, model.p_global], [model.con_h_expr_e])

        g.append(h_e_fun(xtraj[:, -1], ptraj[:, -1], model.p_global))
        lbg.append(constraints.lh_e)
        ubg.append(constraints.uh_e)

        if dims.nphi_e > 0:
            conl_constr_expr_e = ca.substitute(model.con_phi_expr_e, model.con_r_in_phi_e, model.con_r_expr_e)
            conl_constr_e_fun = ca.Function('conl_constr_e_fun', [model.x, model.p, model.p_global], [conl_constr_expr_e])
            g.append(conl_constr_e_fun(xtraj[:, -1], ptraj[:, -1], model.p_global))
            lbg.append(constraints.lphi_e)
            ubg.append(constraints.uphi_e)

        ### Cost
        # initial cost term
        nlp_cost = 0
        cost_expr_0 = ocp.get_initial_cost_expression()
        cost_fun_0 = ca.Function('cost_fun_0', [model.x, model.u, model.p, model.p_global], [cost_expr_0])
        nlp_cost += solver_options.cost_scaling[0] * cost_fun_0(xtraj[:, 0], utraj[:, 0], ptraj[:, 0], model.p_global)

        # intermediate cost term
        cost_expr = ocp.get_path_cost_expression()
        cost_fun = ca.Function('cost_fun', [model.x, model.u, model.p, model.p_global], [cost_expr])
        for i in range(1, solver_options.N_horizon):
            nlp_cost += solver_options.cost_scaling[i] * cost_fun(xtraj[:, i], utraj[:, i], ptraj[:, i], model.p_global)

        # terminal cost term
        cost_expr_e = ocp.get_terminal_cost_expression()
        cost_fun_e = ca.Function('cost_fun_e', [model.x, model.p, model.p_global], [cost_expr_e])
        nlp_cost += solver_options.cost_scaling[-1] * cost_fun_e(xtraj[:, -1], ptraj[:, -1], model.p_global)

        # call w all primal variables
        w = ca.vertcat(ca.vec(xtraj), ca.vec(utraj))
        lbw = ca.vertcat(ca.vec(lb_xtraj), ca.vec(lb_utraj))
        ubw = ca.vertcat(ca.vec(ub_xtraj), ca.vec(ub_utraj))
        p_nlp = ca.vertcat(ca.vec(ptraj), model.p_global)

        # create NLP
        nlp = {"x": w, "p": p_nlp, "g": ca.vertcat(*g), "f": nlp_cost}
        bounds = {"lbx": lbw, "ubx": ubw, "lbg": ca.vertcat(*lbg), "ubg": ca.vertcat(*ubg)}

        return nlp, bounds


    def __init__(self, acados_ocp: AcadosOcp, verbose=True):

        if not isinstance(acados_ocp, AcadosOcp):
            raise TypeError('acados_ocp should be of type AcadosOcp.')

        self.acados_ocp = acados_ocp
        self.casadi_nlp, self.bounds = self.create_casadi_nlp_formulation(acados_ocp)
        self.casadi_solver = ca.nlpsol("nlp_solver", 'ipopt', self.casadi_nlp)
        self.nlp_sol = None


    def solve_for_x0(self, x0_bar, fail_on_nonzero_status=True, print_stats_on_failure=True):
         raise NotImplementedError()


    def solve(self) -> int:
        """
        Solve the ocp with current input.

        :return: status of the solver
        """
        self.nlp_sol = self.casadi_solver(lbx=self.bounds['lbx'], ubx=self.bounds['ubx'], lbg=self.bounds['lbg'], ubg=self.bounds['ubg'])

        # TODO: return correct status
        return 0

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
        :param field: string in ['x', 'u', 'z', 'pi', 'lam', 'sl', 'su', 'p', 'sens_u', 'sens_pi', 'sens_x', 'sens_lam', 'sens_sl', 'sens_su']

        .. note:: regarding lam: \n
                the inequalities are internally organized in the following order: \n
                [ lbu lbx lg lh lphi ubu ubx ug uh uphi; \n
                lsbu lsbx lsg lsh lsphi usbu usbx usg ush usphi]

        .. note:: pi: multipliers for dynamics equality constraints \n
                      lam: multipliers for inequalities \n
                      t: slack variables corresponding to evaluation of all inequalities (at the solution) \n
                      sl: slack variables of soft lower inequality constraints \n
                      su: slack variables of soft upper inequality constraints \n
        """
        if not isinstance(stage, int):
            raise TypeError('stage should be integer.')
        if self.nlp_sol is None:
            raise ValueError('No solution available. Please call solve() first.')
        dims = self.acados_ocp.dims

        if field == 'x':
            sol_w = self.nlp_sol['x']
            return sol_w[stage*dims.nx:(stage+1)*dims.nx].full().flatten()
        elif field == 'u':
            sol_w = self.nlp_sol['x']
            offset_x = dims.nx*(self.acados_ocp.solver_options.N_horizon+1)
            return sol_w[offset_x+stage*dims.nu: offset_x+(stage+1)*dims.nu].full().flatten()
        else:
            raise NotImplementedError(f"Field '{field}' is not implemented in AcadosCasadiOcpSolver")

    def get_flat(self, field_: str) -> np.ndarray:
        """
        Get concatenation of all stages of last solution of the solver.

        :param field: string in ['x', 'u', 'z', 'pi', 'lam', 'sl', 'su', 'p', 'p_global']

        .. note:: The parameter 'p_global' has no stage-wise structure and is processed in a memory saving manner by default. \n
                In order to read the 'p_global' parameter, the option 'save_p_global' must be set to 'True' upon instantiation. \n
        """
        raise NotImplementedError()


    def set_flat(self, field_: str, value_: np.ndarray) -> None:
        """
        Set concatenation solver initialization.

        :param field: string in ['x', 'u', 'z', 'pi', 'lam', 'sl', 'su', 'p']
        """
        raise NotImplementedError()


    def load_iterate(self, filename:str, verbose: bool = True):
        raise NotImplementedError()

    def store_iterate_to_obj(self) -> AcadosOcpIterate:
        raise NotImplementedError()

    def load_iterate_from_obj(self, iterate: AcadosOcpIterate):
        raise NotImplementedError()

    def store_iterate_to_flat_obj(self) -> AcadosOcpFlattenedIterate:
        raise NotImplementedError()

    def load_iterate_from_flat_obj(self, iterate: AcadosOcpFlattenedIterate) -> None:
        raise NotImplementedError()

    def get_stats(self, field_: str) -> Union[int, float, np.ndarray]:
        raise NotImplementedError()

    def get_cost(self) -> float:
        raise NotImplementedError()

    def set(self, stage_: int, field_: str, value_: np.ndarray):
        raise NotImplementedError()

    def cost_get(self, stage_: int, field_: str) -> np.ndarray:
        raise NotImplementedError()

    def cost_set(self, stage_: int, field_: str, value_):
        raise NotImplementedError()