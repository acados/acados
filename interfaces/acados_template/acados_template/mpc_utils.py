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
from array import array
from copy import deepcopy
from typing import Tuple, Optional
import casadi as ca
import numpy as np

from . import AcadosOcpIterate, AcadosOcpSolver
from .acados_model import AcadosModel
from .acados_ocp import AcadosOcp, AcadosOcpConstraints
from .utils import casadi_length, is_empty


class AcadosCostConstraintEvaluator:
    """
    This class provides convenience methods to evaluate the cost and constraints related to an AcadosOcp definition.
    A typical use case would be a closed-loop evaluation with the same standard costs and slack costs as defined
    in the original AcadosOcp.

    The evaluator must be updated by the method 'update_all(acados_solver)', if parameters in the solvers are changed.

    Two methods can be used for evaluation:
    - evaluate(x, u, step): evaluates the cost and constraints at a given stage of the OCP. For a closed-loop evaluation
        the stage is typically 0.
    - evaluate_ocp_cost(acados_ocp_iterate): evaluates the cost of a whole OCP trajectory, as evaluated inside acados.

    Limitation: values of numerical properties, such as bound values, W, zl, zu, Zu, Zl, yref,
    etc. are taken from original AcadosOcp. Their changes during runtime are not taken into account in the evaluator.
    """

    def __init__(self, ocp: AcadosOcp, with_parametric_bounds: bool = False):
        ocp.make_consistent()

        if with_parametric_bounds:
            raise NotImplementedError(
                "AcadosCostConstraintEvaluatator: with_parametric_bounds not implemented.")
            # p_bounds = ca.MX.sym('p_bounds', 2*(ocp.dims.nh+ocp.dims.nbx+dims.nbu+dims.ng+dims.nphi), 1)

        self.__ocp = deepcopy(ocp)

        model = ocp.model
        constraints = ocp.constraints
        dims = ocp.dims
        if casadi_length(ocp.model.z) > 0:
            raise NotImplementedError(
                "AcadosCostConstraintEvaluatator: not implemented for models with z.")

        self.__parameter_values = np.tile(ocp.parameter_values, (ocp.dims.np, ocp.dims.N))
        self.__p_global_values = ocp.p_global_values

        self.time_steps = ocp.solver_options.time_steps

        # setup casadi functions for constraints and cost
        cost_expr = get_path_cost_expression(ocp)
        cost_expr_e = get_terminal_cost_expression(ocp)

        p_global = ocp.model.p_global
        cost_fun_args = [ocp.model.x, ocp.model.u, ocp.model.p, p_global]

        self.cost_fun = ca.Function(
            'cost_fun',
            cost_fun_args,
            [cost_expr]
        )
        cost_fun_args_e = [ocp.model.x, ocp.model.p, p_global]

        self.terminal_cost_fun = ca.Function(
            'cost_fun_e',
            cost_fun_args_e,
            [cost_expr_e]
        )

        # constraints
        bu_expr = model.u[constraints.idxbu]
        bx_expr = model.x[constraints.idxbx]
        bx_expr_e = model.x[constraints.idxbx_e]

        if not is_empty(constraints.C):
            g_expr = constraints.C @ model.x + constraints.D @ model.u
        else:
            g_expr = ca.DM.zeros(0, 1)

        if not is_empty(constraints.C_e):
            g_expr_e = constraints.C_e @ model.x
        else:
            g_expr_e = ca.DM.zeros(0, 1)

        h_expr = model.con_h_expr
        h_expr_e = model.con_h_expr_e

        if ocp.dims.nphi > 0:
            raise NotImplementedError("AcadosCostConstraintEvaluatator: not implemented for nontrivial phi.")

        constraint_args = [bu_expr, bx_expr, g_expr, h_expr]
        constraint_args_e = [bx_expr_e, g_expr_e, h_expr_e]

        constraint_expr = ca.vertcat(*constraint_args)
        upper_bound = ca.vertcat(
            constraints.ubu,
            constraints.ubx,
            constraints.ug,
            constraints.uh,
            constraints.uphi
        )
        lower_bound = ca.vertcat(
            constraints.lbu,
            constraints.lbx,
            constraints.lg,
            constraints.lh,
            constraints.lphi
        )

        constraint_expr_e = ca.vertcat(*constraint_args_e)
        upper_bound_e = ca.vertcat(
            constraints.ubx_e,
            constraints.ug_e,
            constraints.uh_e,
            constraints.uphi_e
        )
        lower_bound_e = ca.vertcat(
            constraints.lbx_e,
            constraints.lg_e,
            constraints.lh_e,
            constraints.lphi_e
        )

        lower_violation = ca.fmax(lower_bound - constraint_expr, 0)
        upper_violation = ca.fmax(constraint_expr - upper_bound, 0)

        lower_violation_e = ca.fmax(lower_bound_e - constraint_expr_e, 0)
        upper_violation_e = ca.fmax(constraint_expr_e - upper_bound_e, 0)

        slack_indices = np.concatenate((
            constraints.idxsbu,
            constraints.idxsbx + dims.nbu,
            constraints.idxsg + dims.nbu + dims.nbx,
            constraints.idxsh + dims.nbu + dims.nbx + dims.ng,
            constraints.idxsphi + dims.nbu + dims.nbx + dims.ng + dims.nh
        ))

        slack_indices_e = np.concatenate((
            constraints.idxsbx_e,
            constraints.idxsg_e + dims.nbx_e,
            constraints.idxsh_e + dims.nbx_e + dims.ng_e,
            constraints.idxsphi_e + dims.nbx_e + dims.ng_e + dims.nh_e
        ))

        self.nonslacked_indices = np.setdiff1d(
            np.arange(constraint_expr.shape[0]),
            slack_indices
        )

        self.nonslacked_indices_e = np.setdiff1d(
            np.arange(constraint_expr_e.shape[0]),
            slack_indices_e
        )

        lower_slack_expression = lower_violation[slack_indices]
        upper_slack_expression = upper_violation[slack_indices]

        lower_slack_expression_e = lower_violation_e[slack_indices_e]
        upper_slack_expression_e = upper_violation_e[slack_indices_e]

        self.constraint_function = ca.Function(
            'constraint_function',
            cost_fun_args,
            [lower_violation, upper_violation, lower_slack_expression, upper_slack_expression]
        )

        self.constraint_function_e = ca.Function(
            'constraint_function_e',
            cost_fun_args_e,
            [lower_violation_e, upper_violation_e, lower_slack_expression_e, upper_slack_expression_e]
        )

    @property
    def parameter_values(self):
        """:math:`p` - initial values for parameter vector - can be updated stage-wise"""
        return self.__parameter_values

    @property
    def p_global_values(self):
        r"""initial values for :math:`p_\text{global}` vector,
        Type: `numpy.ndarray` of shape `(np_global, )`.
        """
        return self.__p_global_values

    def update_all(self, acados_solver: AcadosOcpSolver):
        """
        Update the parameter values and global parameter values from the acados solver.
        ATTENTION: Currently only parameter values are updated. Reference values, bounds, etc. are not updated.
        """
        N = self.__ocp.dims.N
        if self.__ocp.dims.np > 0:
            for i in range(N):
                self.__parameter_values[:, i] = acados_solver.get(i, 'p')

        if acados_solver.save_p_global is False:
            print('\nCan not set \'p_global\', since the solver does not store these values by default. '
                  'In order to update \'p_global\', please set the option \'save_p_global=True\' in the '
                  'instantiation of the acados solver.')
        else:
            self.__p_global_values = acados_solver.get_flat('p_global')

    def _get_check_parameters(
            self, p: np.ndarray = None,
            p_global: np.ndarray = None) -> Tuple[np.ndarray, np.ndarray]:
        parameter_values_return = self.__parameter_values
        p_global_values_return = self.__p_global_values

        if p is not None:
            if p.shape == self.__parameter_values.shape:
                parameter_values_return = p
            else:
                raise ValueError(
                    f"Parameter vector 'p' has wrong shape {p.shape} instead of {self.__parameter_values.shape}.")

        if p_global is not None:
            if p_global.shape == self.__p_global_values.shape:
                p_global_values_return = p_global
            else:
                raise ValueError(
                    f"Global parameter vector 'p_global' has wrong shape {p_global.shape} instead of {self.__p_global_values.shape}.")
        return parameter_values_return, p_global_values_return

    def evaluate(self, x: np.ndarray,
                 u: np.ndarray,
                 step: int = 0,
                 p: np.ndarray = None,
                 p_global: np.ndarray = None,
                 ) -> dict:
        """
        Evaluates the cost and constraints at a given stage of the OCP. For a closed-loop evaluation the stage is
        typically 0. If parameter values are provided, they are also set in the evaluator.
        @param x: state vector
        @param u: control input vector
        @param step: stage index (0 <= stage < N)
        @param p: parameter vector used for evaluation (optional)
        @param p_global: global parameter vector used for evaluation (optional)

        @return: dictionary with the following keys:
            - 'cost': total cost of the transition
            - 'cost_without_slacks': total cost of transition without slack penalties
            - 'cost_slacks': total cost of slack penalties
            - 'violation_soft_constraints': individual violation of soft constraints (equal to slacks)
            - 'violation_hard_constraints': individual violation of hard constraints
        """
        parameter_values, p_global_values = self._get_check_parameters(p, p_global)

        if len(parameter_values) > 0:
            parameter_values = parameter_values[:, step]
        else:
            parameter_values = parameter_values
        cost_fun_args = [x, u, parameter_values, p_global_values]

        # evaluate cost
        cost_without_slacks = self.cost_fun(*cost_fun_args).full() * self.time_steps[step]

        # evaluate constraints
        lower_violation, upper_violation, lower_slack, upper_slack = (
            self.constraint_function(x, u,
                                     parameter_values,
                                     p_global_values))
        violation_hard_constraints = np.concatenate(
            (lower_violation[self.nonslacked_indices], upper_violation[self.nonslacked_indices]))

        # evaluate cost of soft constraints
        lower_slack_cost, upper_slack_cost = 0.,0.

        if self.__ocp.cost.Zl.size > 0:
            lower_slack_cost += 0.5 * self.__ocp.cost.Zl @ (lower_slack.full() * lower_slack.full())
        if self.__ocp.cost.zl.size > 0:
            lower_slack_cost += self.__ocp.cost.zl @ lower_slack.full()

        if self.__ocp.cost.Zu.size > 0:
            upper_slack_cost += 0.5 * self.__ocp.cost.Zu @ (upper_slack.full()* upper_slack.full())
        if self.__ocp.cost.zu.size > 0:
            upper_slack_cost += self.__ocp.cost.zu @ upper_slack.full()

        slack_cost = (lower_slack_cost + upper_slack_cost) * self.time_steps[step]


        if len(slack_cost) == 0:
            cost = cost_without_slacks
        else:
            cost = cost_without_slacks + slack_cost

        # evaluate sum
        result = {
            'cost': cost.item(),
            'cost_without_slacks': cost_without_slacks.item(),
            'cost_slacks': slack_cost.item(),
            'violation_soft_constraints': np.concatenate((lower_slack.full(), upper_slack.full())),
            'violation_hard_constraints': violation_hard_constraints,
        }
        return result

    def evaluate_ocp_cost(
            self,
            acados_ocp_iterate: AcadosOcpIterate,
            p: np.ndarray = None,
            p_global: np.ndarray = None):
        """
        Evaluates the cost of a whole OCP trajectory, as evaluated inside acados.
        If parameter values are provided, they are also set in the evaluator.
        @param acados_ocp_iterate: acados OCP iterate object
        @param p: parameter vector used for evaluation (optional)
        @param p_global: global parameter vector used for evaluation (optional)
        """
        parameter_values, p_global_values = self._get_check_parameters(p, p_global)
        cost = 0

        # the cost on the first step is different in the OCP
        step = 0
        result = self.evaluate(acados_ocp_iterate.x_traj[0], acados_ocp_iterate.u_traj[0], step=step)
        cost += result['cost_without_slacks']
        step += 1

        for x_traj, u_traj in zip(acados_ocp_iterate.x_traj[1:], acados_ocp_iterate.u_traj[1:]):
            result = self.evaluate(x_traj, u_traj, step=step)
            cost += result['cost']
            step += 1

        if len(parameter_values) > 0:
            parameter_values = parameter_values[:, -1]
        else:
            parameter_values = parameter_values

        cost_fun_args = [acados_ocp_iterate.x_traj[-1], parameter_values, p_global_values]
        cost += self.terminal_cost_fun(*cost_fun_args).full()

        lower_violation_e, upper_violation_e, lower_slack_e, upper_slack_e = (
            self.constraint_function_e(acados_ocp_iterate.x_traj[-1],
                                       parameter_values,
                                       p_global_values))

        violation_hard_constraints = np.concatenate(
            (lower_violation_e[self.nonslacked_indices_e],
             upper_violation_e[self.nonslacked_indices_e]))

        # evaluate cost of soft constraints
        lower_slack_cost_e, upper_slack_cost_e = np.array([0.]), np.array([0.])
        if self.__ocp.cost.Zl_e.size > 0:
            lower_slack_cost_e += 0.5 * lower_slack_e.full().transpose()@self.__ocp.cost.Zl_e @ lower_slack_e.full()

        if self.__ocp.cost.zl_e.size > 0:
            lower_slack_cost_e += self.__ocp.cost.zl_e @ lower_slack_e.full()


        if self.__ocp.cost.Zu_e.size > 0:
            upper_slack_cost_e += 0.5 * upper_slack_e.full().transpose()@self.__ocp.cost.Zu_e @  upper_slack_e.full()

        if self.__ocp.cost.zu_e.size > 0:
            upper_slack_cost_e += self.__ocp.cost.zu_e @ upper_slack_e.full()

        if lower_slack_e.full().size > 0:
            cost += lower_slack_cost_e
        if upper_slack_e.full().size > 0:
            cost += upper_slack_cost_e

        return cost[0][0]


def get_path_cost_expression(ocp: AcadosOcp):
    # multiply path costs with td, nonuniform grid
    # multiply with cost scaling
    # cost scaling ueberschreibt tds, falls gesetzt
    model = ocp.model
    if ocp.cost.cost_type == "LINEAR_LS":
        y = ocp.cost.Vx @ model.x + ocp.cost.Vu @ model.u

        if casadi_length(model.z) > 0:
            y += ocp.cost.Vz @ model.z
        residual = y - ocp.cost.yref
        cost_dot = 0.5 * (residual.T @ ocp.cost.W @ residual)

    elif ocp.cost.cost_type == "NONLINEAR_LS":
        residual = model.cost_y_expr - ocp.cost.yref
        cost_dot = 0.5 * (residual.T @ ocp.cost.W @ residual)

    elif ocp.cost.cost_type == "EXTERNAL":
        cost_dot = model.cost_expr_ext_cost

    elif ocp.cost.cost_type == "CONVEX_OVER_NONLINEAR":
        raise NotImplementedError(
            "get_terminal_cost_expression: not implemented for CONVEX_OVER_NONLINEAR.")
        #cost_dot = ca.substitute(
        #    model.cost_psi_expr, model.cost_r_in_psi_expr, model.cost_y_expr)
    else:
        raise Exception("create_model_with_cost_state: Unknown cost type.")

    return cost_dot


def get_terminal_cost_expression(ocp: AcadosOcp):
    model = ocp.model
    if ocp.cost.cost_type_e == "LINEAR_LS":
        y = ocp.cost.Vx_e @ model.x
        residual = y - ocp.cost.yref_e
        cost_dot = 0.5 * (residual.T @ ocp.cost.W_e @ residual)

    elif ocp.cost.cost_type == "NONLINEAR_LS":
        residual = model.cost_y_expr_e - ocp.cost.yref_e
        cost_dot = 0.5 * (residual.T @ ocp.cost.W_e @ residual)

    elif ocp.cost.cost_type == "EXTERNAL":
        cost_dot = model.cost_expr_ext_cost_e

    elif ocp.cost.cost_type == "CONVEX_OVER_NONLINEAR":
        raise NotImplementedError(
            "get_terminal_cost_expression: not implemented for CONVEX_OVER_NONLINEAR.")
        #cost_dot = ca.substitute(
        #    model.cost_psi_expr_e, model.cost_r_in_psi_expr_e, model.cost_y_expr_e)
    else:
        raise Exception("create_model_with_cost_state: Unknown terminal cost type.")

    return cost_dot


def create_model_with_cost_state(ocp: AcadosOcp) -> Tuple[AcadosModel, np.ndarray]:
    """
    Creates a new AcadosModel with an extra state `cost_state`,
    which has the dynamics of the cost function and slack penalties corresponding to the intermediate shooting nodes.

    Note: In contrast cost_discretization='INTEGRATOR', this allows to integrate also the slack penalties.
    Since l1 slack penalties are nondifferentiable, an accurate cost integration with the model created by this function should use many integrator steps, when slack penalties are part of the OCP formulation.

    Returns the augmented model and the parameter values of the given AcadosOcp.
    """

    model = deepcopy(ocp.model)
    symbol = model.get_casadi_symbol()
    cost_state = symbol("cost_state")
    cost_state_dot = symbol("cost_state_dot")

    cost_dot = get_path_cost_expression(ocp)

    i_slack = 0
    for ibu in ocp.constraints.idxsbu:
        iu = ocp.constraints.idxbu[ibu]
        lower_violation = ca.fmax(ocp.constraints.lbu[ibu] - model.u[iu], 0)
        upper_violation = ca.fmax(model.u[iu] - ocp.constraints.ubu[ibu], 0)
        cost_dot += ocp.cost.zl[i_slack] * lower_violation + \
                    ocp.cost.Zl[i_slack] * lower_violation ** 2
        cost_dot += ocp.cost.zu[i_slack] * upper_violation + \
                    ocp.cost.Zu[i_slack] * upper_violation ** 2
        i_slack += 1

    for ibx in ocp.constraints.idxsbx:
        ix = ocp.constraints.idxbx[ibx]
        lower_violation = ca.fmax(ocp.constraints.lbx[ibx] - model.x[ix], 0)
        upper_violation = ca.fmax(model.x[ix] - ocp.constraints.ubx[ibx], 0)
        cost_dot += ocp.cost.zl[i_slack] * lower_violation + \
                    ocp.cost.Zl[i_slack] * lower_violation ** 2
        cost_dot += ocp.cost.zu[i_slack] * upper_violation + \
                    ocp.cost.Zu[i_slack] * upper_violation ** 2
        i_slack += 1

    if not is_empty(ocp.constraints.C):
        g = ocp.constraints.C @ ocp.model.x + ocp.constraints.D @ ocp.model.u
        for ig in ocp.constraints.idxsg:
            lower_violation = ca.fmax(ocp.constraints.lg[ig] - g[ig], 0)
            upper_violation = ca.fmax(g[ig] - ocp.constraints.ug[ig], 0)
            cost_dot += ocp.cost.zl[i_slack] * lower_violation + \
                        ocp.cost.Zl[i_slack] * lower_violation ** 2
            cost_dot += ocp.cost.zu[i_slack] * upper_violation + \
                        ocp.cost.Zu[i_slack] * upper_violation ** 2
            i_slack += 1

    for ih in ocp.constraints.idxsh:
        lower_violation = ca.fmax(
            ocp.constraints.lh[ih] - ocp.model.con_h_expr[ih], 0)
        upper_violation = ca.fmax(
            ocp.model.con_h_expr[ih] - ocp.constraints.uh[ih], 0)
        cost_dot += ocp.cost.zl[i_slack] * lower_violation + \
                    ocp.cost.Zl[i_slack] * lower_violation ** 2
        cost_dot += ocp.cost.zu[i_slack] * upper_violation + \
                    ocp.cost.Zu[i_slack] * upper_violation ** 2
        i_slack += 1

    if not is_empty(ocp.constraints.idxsphi):
        raise NotImplementedError(
            f"Not implemented for nontrivial ocp.constraints.idxsphi = {ocp.constraints.idxsphi}")

    model.x = ca.vertcat(model.x, cost_state)
    model.xdot = ca.vertcat(model.xdot, cost_state_dot)
    model.f_expl_expr = ca.vertcat(model.f_expl_expr, cost_dot)
    model.f_impl_expr = ca.vertcat(model.f_impl_expr, cost_state_dot - cost_dot)

    return model, ocp.parameter_values


def detect_constraint_structure(model: AcadosModel, constraints: AcadosOcpConstraints, stage_type: str):
    """
    - stage_type: allowed values "initial", "path", "terminal"
    """
    x = model.x
    u = model.u
    z = model.z
    nx = x.shape[0]
    nu = u.shape[0]

    p_global = model.p_global
    casadi_var = model.get_casadi_symbol()
    casadi_dm_zeros = ca.DM.zeros

    if stage_type == 'initial':
        expr_constr = model.con_h_expr_0
        lb = constraints.lh_0
        ub = constraints.uh_0
        print('\nConstraint detection for initial constraints.')
    elif stage_type == 'path':
        expr_constr = model.con_h_expr
        lb = constraints.lh
        ub = constraints.uh
        print('\nConstraint detection for path constraints.')
    elif stage_type == 'terminal':
        expr_constr = model.con_h_expr_e
        lb = constraints.lh_e
        ub = constraints.uh_e
        print('\nConstraint detection for terminal constraints.')
    else:
        raise ValueError('Constraint detection: Wrong stage_type.')

    if is_empty(expr_constr):
        expr_constr = casadi_var('con_h_expr', 0, 0)

    # Initialize
    constr_expr_h = casadi_var('con_h_expr', 0, 0)
    lh = []
    uh = []

    c_lin = casadi_dm_zeros(0, nx)
    d_lin = casadi_dm_zeros(0, nu)
    lg = []
    ug = []

    Jbx = casadi_dm_zeros(0, nx)
    lbx = []
    ubx = []

    Jbu = casadi_dm_zeros(0, nu)
    lbu = []
    ubu = []

    # Loop over CasADi formulated constraints
    for ii in range(expr_constr.shape[0]):
        c = expr_constr[ii]
        if any(ca.which_depends(c, z)) or \
                not ca.is_linear(c, ca.vertcat(x, u)) or \
                any(ca.which_depends(c, model.p)) or \
                any(ca.which_depends(c, p_global)):

            # External constraint
            constr_expr_h = ca.vertcat(constr_expr_h, c)
            lh.append(lb[ii])
            uh.append(ub[ii])
            print(f'constraint {ii + 1} is kept as nonlinear constraint.')
            print(c)
            print(' ')
        else:  # c is linear in x and u
            Jc_fun = ca.Function('Jc_fun', [x, u], [
                ca.jacobian(c, ca.vertcat(x, u))])
            Jc = Jc_fun(0, 0)
            if np.sum(Jc != 0) == 1:
                # c is bound
                idb = Jc.full().squeeze().nonzero()[0][0]
                if idb < nx:
                    # Bound on x
                    Jbx = ca.vertcat(Jbx, casadi_dm_zeros(1, nx))
                    Jbx[-1, idb] = 1
                    lbx.append(lb[ii] / Jc[idb])
                    ubx.append(ub[ii] / Jc[idb])
                    print(f'constraint {ii + 1} is reformulated as bound on x.')
                    print(c)
                    print(' ')
                else:
                    # Bound on u
                    Jbu = ca.vertcat(Jbu, casadi_dm_zeros(1, nu))
                    Jbu[-1, idb - nx] = 1
                    lbu.append(lb[ii] / Jc[idb])
                    ubu.append(ub[ii] / Jc[idb])
                    print(f'constraint {ii + 1} is reformulated as bound on u.')
                    print(c)
                    print(' ')
            else:
                # c is general linear constraint
                c_lin = ca.vertcat(c_lin, Jc[:nx])
                d_lin = ca.vertcat(d_lin, Jc[nx:])
                lg.append(lb[ii])
                ug.append(ub[ii])
                print(
                    f'constraint {ii + 1} is reformulated as general linear constraint.')
                print(c)
                print(' ')

    if stage_type == 'terminal':
        # Checks
        if any(ca.which_depends(expr_constr, u)) or lbu or (d_lin.size()[0] > 0 and any(d_lin)):
            raise ValueError(
                'Terminal constraint may not depend on control input.')
        # h
        constraints.constr_type_e = 'BGH'
        if lh:
            model.con_h_expr_e = constr_expr_h
            constraints.lh_e = np.array(lh)
            constraints.uh_e = np.array(uh)
        else:
            model.con_h_expr_e = None
            constraints.lh_e = np.array([])
            constraints.uh_e = np.array([])
        # linear constraint g
        if lg:
            constraints.C_e = c_lin
            constraints.lg_e = np.array(lg)
            constraints.ug_e = np.array(ug)
        # Bounds x
        if lbx:
            constraints.idxbx_e = J_to_idx(Jbx)
            constraints.lbx_e = np.array(lbx)
            constraints.ubx_e = np.array(ubx)

    elif stage_type == 'initial':
        print("At initial stage, only h constraints are detected.")
        constraints.constr_type_0 = 'BGH'
        # h
        if lh:
            model.con_h_expr_0 = constr_expr_h
            constraints.lh_0 = np.array(lh)
            constraints.uh_0 = np.array(uh)
        else:
            model.con_h_expr_0 = None
            constraints.lh_0 = np.array([])
            constraints.uh_0 = np.array([])
        # linear constraint g
        if lg:
            constraints.C = np.array(c_lin)
            constraints.D = np.array(d_lin)
            constraints.lg = np.array(lg)
            constraints.ug = np.array(ug)
        # Bounds x
        if lbx:
            constraints.idxbx_0 = J_to_idx(Jbx)
            constraints.lbx_0 = np.array(lbx)
            constraints.ubx_0 = np.array(ubx)
        # Bounds u
        if lbu:
            constraints.idxbu = J_to_idx(Jbu)
            constraints.lbu = np.array(lbu)
            constraints.ubu = np.array(ubu)
    else:  # path
        constraints.constr_type = 'BGH'
        # nonlinear constraint
        if lh:
            model.con_h_expr = constr_expr_h
            constraints.lh = np.array(lh)
            constraints.uh = np.array(uh)
        else:
            model.con_h_expr = None
            constraints.lh = np.array([])
            constraints.uh = np.array([])
        # linear constraint g
        if lg:
            constraints.C = np.array(c_lin)
            constraints.D = np.array(d_lin)
            constraints.lg = np.array(lg)
            constraints.ug = np.array(ug)
        # Bounds x
        if lbx:
            constraints.idxbx = J_to_idx(Jbx)
            constraints.lbx = np.array(lbx)
            constraints.ubx = np.array(ubx)
        # Bounds u
        if lbu:
            constraints.idxbu = J_to_idx(Jbu)
            constraints.lbu = np.array(lbu)
            constraints.ubu = np.array(ubu)


def J_to_idx(J):
    J = ca.sparsify(J)
    nrows = J.size()[0]
    idx = []
    for i in range(nrows):
        this_idx = ca.DM(J[i, :].sparsity()).full().flatten().nonzero()[0]
        if len(this_idx) != 1:
            raise ValueError(
                f'J_to_idx: Invalid J matrix. Exiting. Found more than one nonzero in row {i + 1}.')
        if J[i, this_idx] != 1:
            raise ValueError(
                f'J_to_idx: J matrices can only contain 1s, got J({i + 1}, {this_idx}) = {J[i, this_idx]}')
        idx.append(this_idx[0])
    return np.array(idx)
