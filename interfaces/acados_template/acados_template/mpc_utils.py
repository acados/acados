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

from .acados_model import AcadosModel
from .acados_ocp import AcadosOcp, AcadosOcpConstraints
from .utils import casadi_length, is_empty, array_to_float


class AcadosCostConstraintEvaluator:
    """
    Limitation: values of numerical properties, such as bound values, W, zl, zu, Zu, Zl, yref, etc. are taken from original AcadosOcp;
    If they are changed in the solve this is not taken into account here.
    Could be generalized later, by making the casadi functions parametric in bounds.
    """

    def __init__(self, ocp: AcadosOcp, with_parametric_bounds: bool):
        ocp.make_consistent()

        if with_parametric_bounds:
            raise NotImplementedError(
                "AcadosCostConstraintEvaluatator: with_parametric_bounds not implemented.")
            # p_bounds = ca.MX.sym('p_bounds', 2*(ocp.dims.nh+ocp.dims.nbx+dims.nbu+dims.ng+dims.nphi), 1)

        self.ocp = ocp

        model = ocp.model
        constraints = ocp.constraints
        dims = ocp.dims
        if casadi_length(ocp.model.z) > 0:
            raise NotImplementedError(
                "AcadosCostConstraintEvaluatator: not implemented for models with z.")

        self.__parameter_values = ocp.parameter_values
        self.__p_global_values = ocp.p_global_values

        self.dt = ocp.solver_options.time_steps[0]

        # setup casadi functions for constraints and cost
        cost_expr = get_path_cost_expression(ocp)

        # build function. fields may be None or empty casadi expressions
        if ocp.model.p_global is None:
            p_global = ca.SX(0,0)
        else:
            p_global = ocp.model.p_global

        cost_fun_args = [ocp.model.x, ocp.model.u, ocp.model.p, p_global]
        cost_fun_args = [arg for arg in cost_fun_args if arg is not None]
        self.cost_fun = ca.Function(
            'cost_fun',
            cost_fun_args,
            [cost_expr]
        )

        # constraints
        bu_expr = model.u[constraints.idxbu]
        bx_expr = model.x[constraints.idxbx]
        if not is_empty(constraints.C):
            g_expr = constraints.C @ model.x + constraints.D @ model.u
        else:
            g_expr = ca.DM.zeros(0, 1)

        h_expr = model.con_h_expr

        if ocp.dims.nphi > 0:
            raise NotImplementedError("AcadosCostConstraintEvaluatator: not implemented for nontrivial phi.")

        # build function. fields may be None, which creates problems for casadi function
        constraint_args = [bu_expr, bx_expr, g_expr, h_expr]
        constraint_args = [arg for arg in constraint_args if arg is not None]

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

        lower_violation = ca.fmax(lower_bound - constraint_expr, 0)
        upper_violation = ca.fmax(constraint_expr - upper_bound, 0)

        slack_indices = np.concatenate((
            constraints.idxsbu,
            constraints.idxsbx + dims.nbu,
            constraints.idxsg + dims.nbu + dims.nbx,
            constraints.idxsh + dims.nbu + dims.nbx + dims.ng,
            constraints.idxsphi + dims.nbu + dims.nbx + dims.ng + dims.nh
        ))

        self.nonslacked_indices = np.setdiff1d(
            np.arange(constraint_expr.shape[0]),
            slack_indices
        )

        lower_slack_expression = lower_violation[slack_indices]
        upper_slack_expression = upper_violation[slack_indices]

        self.constraint_function = ca.Function(
            'constraint_function',
            cost_fun_args,
            [lower_violation, upper_violation, lower_slack_expression, upper_slack_expression]
        )

    @property
    def parameter_values(self):
        """:math:`p` - initial values for parameter vector - can be updated stagewise"""
        return self.__parameter_values

    @parameter_values.setter
    def parameter_values(self, parameter_values):
        if isinstance(parameter_values, np.ndarray):
            self.__parameter_values = parameter_values
        else:
            raise Exception('Invalid parameter_values value. ' +
                            f'Expected numpy array, got {type(parameter_values)}.')
        if parameter_values.shape[0] != self.ocp.dims.np:
            raise Exception('inconsistent dimension np, regarding model.p and parameter_values.' + \
                            f'\nGot np = {self.ocp.dims.np}, ' + \
                            f'self.parameter_values.shape = {parameter_values.shape[0]}\n')

    @property
    def p_global_values(self):
        r"""initial values for :math:`p_\text{global}` vector,
        Type: `numpy.ndarray` of shape `(np_global, )`.
        """
        return self.__p_global_values

    @p_global_values.setter
    def p_global_values(self, p_global_values):
        if isinstance(p_global_values, np.ndarray):
            self.__p_global_values = p_global_values
        else:
            raise Exception('Invalid p_global_values value. ' +
                            f'Expected numpy array, got {type(p_global_values)}.')
        if p_global_values.shape[0] != self.ocp.dims.np_global:
            raise Exception('inconsistent dimension np_global, regarding model.p_global and p_global_values.' + \
                f'\nGot np_global = {self.ocp.dims.np_global}, '
                f'self.p_global_values.shape = {p_global_values.shape[0]}\n')

    def evaluate(self,
                 x: np.ndarray, u: np.ndarray,
                 p: Optional[np.ndarray] = None,
                 p_global: Optional[np.ndarray] = None,
                 dt: Optional[float] = None) -> dict:
        if p is not None:
            self.parameter_values(p)
        if p_global is not None:
            self.p_global_values(p_global)
        if dt is None:
            dt = self.dt

        cost_fun_args = [x, u, self.__parameter_values, self.__p_global_values]
        # cost_fun_args = [arg for arg in cost_fun_args if arg is not None]

        # evaluate cost
        cost_without_slacks = self.cost_fun(*cost_fun_args).full()

        # evaluate constraints
        lower_violation, upper_violation, lower_slack, upper_slack = (
            self.constraint_function(x, u,
                                     self.__parameter_values,
                                     self.__p_global_values))
        violation_hard_constraints = np.concatenate(
            (lower_violation[self.nonslacked_indices], upper_violation[self.nonslacked_indices]))

        # evaluate cost of soft constraints
        lower_slack_cost = 0.5 * self.ocp.cost.Zl @ (lower_slack * lower_slack) + self.ocp.cost.zl @ lower_slack
        upper_slack_cost = 0.5 * self.ocp.cost.Zu @ (upper_slack * upper_slack) + self.ocp.cost.zu @ upper_slack
        slack_cost = (lower_slack_cost + upper_slack_cost).full() * dt
        if len(slack_cost)==0:
            cost = cost_without_slacks
        else:
            cost = cost_without_slacks + slack_cost


        # evaluate sum
        result = {
            'cost': array_to_float(cost),
            'cost_without_slacks': array_to_float(cost_without_slacks),
            'slack_cost': array_to_float(slack_cost),
            'violation_hard_constraints': violation_hard_constraints,
        }
        return result


def get_path_cost_expression(ocp: AcadosOcp):
    model = ocp.model
    if ocp.cost.cost_type == "LINEAR_LS":
        y = ocp.cost.Vx @ model.x + ocp.cost.Vu @ model.u

        if casadi_length(model.z) > 0:
            y += ocp.cost.Vz @ model.z
        residual = y - ocp.cost.yref
        cost_dot = 0.5*(residual.T @ ocp.cost.W @ residual)

    elif ocp.cost.cost_type == "NONLINEAR_LS":
        residual = model.cost_y_expr - ocp.cost.yref
        cost_dot = 0.5*(residual.T @ ocp.cost.W @ residual)

    elif ocp.cost.cost_type == "EXTERNAL":
        cost_dot = model.cost_expr_ext_cost

    elif ocp.cost.cost_type == "CONVEX_OVER_NONLINEAR":
        cost_dot = ca.substitute(
            model.cost_psi_expr, model.cost_r_in_psi_expr, model.cost_y_expr)
    else:
        raise Exception("create_model_with_cost_state: Unknown cost type.")

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
    model.f_impl_expr = ca.vertcat(model.f_impl_expr, cost_state_dot-cost_dot)

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

    if model.p_global is None:
        p_global = []  # to have same structure of model.p

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

    if expr_constr is None:
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
            print(f'constraint {ii+1} is kept as nonlinear constraint.')
            print(c)
            print(' ')
        else:  # c is linear in x and u
            Jc_fun = ca.Function('Jc_fun', [x[0]], [
                                 ca.jacobian(c, ca.vertcat(x, u))])
            Jc = Jc_fun(0)

            if np.sum(Jc != 0) == 1:
                # c is bound
                idb = Jc.full().squeeze().nonzero()[0][0]
                if idb < nx:
                    # Bound on x
                    Jbx = ca.vertcat(Jbx, casadi_dm_zeros(1, nx))
                    Jbx[-1, idb] = 1
                    lbx.append(lb[ii] / Jc[idb])
                    ubx.append(ub[ii] / Jc[idb])
                    print(f'constraint {ii+1} is reformulated as bound on x.')
                    print(c)
                    print(' ')
                else:
                    # Bound on u
                    Jbu = ca.vertcat(Jbu, casadi_dm_zeros(1, nu))
                    Jbu[-1, idb - nx] = 1
                    lbu.append(lb[ii] / Jc[idb])
                    ubu.append(ub[ii] / Jc[idb])
                    print(f'constraint {ii+1} is reformulated as bound on u.')
                    print(c)
                    print(' ')
            else:
                # c is general linear constraint
                c_lin = ca.vertcat(c_lin, Jc[:nx])
                d_lin = ca.vertcat(d_lin, Jc[nx:])
                lg.append(lb[ii])
                ug.append(ub[ii])
                print(
                    f'constraint {ii+1} is reformulated as general linear constraint.')
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
    nrows = J.size()[0]
    idx = []
    for i in range(nrows):
        this_idx = ca.DM(J[i, :].sparsity()).full().nonzero()[0]
        if len(this_idx) != 1:
            raise ValueError(
                f'J_to_idx: Invalid J matrix. Exiting. Found more than one nonzero in row {i+1}.')
        if J[i, this_idx] != 1:
            raise ValueError(
                f'J_to_idx: J matrices can only contain 1s, got J({i+1}, {this_idx}) = {J[i, this_idx]}')
        idx.append(this_idx[0])
    return np.array(idx)
