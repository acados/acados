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

from copy import deepcopy
from typing import Tuple
import casadi as ca
import numpy as np

from .acados_model import AcadosModel
from .acados_ocp import AcadosOcp, AcadosOcpConstraints
from .utils import casadi_length, is_empty


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

    if ocp.cost.cost_type == "LINEAR_LS":
        y = ocp.cost.Vx @ model.x + ocp.cost.Vu @ model.u
        if casadi_length(model.z) > 0:
            ocp.cost.Vz @ model.z
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
