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
from .acados_ocp import AcadosOcp
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
        cost_dot = ca.substitute(model.cost_psi_expr, model.cost_r_in_psi_expr, model.cost_y_expr)

    else:
        raise Exception("create_model_with_cost_state: Unknown cost type.")

    i_slack = 0
    for ibu in ocp.constraints.idxsbu:
        iu = ocp.constraints.idxbu[ibu]
        lower_violation = ca.fmax(ocp.constraints.lbu[ibu] - model.u[iu], 0)
        upper_violation = ca.fmax(model.u[iu] - ocp.constraints.ubu[ibu], 0)
        cost_dot += ocp.cost.zl[i_slack] * lower_violation + ocp.cost.Zl[i_slack] * lower_violation ** 2
        cost_dot += ocp.cost.zu[i_slack] * upper_violation + ocp.cost.Zu[i_slack] * upper_violation ** 2
        i_slack += 1

    for ibx in ocp.constraints.idxsbx:
        ix = ocp.constraints.idxbx[ibx]
        lower_violation = ca.fmax(ocp.constraints.lbx[ibx] - model.x[ix], 0)
        upper_violation = ca.fmax(model.x[ix] - ocp.constraints.ubx[ibx], 0)
        cost_dot += ocp.cost.zl[i_slack] * lower_violation + ocp.cost.Zl[i_slack] * lower_violation ** 2
        cost_dot += ocp.cost.zu[i_slack] * upper_violation + ocp.cost.Zu[i_slack] * upper_violation ** 2
        i_slack += 1


    if not is_empty(ocp.constraints.C):
        g = ocp.constraints.C @ ocp.model.x + ocp.constraints.D @ ocp.model.u
        for ig in ocp.constraints.idxsg:
            lower_violation = ca.fmax(ocp.constraints.lg[ig] - g[ig], 0)
            upper_violation = ca.fmax(g[ig] - ocp.constraints.ug[ig], 0)
            cost_dot += ocp.cost.zl[i_slack] * lower_violation + ocp.cost.Zl[i_slack] * lower_violation ** 2
            cost_dot += ocp.cost.zu[i_slack] * upper_violation + ocp.cost.Zu[i_slack] * upper_violation ** 2
            i_slack += 1

    for ih in ocp.constraints.idxsh:
        lower_violation = ca.fmax(ocp.constraints.lh[ih] - ocp.model.con_h_expr[ih], 0)
        upper_violation = ca.fmax(ocp.model.con_h_expr[ih] - ocp.constraints.uh[ih], 0)
        cost_dot += ocp.cost.zl[i_slack] * lower_violation + ocp.cost.Zl[i_slack] * lower_violation ** 2
        cost_dot += ocp.cost.zu[i_slack] * upper_violation + ocp.cost.Zu[i_slack] * upper_violation ** 2
        i_slack += 1

    if not is_empty(ocp.constraints.idxsphi):
        raise NotImplementedError(f"Not implemented for nontrivial ocp.constraints.idxsphi = {ocp.constraints.idxsphi}")

    model.x = ca.vertcat(model.x, cost_state)
    model.xdot = ca.vertcat(model.xdot, cost_state_dot)
    model.f_expl_expr = ca.vertcat(model.f_expl_expr, cost_dot)
    model.f_impl_expr = ca.vertcat(model.f_impl_expr, cost_state_dot-cost_dot)

    return model, ocp.parameter_values
