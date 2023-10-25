from typing import Optional

from .acados_ocp import AcadosOcp
from scipy.linalg import block_diag

import casadi as ca
import numpy as np


def formulate_constraint_as_L2_penalty(
    ocp: AcadosOcp,
    constr_expr: ca.SX,
    weight: float,
    upper_bound: Optional[float],
    lower_bound: Optional[float],
    residual_name: str = "new_residual",
) -> AcadosOcp:

    if upper_bound is None and lower_bound is None:
        raise ValueError("Either upper or lower bound must be provided.")

    # compute violation expression
    violation_expr = 0.0
    if upper_bound is not None:
        violation_expr = ca.fmax(violation_expr, (constr_expr - upper_bound))
    if lower_bound is not None:
        violation_expr = ca.fmax(violation_expr, (lower_bound - constr_expr))

    # add penalty as cost
    ocp.cost.yref = np.concatenate((ocp.cost.yref, np.zeros(1)))
    ocp.model.cost_y_expr = ca.vertcat(ocp.model.cost_y_expr, violation_expr)
    if ocp.cost.cost_type == "NONLINEAR_LS":
        ocp.cost.W = block_diag(ocp.cost.W, weight)
    elif ocp.cost.cost_type == "CONVEX_OVER_NONLINEAR":
        new_residual = ca.SX.sym(residual_name, constr_expr.shape)
        ocp.model.cost_r_in_psi_expr = ca.vertcat(ocp.model.cost_r_in_psi_expr, new_residual)
        ocp.model.cost_psi_expr += .5 * weight * new_residual**2

    return ocp


def formulate_constraint_as_Huber_penalty(
    ocp: AcadosOcp,
    constr_expr: ca.SX,
    weight: float,
    upper_bound: Optional[float],
    lower_bound: Optional[float],
    residual_name: str = "new_residual",
    huber_tau: float = 1.0,
    huber_delta: float = 1.0,
    use_xgn = True,
) -> AcadosOcp:

    if upper_bound is None and lower_bound is None:
        raise ValueError("Either upper or lower bound must be provided.")

    if ocp.cost.cost_type != "CONVEX_OVER_NONLINEAR":
        raise Exception("Huber penalty is only supported for CONVEX_OVER_NONLINEAR cost type.")

    if (upper_bound is None or lower_bound is None):
        raise NotImplementedError("only symmetric Huber for now")

    if use_xgn and ocp.model.cost_conl_custom_outer_hess is None:
        # switch to XGN Hessian start with exact Hessian of previously defined cost
        exact_cost_hess = ca.hessian(ocp.model.cost_psi_expr, ocp.model.cost_r_in_psi_expr)[0]
        ocp.model.cost_conl_custom_outer_hess = exact_cost_hess

    # define residual
    new_residual = ca.SX.sym(residual_name, constr_expr.shape)
    ocp.model.cost_r_in_psi_expr = ca.vertcat(ocp.model.cost_r_in_psi_expr, new_residual)

    # define penalty
    penalty, penalty_grad, penalty_hess, penalty_hess_xgn = \
            symmetric_huber_penalty(new_residual, delta=huber_delta, tau=huber_tau)

    # add penalty to cost
    ocp.model.cost_psi_expr += weight * penalty
    ocp.model.cost_y_expr = ca.vertcat(ocp.model.cost_y_expr, constr_expr)
    ocp.cost.yref = np.concatenate((ocp.cost.yref, np.zeros(1)))

    # add Hessian term
    if use_xgn:
        zero_offdiag = ca.SX.zeros(ocp.model.cost_conl_custom_outer_hess.shape[0], penalty_hess_xgn.shape[1])
        ocp.model.cost_conl_custom_outer_hess = ca.blockcat(ocp.model.cost_conl_custom_outer_hess, zero_offdiag, zero_offdiag.T, weight * penalty_hess_xgn)
    elif ocp.model.cost_conl_custom_outer_hess is not None:
        zero_offdiag = ca.SX.zeros(ocp.model.cost_conl_custom_outer_hess.shape[0], penalty_hess_xgn.shape[1])
        # add penalty Hessian to existing Hessian
        ocp.model.cost_conl_custom_outer_hess = ca.blockcat(ocp.model.cost_conl_custom_outer_hess, zero_offdiag, zero_offdiag.T, weight * penalty_hess)

    return ocp


def huber_loss(var: ca.SX, delta: float, tau: float):
    loss = tau/delta * ca.if_else(
        ca.fabs(var) < delta,
        0.5*var**2,
        delta*(ca.fabs(var) - 0.5*delta))

    loss_hess, loss_grad = ca.hessian(loss, var)
    loss_hess_XGN = ca.if_else(
        var == 0,
        loss_hess,
        ca.diag(loss_grad / var))

    return loss, loss_grad, loss_hess, loss_hess_XGN


def symmetric_huber_penalty(u: ca.SX, delta: float, tau: Optional[float] = None, w: Optional[float] = None):
    """
    Symmetric Huber penalty for a constraint -1 <= u <= 1.
    delta: the length of the quadratic behavior
    w: hessian in quadratic region
    """

    assert delta >= 0 and delta <= 1

    if tau is None:
        if w is None:
            raise Exception('Either specify w or tau')
        tau = 2*w*delta
    elif w is not None:
        raise Exception('Either specify w or tau')

    loss, _, _, loss_hess_XGN = huber_loss(u, delta, tau)

    # barrier huber
    # penalty = 0.5*(ca.substitute(loss, u, u - 1) + ca.substitute(loss, u, u + 1) - ca.substitute(loss, u, -1) - ca.substitute(loss, u, 1))

    # shifted by delta to get a penalty
    penalty = 0.5 * (ca.substitute(loss, u, u - (1+delta)) + ca.substitute(loss, u, u + (1+delta)))
    penalty += 0.5 * (-ca.substitute(loss, u, -(1+delta)) - ca.substitute(loss, u, 1-delta)) - tau*delta

    penalty_hess, penalty_grad = ca.hessian(penalty, u)

    penalty_hess_xgn = 0.5*ca.if_else(u < 0, ca.substitute(loss_hess_XGN, u, u+1+delta), ca.substitute(loss_hess_XGN, u, u-1-delta))

    return penalty, penalty_grad, penalty_hess, penalty_hess_xgn


