#
# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# Licensed under the 2-Clause BSD License.

from typing import Optional, Tuple, Union
import casadi as ca
import numpy as np

def huber_loss(var: ca.SX, delta: float, tau: float) -> Tuple[ca.SX, ca.SX, ca.SX, ca.SX]:
    """
    Huber loss function.
    Return the loss, gradient, hessian and the XGN Hessian
    """
    loss = (tau / delta) * ca.if_else(
        ca.fabs(var) < delta, 0.5 * var**2, delta * (ca.fabs(var) - 0.5 * delta)
    )

    loss_hess, loss_grad = ca.hessian(loss, var)
    loss_hess_XGN = ca.if_else(var == 0, loss_hess, ca.diag(loss_grad / var))

    return loss, loss_grad, loss_hess, loss_hess_XGN


def one_sided_huber_penalty(
    u: ca.SX,
    delta: float,
    tau: Optional[float] = None,
    w: Optional[float] = None,
    min_hess: float = 0.0,
):
    """
    One-sided Huber penalty for a constraint u <= 0.
    Note: either tau or w need to be specified.
    delta: the length of the quadratic behavior
    tau: gradient in linear region
    w: hessian in quadratic region
    min_hess: provide a minimum value for the hessian
    """

    if delta < 0:
        raise ValueError("delta must be positive")

    if tau is None:
        if w is None:
            raise ValueError("Either specify w or tau")
        tau = 2 * w * delta
    elif w is not None:
        raise ValueError("Either specify w or tau")

    loss, _, _, loss_hess_XGN = huber_loss(u, delta, tau)
    # shifted by delta to get a penalty
    penalty = 0.5 * (ca.substitute(loss, u, u - delta) + tau * u)

    penalty_0 = ca.substitute(penalty, u, 0)
    penalty = penalty - penalty_0

    penalty_hess, penalty_grad = ca.hessian(penalty, u)

    penalty_hess_xgn = 0.5 * ca.substitute(loss_hess_XGN, u, u - delta)

    if min_hess > 0.0:
        penalty_hess = ca.fmax(min_hess, penalty_hess)
        penalty_hess_xgn = ca.fmax(min_hess, penalty_hess_xgn)

    return penalty, penalty_grad, penalty_hess, penalty_hess_xgn


def symmetric_huber_penalty(
    u: ca.SX,
    delta: float,
    tau: Optional[float] = None,
    w: Optional[float] = None,
    min_hess: float = 0.0,
):
    """
    Symmetric Huber penalty for a constraint -1 <= u <= 1.
    Note: either tau or w need to be specified.
    delta: the length of the quadratic behavior
    tau: gradient in linear region
    w: hessian in quadratic region
    min_hess: provide a minimum value for the hessian
    """

    if delta < 0:
        raise ValueError("delta must be positive")

    if tau is None:
        if w is None:
            raise ValueError("Either specify w or tau")
        tau = 2 * w * delta
    elif w is not None:
        raise ValueError("Either specify w or tau")

    loss, _, _, loss_hess_XGN = huber_loss(u, delta, tau)

    # barrier huber
    # penalty = 0.5*(ca.substitute(loss, u, u - 1) + ca.substitute(loss, u, u + 1) - ca.substitute(loss, u, -1) - ca.substitute(loss, u, 1))

    # shifted by delta to get a penalty
    penalty = 0.5 * (
        ca.substitute(loss, u, u - (1 + delta))
        + ca.substitute(loss, u, u + (1 + delta))
    )
    penalty += 0.5 * (
        -ca.substitute(loss, u, -(1 + delta)) - ca.substitute(loss, u, 1 - delta)
    )

    penalty_0 = ca.substitute(penalty, u, 0)
    penalty = penalty - penalty_0

    penalty_hess, penalty_grad = ca.hessian(penalty, u)

    penalty_hess_xgn = 0.5 * ca.if_else(
        u < 0,
        ca.substitute(loss_hess_XGN, u, u + 1 + delta),
        ca.substitute(loss_hess_XGN, u, u - 1 - delta),
    )

    if min_hess > 0.0:
        penalty_hess = ca.fmax(min_hess, penalty_hess)
        penalty_hess_xgn = ca.fmax(min_hess, penalty_hess_xgn)

    return penalty, penalty_grad, penalty_hess, penalty_hess_xgn


def get_quadratic_penalty_expression(h_expr: Union[ca.SX, ca.MX], lh: np.ndarray, uh: np.ndarray, Z_l: np.ndarray, Z_u: np.ndarray):
    """
    Returns a CasADi expression corresponding to a quadratic penalty on the constraint violation with quadratic weight diag(Z_l) for lower bound violations and diag(Z_u) for upper bound violations..
    """
    lower_violation = ca.fmax(lh - h_expr, 0)
    upper_violation = ca.fmax(h_expr - uh, 0)

    return 0.5 * ca.sum1(Z_l * lower_violation**2 + Z_u * upper_violation**2)