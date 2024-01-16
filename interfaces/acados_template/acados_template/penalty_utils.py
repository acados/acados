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

from typing import Optional
import casadi as ca

def huber_loss(var: ca.SX, delta: float, tau: float):
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
            raise Exception("Either specify w or tau")
        tau = 2 * w * delta
    elif w is not None:
        raise Exception("Either specify w or tau")

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
            raise Exception("Either specify w or tau")
        tau = 2 * w * delta
    elif w is not None:
        raise Exception("Either specify w or tau")

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

