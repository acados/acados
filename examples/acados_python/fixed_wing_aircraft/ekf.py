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


def build_ekf(model, integrator, W, V):
    """
    Build EKF functions.

    Args:
        model: AcadosModel from export_maverix_model()
        integrator: AcadosSimSolver, built with sens_forw and sens_forw_p
        W: process noise covariance of [x, wind], over one integrator step
        V: measurement noise covariance

    Returns:
        ekf: function(x, wind, P, y, u) -> (x, wind, P)
        measure: function(x_true, wind_true) -> y
    """

    # Build measurement g(x) and Jacobian C(x)
    x_sym = ca.vertcat(model.x, model.p)
    g = ca.DM([0, 0, 9.81])
    R = ca.reshape(x_sym[6:15], 3, 3)
    v_b = R.T @ (x_sym[3:6] - x_sym[22:25])
    ddq = model.f_expl_expr[3:6]

    y_sym = ca.vertcat(x_sym[0:3], x_sym[3:6], x_sym[15:18], ca.norm_2(v_b), R.T @ (ddq - g))
    g_func = ca.Function('g', [x_sym], [y_sym])
    C_func = ca.Function('C', [x_sym], [ca.jacobian(y_sym, x_sym)])

    def ekf(x, wind, P, y, u):
        x_aug = np.concatenate([x, wind])

        # Predict, over the estimator step and with the estimated wind
        x_aug[:22] = integrator.simulate(x=x, u=u, p=wind)
        x_pred = x_aug.copy()

        # state transition Jacobian from integrator
        A = np.eye(25)
        A[:22, :22] = integrator.get("Sx")
        A[:22, 22:] = integrator.get("S_p")
        P = A @ P @ A.T + W

        # Update
        C = np.array(C_func(x_aug))
        S = C @ P @ C.T + V
        K = P @ C.T @ np.linalg.pinv(S)
        x_aug = x_aug + K @ (y - np.array(g_func(x_aug)).flatten())

        # try re-orthonormalization
        try:
            x_aug[6:15] = re_ortho(x_aug[6:15])
        except:
            x_aug[6:15] = x_pred[6:15] # fallback to prediction
            print("Warning: SVD did not converge.")

        # Update covariance
        I_KC = np.eye(25) - K @ C
        P = I_KC @ P @ I_KC.T + K @ V @ K.T
        P = 0.5 * (P + P.T)

        return x_aug[:22].copy(), x_aug[22:25].copy(), P

    def measure(x_true, wind_true, noise=True):
        xf = np.concatenate([x_true, wind_true])
        y = np.array(g_func(xf)).flatten()
        if noise:
            y += np.random.randn(V.shape[0]) * np.sqrt(np.diag(V))
        return y

    return ekf, measure


def re_ortho(R_flat):
    """SVD re-orthonormalization of rotation matrix."""
    U, _, Vt = np.linalg.svd(R_flat.reshape(3, 3))
    return (U @ Vt).flatten()
