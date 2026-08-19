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

from acados_template import AcadosModel
from casadi import DM, SX, atan2, cross, diag, inv, reshape, skew, vertcat

from utils import buildForcesanMomentsFunction


def export_maverix_model(aero_coeff_file: str = 'data/aerodynamic_coefficients.json') -> AcadosModel:
    model_name = 'aircraft'

    # constants
    m = 2.4
    g = DM([0, 0, 9.81])
    I = diag(DM([0.0920, 0.0980, 0.1720]))
    Iinv = inv(I)
    aero_func = buildForcesanMomentsFunction(aero_coeff_file)

    # states
    q      = SX.sym('q', 3, 1)
    dq     = SX.sym('dq', 3, 1)
    R_flat = SX.sym('R', 9, 1)  # column-major flattened
    omega  = SX.sym('omega', 3, 1)
    thrust = SX.sym('thrust', 2, 1)
    flaps  = SX.sym('flaps', 2, 1)
    x = vertcat(q, dq, R_flat, omega, thrust, flaps)

    # controls
    dthrust = SX.sym('dthrust', 2, 1)
    dflaps  = SX.sym('dflaps', 2, 1)
    u = vertcat(dthrust, dflaps)

    # parameter
    wind = SX.sym('wind', 3, 1)  # wind in inertial frame
    p = wind

    # dynamics
    R = reshape(R_flat, 3, 3)  # rotation matrix from body to inertial frame
    v_b = R.T @ (dq - wind)    # airspeed vector, velocity relative to air in body frame
    F_b, M_b = aero_func(v_b, omega, thrust, flaps)

    F = (R @ F_b)
    ddq = F / m + g
    domega = Iinv @ (M_b - cross(omega, I @ omega))

    # baumgarte stabilization
    roh = 1
    A = roh / 4 * (DM.eye(3) - R.T @ R)
    dR = R @ (skew(omega) + A)
    dR_flat = reshape(dR, 9, 1)

    # expressions
    xdot = SX.sym('xdot', x.numel(), 1)
    f_expl = vertcat(dq, ddq, dR_flat, domega, dthrust, dflaps)
    f_impl = xdot - f_expl

    # nonlinear constraint expressions
    angle_of_attack = atan2(v_b[2], v_b[0])
    sideslip_angle  = atan2(v_b[1], v_b[0])  # approx. for small aoa
    con_h_expr = vertcat(v_b, angle_of_attack, sideslip_angle)

    # build model
    model = AcadosModel()
    model.x           = x
    model.xdot        = xdot
    model.u           = u
    model.p           = p
    model.f_expl_expr = f_expl
    model.f_impl_expr = f_impl
    model.con_h_expr  = con_h_expr
    model.name        = model_name

    return model
