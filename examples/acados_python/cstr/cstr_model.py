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

# authors: Katrin Baumgaertner, Jonathan Frey

import numpy as np
from dataclasses import dataclass
from acados_template import AcadosModel
from casadi import SX, vertcat, exp
from cstr_utils import compute_lqr_gain


@dataclass
class CSTRParameters:
    # nominal parameter values
    F0: float = 0.1  # m^3/min
    T0: float = 350.0  # K
    c0: float = 1.0  # kmol/m^3
    r: float = 0.219  # m
    k0: float = 7.2 * 1e10  # 1/min
    EbR: float = 8750  # K
    U: float = 54.94  # kJ / (min*m^2*K)
    rho: float = 1000  # kg / m^3
    Cp: float = 0.239  # kJ / (kg*K)
    dH: float = -5 * 1e4  # kJ / kmol
    # to avoid division by zero
    eps: float = 1e-5  # m
    xs: np.ndarray = np.array([0.878, 324.5, 0.659])
    us: np.ndarray = np.array([300, 0.1])


def setup_cstr_model(params: CSTRParameters):

    model_name = "cstr_ode"

    # set up states
    c = SX.sym("c")  # molar concentration of species A
    T = SX.sym("T")  # reactor temperature
    h = SX.sym("h")  # level of the tank

    x = vertcat(c, T, h)

    # controls
    Tc = SX.sym("Tc")  # temperature of coolant liquid
    F = SX.sym("F")  # outlet flowrate

    u = vertcat(Tc, F)

    # xdot
    c_dot = SX.sym("c_dot")
    T_dot = SX.sym("T_dot")
    h_dot = SX.sym("h_dot")

    xdot = vertcat(c_dot, T_dot, h_dot)

    # parameters
    F0 = SX.sym("F0")  # inlet flowrate
    p = F0

    # dynamics
    A_const = np.pi * params.r**2
    denom = A_const * (h + params.eps)
    k = params.k0 * exp(-params.EbR / T)
    rate = k * c

    f_expl = vertcat(
        F0 * (params.c0 - c) / denom - rate,
        F0 * (params.T0 - T) / denom
        - params.dH / (params.rho * params.Cp) * rate
        + 2 * params.U / (params.r * params.rho * params.Cp) * (Tc - T),
        (F0 - F) / A_const,
    )

    f_impl = xdot - f_expl

    model = AcadosModel()

    model.f_impl_expr = f_impl
    model.f_expl_expr = f_expl
    model.x = x
    model.xdot = xdot
    model.u = u
    model.p = p
    model.name = model_name

    return model


def setup_linearized_model(model, cstr_params, mpc_params):
    # linearized dynamics
    A, B, P = compute_lqr_gain(model, cstr_params, mpc_params)
    f_discrete = (
        cstr_params.xs + A @ (model.x - cstr_params.xs) + B @ (model.u - cstr_params.us)
    )

    model.disc_dyn_expr = f_discrete

    # print(f"P_mat {P}")
    return model
