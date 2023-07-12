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

from acados_template import AcadosOcp, AcadosOcpSolver
from scipy.linalg import block_diag
import numpy as np
from dataclasses import dataclass
from casadi import vertcat


@dataclass
class MpcCSTRParameters:
    umin: np.ndarray  # lower bound on u
    umax: np.ndarray  # upper bound on u
    Q: np.ndarray
    R: np.ndarray
    Tf: float = 0.25 * 15  # horizon length
    N: int = 15
    dt: float = 0.25
    linear_mpc: bool = False

    # NOTE: computed with setup_linearized_model()
    P: np.ndarray = np.array(
        [
            [5.92981953e-01, -8.40033347e-04, -1.54536980e-02],
            [-8.40033347e-04, 7.75225208e-06, 2.30677411e-05],
            [-1.54536980e-02, 2.30677411e-05, 2.59450075e00],
        ]
    )

    def __init__(self, xs, us, dt=0.25, linear_mpc=False, N=16, Tf=4):
        self.Q = np.diag(1.0 / xs**2)
        self.R = np.diag(1.0 / us**2)
        # from slide
        # self.umin = np.array([0.975, 0.75]) * us
        # self.umax = np.array([1.025, 1.25]) * us
        # from figure code
        self.umin = np.array([0.95, 0.85]) * us
        self.umax = np.array([1.05, 1.15]) * us


def setup_acados_ocp_solver(
    model, mpc_params: MpcCSTRParameters, cstr_params, use_rti=False
):

    ocp = AcadosOcp()

    # set model
    ocp.model = model
    x = model.x
    u = model.u
    nx = x.shape[0]
    nu = u.shape[0]

    # number of shooting intervals
    ocp.dims.N = mpc_params.N

    # set prediction horizon
    ocp.solver_options.tf = mpc_params.Tf

    # nominal parameter values
    ocp.parameter_values = np.array([cstr_params.F0])

    # set cost
    ocp.cost.W_e = mpc_params.P
    ocp.cost.W = block_diag(mpc_params.Q, mpc_params.R)

    ocp.cost.cost_type = "NONLINEAR_LS"
    ocp.cost.cost_type_e = "NONLINEAR_LS"

    ocp.model.cost_y_expr = vertcat(x, u)
    ocp.model.cost_y_expr_e = x

    ocp.cost.yref = np.zeros((nx + nu,))
    ocp.cost.yref_e = np.zeros((nx,))

    # set constraints
    ocp.constraints.lbu = mpc_params.umin
    ocp.constraints.ubu = mpc_params.umax
    ocp.constraints.idxbu = np.arange(nu)

    ocp.constraints.x0 = cstr_params.xs

    # set options
    ocp.solver_options.qp_solver = "PARTIAL_CONDENSING_HPIPM"  # FULL_CONDENSING_QPOASES
    # PARTIAL_CONDENSING_HPIPM, FULL_CONDENSING_QPOASES, FULL_CONDENSING_HPIPM,
    # PARTIAL_CONDENSING_QPDUNES, PARTIAL_CONDENSING_OSQP
    # ocp.solver_options.qp_solver = 'FULL_CONDENSING_QPOASES'
    ocp.solver_options.qp_solver_cond_N = mpc_params.N  # for partial condensing

    ocp.solver_options.hessian_approx = "GAUSS_NEWTON"
    # ocp.solver_options.print_level = 1
    if use_rti:
        ocp.solver_options.nlp_solver_type = "SQP_RTI"  # SQP_RTI, SQP
    else:
        ocp.solver_options.nlp_solver_type = "SQP"  # SQP_RTI, SQP

    if mpc_params.linear_mpc:
        ocp.solver_options.integrator_type = "DISCRETE"
    else:
        ocp.solver_options.integrator_type = "IRK"
        ocp.solver_options.sim_method_num_stages = 4
        ocp.solver_options.sim_method_num_steps = 1  # 5

    ocp.solver_options.levenberg_marquardt = 1e-5
    # ocp.solver_options.tol = 1e-3
    ocp.solver_options.line_search_use_sufficient_descent

    # create
    ocp_solver = AcadosOcpSolver(ocp, json_file="acados_ocp.json")

    return ocp_solver


if __name__ == "__main__":
    setup_acados_ocp_solver()
