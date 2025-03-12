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

import numpy as np
import casadi as ca
from acados_template import AcadosOcp, AcadosModel, AcadosOcpSolver, latexify_plot
import matplotlib.pyplot as plt
latexify_plot()

def export_parametric_ocp() -> AcadosOcp:
    nx = 4
    nu = 2
    model = AcadosModel()
    ocp = AcadosOcp()
    model.x = ca.SX.sym("x", nx)
    model.u = ca.SX.sym("u", nu)

    x_next = model.x
    x_next[:nu] += model.u
    model.disc_dyn_expr = x_next

    # define cost
    ny = nx+nu
    Vx = np.zeros((ny, nx))
    Vx[:nx, :] = np.eye(nx)
    Vu = np.zeros((ny, nu))
    Vu[nx:, :] = np.eye(nu)
    Vx_e = np.eye(nx)

    ocp.cost.Vx = Vx
    ocp.cost.Vx_e = Vx_e
    ocp.cost.Vu = Vu
    ocp.cost.cost_type = "LINEAR_LS"
    ocp.cost.cost_type_e = "LINEAR_LS"
    ocp.cost.W = np.eye(ny)
    ocp.cost.W_e = np.eye(nx)
    ocp.cost.yref = np.zeros((ny, ))
    ocp.cost.yref_e = np.zeros((nx, ))

    model.name = "non_ocp"
    ocp.model = model

    ocp.constraints.lbx_0 = np.array([-1.0])
    ocp.constraints.ubx_0 = np.array([1.0])
    ocp.constraints.idxbx_0 = np.array([0])

    ocp.solver_options.integrator_type = "DISCRETE"
    ocp.solver_options.qp_solver = "FULL_CONDENSING_HPIPM"
    ocp.solver_options.hessian_approx = "EXACT"
    ocp.solver_options.N_horizon = 1
    ocp.solver_options.tf = 1.0

    ocp.p_global_values = np.zeros((1,))
    ocp.solver_options.with_solution_sens_wrt_params = True
    ocp.solver_options.with_value_sens_wrt_params = True
    ocp.solver_options.nlp_solver_ext_qp_res = 1

    return ocp

def test_reg():

    ocp = export_parametric_ocp()
    ocp.solver_options.qp_solver_t0_init = 0
    ocp.solver_options.nlp_solver_max_iter = 2 # QP should converge in one iteration
    # TODO: set regularization options

    ocp_solver = AcadosOcpSolver(ocp, json_file="parameter_augmented_acados_ocp.json", verbose=False)

    # TODO: test different matrices
    W_mat = np.zeros()
    # ...


if __name__ == "__main__":
    test_reg()
