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
from scipy.linalg import block_diag
from math import isclose
latexify_plot()

def export_parametric_ocp() -> AcadosOcp:
    nx = 4
    nu = 2
    model = AcadosModel()
    ocp = AcadosOcp()
    model.x = ca.SX.sym("x", nx)
    model.u = ca.SX.sym("u", nu)

    x_next = ca.vertcat(model.x)
    x_next[:nu] += model.u
    model.disc_dyn_expr = x_next

    # define cost
    p_W = ca.SX.sym("W", nx+nu, nx+nu)
    p_W_e = ca.SX.sym("W_e", nx, nx)
    model.p = ca.vertcat(p_W[:], p_W_e[:])
    ocp.parameter_values = np.ones((model.p.size()[0], ))
    ny = nx+nu
    xu = ca.vertcat(model.x, model.u)
    model.cost_expr_ext_cost = 0.5*xu.T @ p_W @ xu
    model.cost_expr_ext_cost_e = 0.5*model.x.T @ p_W_e @ model.x
    ocp.cost.cost_type = "EXTERNAL"
    ocp.cost.cost_type_e = "EXTERNAL"

    model.name = "non_ocp"
    ocp.model = model

    ocp.constraints.x0 = np.ones((nx, ))

    ocp.solver_options.integrator_type = "DISCRETE"
    ocp.solver_options.qp_solver = "FULL_CONDENSING_HPIPM"
    ocp.solver_options.hessian_approx = "EXACT"
    ocp.solver_options.regularize_method = "MIRROR"
    ocp.solver_options.nlp_solver_type = "SQP"
    ocp.solver_options.N_horizon = 1
    ocp.solver_options.tf = 1.0
    ocp.solver_options.print_level = 1
    ocp.solver_options.nlp_solver_ext_qp_res = 1
    ocp.solver_options.nlp_solver_max_iter = 2
    ocp.solver_options.eval_residual_at_max_iter = False
    ocp.solver_options.reg_adaptive_eps = True
    ocp.solver_options.reg_max_cond_block = 1e3

    return ocp

def set_cost_matrix(solver, W_mat, W_mat_e):
    p_val = np.concatenate([W_mat.flatten(), W_mat_e.flatten()])
    solver.set(0, "p", p_val)
    solver.set(1, "p", p_val)

def test_reg_adaptive_eps(regularize_method='MIRROR'):
    ocp = export_parametric_ocp()
    ocp.solver_options.qp_solver_t0_init = 0
    ocp.solver_options.nlp_solver_max_iter = 2 # QP should converge in one iteration
    ocp.solver_options.regularize_method = regularize_method

    ocp_solver = AcadosOcpSolver(ocp, json_file="parameter_augmented_acados_ocp.json", verbose=False)

    nx = ocp.dims.nx
    nu = ocp.dims.nu

    W_mat2 = np.zeros((nx+nu, nx+nu))
    W_mat2[0,0] = 1e6
    W_mat2[nx+nu-1, nx+nu-1] = 1e-4

    p = np.pi/2
    A_u = np.array([[np.cos(p), -np.sin(p)], [np.sin(p), np.cos(p)]])
    mat_u = A_u.T @ np.diag([-1, -1e-3]) @ A_u

    # cf. https://stackoverflow.com/questions/65190660/orthogonality-of-a-4x4-matrix
    A_x = np.array([[0.5000,   0.5000,   0.5000,   0.5000],
                    [0.6533,   0.2706,  -0.2706,  -0.6533],
                    [0.5000,  -0.5000,  -0.5000,   0.5000],
                    [0.2706,  -0.6533,   0.6533,  -0.2706]])

    mat_x = A_x.T @ np.diag([15, 4.0, -2e5, 1e-6]) @ A_x

    W_mat3 = block_diag(mat_x, mat_u)

    W_mats = [np.zeros((nx+nu, nx+nu)), W_mat2, W_mat3]
    W_mat_e = np.zeros((nx, nx))

    for i, W_mat in enumerate(W_mats):
        print(f"{regularize_method} i={i}")
        print("---------------------")

        # Test zero matrix
        set_cost_matrix(ocp_solver, W_mat, W_mat_e)

        status = ocp_solver.solve()
        ocp_solver.print_statistics()
        nlp_iter = ocp_solver.get_stats("nlp_iter")

        qp_diagnostics = ocp_solver.qp_diagnostics()
        assert qp_diagnostics['condition_number_stage'][0] <= ocp.solver_options.reg_max_cond_block +1e-8, f"Condition number must be <= {ocp.solver_options.reg_max_cond_block} per stage, got {qp_diagnostics['condition_number_stage'][0]} for i = {i}"
        assert qp_diagnostics['condition_number_stage'][1] <= ocp.solver_options.reg_max_cond_block +1e-8, f"Condition number must be <= {ocp.solver_options.reg_max_cond_block} per stage, got {qp_diagnostics['condition_number_stage'][1]} for i = {i}"

        assert nlp_iter == 1, f"Number of NLP iterations should be 1, got {nlp_iter} for i = {i}"
        assert status == 0, f"acados returned status {status} for i = {i}"

        hessian_0 = ocp_solver.get_hessian_block(0)
        if i == 0:
            assert np.equal(hessian_0, np.eye(nx+nu)).all(), f"Zero Hessian matrix should be transformed into identity for {regularize_method} for i = {i}"
        elif i == 1:
            assert np.equal(hessian_0, np.diag([1e3, 1e3, 1e6, 1e3, 1e3, 1e3])).all(), f"Something in adaptive {regularize_method} went wrong for i = {i}!"
        elif i == 2:
            print(np.linalg.eigvals(W_mat))
            print(hessian_0)
            print(np.real(np.linalg.eigvals(hessian_0)))
            if regularize_method == 'MIRROR':
                reg_eps = 2e5/ocp.solver_options.reg_max_cond_block
                assert np.allclose(np.linalg.eigvals(hessian_0), np.array([reg_eps, 2e5, reg_eps, reg_eps, reg_eps, reg_eps])), f"Something in adaptive {regularize_method} went wrong for i = {i}!"
            elif regularize_method == 'PROJECT':
                reg_eps = 15/ocp.solver_options.reg_max_cond_block
                assert np.allclose(np.real(np.linalg.eigvals(hessian_0)), np.array([15, 4, reg_eps, reg_eps, reg_eps, reg_eps]), rtol=1e-03, atol=1e-3), f"Something in adaptive {regularize_method} went wrong for i = {i}!"

        hessian_1 = ocp_solver.get_hessian_block(1)
        assert np.equal(hessian_1, np.eye(nx)).all(), f"Zero Hessian matrix should be transformed into identity for {regularize_method} for i = {i}"


if __name__ == "__main__":
    test_reg_adaptive_eps("MIRROR")
    test_reg_adaptive_eps("PROJECT")
