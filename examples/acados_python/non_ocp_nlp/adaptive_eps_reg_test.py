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
    ocp.solver_options.regularize_method = "MIRROR"
    ocp.solver_options.nlp_solver_type = "SQP"
    ocp.solver_options.N_horizon = 1
    ocp.solver_options.tf = 1.0
    ocp.solver_options.print_level = 1
    ocp.solver_options.nlp_solver_ext_qp_res = 1
    ocp.solver_options.nlp_solver_max_iter = 1
    ocp.solver_options.eval_residual_at_max_iter = False
    ocp.solver_options.reg_adaptive_eps = True
    ocp.solver_options.reg_max_cond_block = 1e3

    return ocp

def test_reg_adaptive_mirror():

    ocp = export_parametric_ocp()
    ocp.solver_options.qp_solver_t0_init = 0
    ocp.solver_options.nlp_solver_max_iter = 2 # QP should converge in one iteration
    # TODO: set regularization options

    ocp_solver = AcadosOcpSolver(ocp, json_file="parameter_augmented_acados_ocp.json", verbose=False)

    nx = ocp.dims.nx
    nu = ocp.dims.nu

    # Test zero matrix
    W_mat = np.zeros((nx+nu, nx+nu))
    W_mat_e = np.zeros((nx, nx))
    ocp_solver.cost_set(0, 'W', W_mat)
    ocp_solver.cost_set(1, 'W', W_mat_e)

    _ = ocp_solver.solve()

    hessian_0 = ocp_solver.get_hessian_block(0)
    assert np.equal(hessian_0, np.eye(nx+nu)).all(), "Zero Hessian matrix should be transformed into identity"
    hessian_1 = ocp_solver.get_hessian_block(1)
    assert np.equal(hessian_1, np.eye(nx)).all(), "Zero Hessian matrix should be transformed into identity"
    qp_diagnostics = ocp_solver.qp_diagnostics()
    assert qp_diagnostics['condition_number_stage'][0] <= ocp.solver_options.reg_max_cond_block, "Condition number must be <= ocp.solver_options.reg_max_cond_block per stage"
    assert qp_diagnostics['condition_number_stage'][1] <= ocp.solver_options.reg_max_cond_block, "Condition number must be <= ocp.solver_options.reg_max_cond_block per stage"


    # Second test
    W_mat = np.zeros((nx+nu, nx+nu))
    W_mat[0,0] = 1e6
    W_mat[nx+nu-1, nx+nu-1] = 1e-4
    W_mat_e = np.zeros((nx, nx))
    ocp_solver.cost_set(0, 'W', W_mat)
    ocp_solver.cost_set(1, 'W', W_mat_e)
    _ = ocp_solver.solve()

    hessian_0 = ocp_solver.get_hessian_block(0)
    assert np.equal(hessian_0, np.diag([1e3, 1e3, 1e6, 1e3, 1e3, 1e3])).all(), "Something in adaptive mirror went wrong!"
    qp_diagnostics = ocp_solver.qp_diagnostics()
    assert qp_diagnostics['condition_number_stage'][0] <= ocp.solver_options.reg_max_cond_block, "Condition number must be <= ocp.solver_options.reg_max_cond_block per stage"

    # Third test
    p = np.pi/2
    A_u = np.array([[np.cos(p), -np.sin(p)], [np.sin(p), np.cos(p)]])
    mat_u = A_u.T @ np.diag([-1, -1e-3]) @ A_u

    # cf. https://stackoverflow.com/questions/65190660/orthogonality-of-a-4x4-matrix
    A_x = np.array([[0.5000,   0.5000,   0.5000,   0.5000],
                    [0.6533,   0.2706,  -0.2706,  -0.6533],
                    [0.5000,  -0.5000,  -0.5000,   0.5000],
                    [0.2706,  -0.6533,   0.6533,  -0.2706]])

    mat_x = A_x.T @ np.diag([15, 4.0, -2e5, 1e-6]) @ A_x

    W_mat = block_diag(mat_x, mat_u)
    print(np.linalg.eigvals(W_mat))
    # print(W_mat)
    W_mat_e = np.zeros((nx, nx))
    ocp_solver.cost_set(0, 'W', W_mat)
    ocp_solver.cost_set(1, 'W', W_mat_e)
    _ = ocp_solver.solve()
    hessian_0 = ocp_solver.get_hessian_block(0)
    # print(hessian_0)
    # print(np.linalg.eigvals(hessian_0))
    assert np.allclose(np.linalg.eigvals(hessian_0), np.diag([2e5, 2e2, 2e2, 2e2, 2e2, 2e2])), "Something in adaptive mirror went wrong!"
    qp_diagnostics = ocp_solver.qp_diagnostics()
    assert qp_diagnostics['condition_number_stage'][0] <= ocp.solver_options.reg_max_cond_block, "Condition number must be <= ocp.solver_options.reg_max_cond_block per stage"


def test_reg_adaptive_project():

    ocp = export_parametric_ocp()
    ocp.solver_options.qp_solver_t0_init = 0
    ocp.solver_options.nlp_solver_max_iter = 2 # QP should converge in one iteration
    ocp.solver_options.regularize_method = 'PROJECT'
    # TODO: set regularization options

    ocp_solver = AcadosOcpSolver(ocp, json_file="parameter_augmented_acados_ocp.json", verbose=False)

    nx = ocp.dims.nx
    nu = ocp.dims.nu

    # Test zero matrix
    W_mat = np.zeros((nx+nu, nx+nu))
    W_mat_e = np.zeros((nx, nx))
    ocp_solver.cost_set(0, 'W', W_mat)
    ocp_solver.cost_set(1, 'W', W_mat_e)

    _ = ocp_solver.solve()

    hessian_0 = ocp_solver.get_hessian_block(0)
    assert np.equal(hessian_0, np.eye(nx+nu)).all(), "Zero Hessian matrix should be transformed into identity"
    hessian_1 = ocp_solver.get_hessian_block(1)
    assert np.equal(hessian_1, np.eye(nx)).all(), "Zero Hessian matrix should be transformed into identity"
    qp_diagnostics = ocp_solver.qp_diagnostics()
    assert qp_diagnostics['condition_number_stage'][0] <= ocp.solver_options.reg_max_cond_block, "Condition number must be <= ocp.solver_options.reg_max_cond_block per stage"
    assert qp_diagnostics['condition_number_stage'][1] <= ocp.solver_options.reg_max_cond_block, "Condition number must be <= ocp.solver_options.reg_max_cond_block per stage"


    # Second test
    W_mat = np.zeros((nx+nu, nx+nu))
    W_mat[0,0] = 1e6
    W_mat[nx+nu-1, nx+nu-1] = 1e-4
    W_mat_e = np.zeros((nx, nx))
    ocp_solver.cost_set(0, 'W', W_mat)
    ocp_solver.cost_set(1, 'W', W_mat_e)
    _ = ocp_solver.solve()

    hessian_0 = ocp_solver.get_hessian_block(0)
    assert np.equal(hessian_0, np.diag([1e3, 1e3, 1e6, 1e3, 1e3, 1e3])).all(), "Something in adaptive mirror went wrong!"
    qp_diagnostics = ocp_solver.qp_diagnostics()
    assert qp_diagnostics['condition_number_stage'][0] <= ocp.solver_options.reg_max_cond_block, "Condition number must be <= ocp.solver_options.reg_max_cond_block per stage"

    # Third test
    p = np.pi/2
    A_u = np.array([[np.cos(p), -np.sin(p)], [np.sin(p), np.cos(p)]])
    mat_u = A_u.T @ np.diag([-1, -1e-3]) @ A_u

    # cf. https://stackoverflow.com/questions/65190660/orthogonality-of-a-4x4-matrix
    A_x = np.array([[0.5000,   0.5000,   0.5000,   0.5000],
                    [0.6533,   0.2706,  -0.2706,  -0.6533],
                    [0.5000,  -0.5000,  -0.5000,   0.5000],
                    [0.2706,  -0.6533,   0.6533,  -0.2706]])

    mat_x = A_x.T @ np.diag([15, 4.0, -2e5, 1e-6]) @ A_x

    W_mat = block_diag(mat_x, mat_u)
    # print(np.linalg.eigvals(W_mat))
    # print(W_mat)
    W_mat_e = np.zeros((nx, nx))
    ocp_solver.cost_set(0, 'W', W_mat)
    ocp_solver.cost_set(1, 'W', W_mat_e)
    _ = ocp_solver.solve()
    hessian_0 = ocp_solver.get_hessian_block(0)
    # print(hessian_0)
    # print(np.linalg.eigvals(hessian_0))
    reg_eps = 15/1e3
    assert np.allclose(np.linalg.eigvals(hessian_0), np.diag([reg_eps, 15, 4, reg_eps, reg_eps, reg_eps])), "Something in adaptive project went wrong!"
    qp_diagnostics = ocp_solver.qp_diagnostics()
    assert qp_diagnostics['condition_number_stage'][0] <= ocp.solver_options.reg_max_cond_block, "Condition number must be <= ocp.solver_options.reg_max_cond_block per stage"


if __name__ == "__main__":
    test_reg_adaptive_mirror()
    test_reg_adaptive_project()
