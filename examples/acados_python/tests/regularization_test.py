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

import sys
sys.path.insert(0, '../pendulum_on_cart/common')

from acados_template import AcadosOcp, AcadosOcpSolver, AcadosModel

import casadi as ca
import numpy as np


def create_linear_model() -> AcadosModel:
    model = AcadosModel()
    nx = 4
    model.x = ca.SX.sym('x', nx)
    model.u = ca.SX.sym('u', 4)
    model.name = 'linear_model'
    A_mat = np.eye(nx)
    A_mat += 0.1 * np.random.rand(nx, nx)
    model.disc_dyn_expr = A_mat @ model.x + 0.5 * model.u
    return model

def main(regularize_method: str):
    ocp = AcadosOcp()

    # set model
    model = create_linear_model()
    ocp.model = model

    Tf = 1.0
    nx = model.x.rows()

    ocp.solver_options.N_horizon = 1

    # set cost
    Q_mat = np.array([[8.332636379960917, -0.2025707437550449, 65.20466910278751, 0],
                      [-0.2025707437550449, 2.096166574366247e-05, -0.01903101665209866, 0],
                      [65.20466910278751, -0.01903101665209866, 22.29791833831988, 0],
                      [0, 0, 0, 2e-5],])
    # Q_mat = np.eye(nx)
    Q_mat_e = np.eye(nx) # Q_mat
    R_mat = Q_mat

    x = model.x
    u = model.u

    ocp.cost.cost_type = 'EXTERNAL'
    ocp.cost.cost_type_e = 'EXTERNAL'

    ocp.model.cost_expr_ext_cost = .5*x.T @ Q_mat @ x + .5*u.T @ R_mat @ u
    ocp.model.cost_expr_ext_cost_e = .5*x.T @ Q_mat_e @ x

    # set constraints
    Fmax = 80
    ocp.constraints.lbu = np.array([-Fmax])
    ocp.constraints.ubu = np.array([+Fmax])
    ocp.constraints.idxbu = np.array([0])

    ocp.constraints.x0 = np.ones((nx,))

    # options
    ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM'
    ocp.solver_options.hessian_approx = 'GAUSS_NEWTON'
    ocp.solver_options.regularize_method = regularize_method
    ocp.solver_options.integrator_type = 'DISCRETE'
    ocp.solver_options.nlp_solver_type = 'SQP'
    ocp.solver_options.nlp_solver_max_iter = 1000
    ocp.solver_options.globalization_full_step_dual = 0
    ocp.solver_options.qp_solver_iter_max = 2000
    ocp.solver_options.tol = 1e-4
    ocp.solver_options.qp_tol = 1e-6

    ocp.solver_options.tf = Tf

    # create solver
    ocp_solver = AcadosOcpSolver(ocp, verbose=False)
    _ = ocp_solver.solve()

    ocp_solver.print_statistics()

    # get solution
    sol = ocp_solver.get_iterate()

    Q_0_mat = ocp_solver.get_from_qp_in(0, "Q")
    hess_block = ocp_solver.get_hessian_block(0)

    # check symmetry
    np.testing.assert_allclose(hess_block, hess_block.T)

    # printing
    print(f"Solver with {regularize_method} returned \n{sol}")
    print(f"regularized Hessian is {hess_block}\n")

    # check eigenvalues
    if regularize_method == 'NO_REGULARIZE':
        print(np.max(np.abs(Q_0_mat - Q_mat)))
        print(f"Q_0_mat = {Q_0_mat}")
        print(f"Q_mat = {Q_mat}")
        np.testing.assert_allclose(Q_0_mat, Q_mat)
    else:
        eigvals = np.linalg.eigvals(hess_block)
        assert np.min(eigvals) > 0
        print(f"eigenvalues after regularization {eigvals}\n")

    # clear solver
    ocp_solver = None

    return sol


if __name__ == '__main__':
    main(regularize_method='NO_REGULARIZE')
    main(regularize_method='MIRROR')
    main(regularize_method='CONVEXIFY')
    main(regularize_method='PROJECT')
