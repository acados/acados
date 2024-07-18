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

from acados_template import AcadosOcp, AcadosOcpSolver, AcadosModel
import numpy as np
import casadi as ca
import scipy.linalg

def create_acados_solver_and_solve_problem(globalization='FIXED_STEP', solver_type='SQP', eval_residual_at_max_iter=None):
    # create ocp object to formulate the OCP
    ocp = AcadosOcp()

    # parameters
    Tf = 0.15
    N = 3
    nx = 2
    nu = 2
    ny = nx + nu
    ny_e = nx
    x = ca.SX.sym("x", nx)
    u = ca.SX.sym("u", nu)

    A = np.eye(nx)*0.9
    B = np.tril(np.ones((nx, nu)))

    # set model
    model = AcadosModel()
    model.disc_dyn_expr = A @ x + B @ u
    model.x = x
    model.u = u
    model.name = "LinearModel"
    ocp.model = model

    # set dimensions
    ocp.dims.N = N

    # set cost
    ocp.cost.cost_type = "LINEAR_LS"
    ocp.cost.cost_type_e = "LINEAR_LS"
    Q_mat = 0.05*np.diag([1e3, 1e-2])
    R_mat = 0.05*np.diag([1.0, 1.0])

    ocp.cost.Vx = np.zeros((ny, nx))
    ocp.cost.Vx[:nx,:nx] = np.eye(nx)

    Vu = np.zeros((ny, nu))
    Vu[ny-1,0] = 1.0
    ocp.cost.Vu = Vu
    ocp.cost.Vx_e = np.eye(nx)

    ocp.cost.yref = np.zeros((ny, ))
    ocp.cost.yref_e = np.zeros((ny_e, ))
    ocp.cost.W = scipy.linalg.block_diag(Q_mat, R_mat)
    ocp.cost.W_e = Q_mat

    # set constraints
    ocp.constraints.x0 = np.array([2.0, 0.0])

    # set options
    ocp.solver_options.qp_solver = "PARTIAL_CONDENSING_HPIPM"
    ocp.solver_options.qp_solver_cond_N = N
    ocp.solver_options.hessian_approx = "EXACT"
    ocp.solver_options.integrator_type = "DISCRETE"
    ocp.solver_options.globalization = globalization
    ocp.solver_options.nlp_solver_max_iter = 1
    ocp.solver_options.print_level = 1
    if eval_residual_at_max_iter != None:
        ocp.solver_options.eval_residual_at_max_iter = eval_residual_at_max_iter
    ocp.solver_options.nlp_solver_type = solver_type

    # set prediction horizon
    ocp.solver_options.tf = Tf

    json_name = "acados_sqp_ocp.json"
    ocp_solver = AcadosOcpSolver(ocp, json_file = json_name)

    sol_X = np.zeros((N+1, nx))
    sol_U = np.zeros((N, nu))

    status = ocp_solver.solve()

    # get solution
    for i in range(N):
        sol_X[i,:] = ocp_solver.get(i, "x")
        sol_U[i,:] = ocp_solver.get(i, "u")
    sol_X[N,:] = ocp_solver.get(N, "x")

    return status, sol_X, sol_U

def main():
    # do the test with SQP
    funnel_status, funnel_X, funnel_U = create_acados_solver_and_solve_problem(globalization='FUNNEL_L1PEN_LINESEARCH')
    assert funnel_status == 0, "Funnel could not find solution after 1 iteration!"
    print("Funnel finds solution")

    full_step_no_eval, fs_no_eval_X, fs_no_eval_U = create_acados_solver_and_solve_problem()
    assert np.allclose(funnel_X, fs_no_eval_X), "Funnel and max iter SQP version did not terminate at same X."
    assert np.allclose(funnel_U, fs_no_eval_U), "Funnel and max iter SQP version did not terminate at same U."
    assert full_step_no_eval == 2, "Full step SQP did not terminate with max iter!"
    print("Full step exits with max iter even though being at solution (as planned)")

    full_step_eval, fs_X, fs_U = create_acados_solver_and_solve_problem(eval_residual_at_max_iter=True)
    assert full_step_eval == 0, "Full step SQP did not find solution after 1 iteration!"
    assert np.allclose(funnel_X, fs_X), "Funnel and full step SQP version did not terminate at same X."
    assert np.allclose(funnel_U, fs_U), "Funnel and full step SQP version did not terminate at same U."
    print("Full step finds solution")

    merit_status, sol_X_sqp, sol_U_sqp = create_acados_solver_and_solve_problem(globalization='MERIT_BACKTRACKING', eval_residual_at_max_iter=True)
    assert merit_status == 0, "Merit function could not find solution after 1 iteration!"
    assert np.allclose(funnel_X, sol_X_sqp), "Funnel and merit function SQP version did not terminate at same X."
    assert np.allclose(funnel_U, sol_U_sqp), "Funnel and merit function SQP version did not terminate at same U."
    print("Merit function finds solution")

    # Do the test with DDP as well
    funnel_status, ddp_X, ddp_U = create_acados_solver_and_solve_problem(globalization='MERIT_BACKTRACKING', solver_type='DDP', eval_residual_at_max_iter=True)
    assert funnel_status == 0, "DDP MERIT_BACKTRACKING could not find solution after 1 iteration!"
    assert np.allclose(ddp_X, sol_X_sqp), "Funnel and DDP did not terminate at same X."
    assert np.allclose(ddp_U, sol_U_sqp), "Funnel and DDP did not terminate at same U."
    print("DDP merit backtracking finds solution")

    # Do the test with DDP as well
    funnel_status, no_eval_ddp_X, no_eval_ddp_U = create_acados_solver_and_solve_problem(globalization='MERIT_BACKTRACKING', solver_type='DDP', eval_residual_at_max_iter=False)
    assert funnel_status == 2, "DDP MERIT_BACKTRACKING without eval last iteration could not find solution after 1 iteration!"
    assert np.allclose(no_eval_ddp_X, sol_X_sqp), "Funnel and DDP without eval last iteration did not terminate at same X."
    assert np.allclose(no_eval_ddp_U, sol_U_sqp), "Funnel and DDP without eval last iteration did not terminate at same U."
    print("DDP merit backtracking exits with max iter even though being at solution (as planned)")

if __name__ == "__main__":
    main()
