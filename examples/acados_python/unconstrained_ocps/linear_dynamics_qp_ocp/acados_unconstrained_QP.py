#
# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# Licensed under the 2-Clause BSD License.

from acados_template import AcadosOcp, AcadosOcpSolver, AcadosModel
import numpy as np
import casadi as ca
import scipy.linalg

def create_acados_solver_and_solve_problem(method="SQP"):
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
    ocp.solver_options.N_horizon = N

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
    ocp.solver_options.print_level = 1
    ocp.solver_options.nlp_solver_type = method # SQP_RTI, SQP

    # set prediction horizon
    ocp.solver_options.tf = Tf

    json_name = "acados_" + method + "_ocp.json"
    ocp_solver = AcadosOcpSolver(ocp, json_file = json_name)

    status = ocp_solver.solve()
    iter = ocp_solver.get_stats("nlp_iter")

    assert iter in [0,1], "Solver should converge within 1 iteration!"

    if status != 0:
        raise Exception(f"acados returned status {status}.")

    # get solution
    sol = ocp_solver.get_flat_iterate()

    return sol.x, sol.u

def main():
    sol_X_sqp, sol_U_sqp = create_acados_solver_and_solve_problem(method="SQP")
    sol_X_ddp, sol_U_ddp = create_acados_solver_and_solve_problem(method="DDP")

    np.testing.assert_allclose(sol_X_ddp, sol_X_sqp, atol=1e-8, err_msg="solution x of ddp and sqp do not coincide")
    np.testing.assert_allclose(sol_U_ddp, sol_U_sqp, atol=1e-8, err_msg="solution u of ddp and sqp do not coincide")

    print("Experiment was succesful!")

if __name__ == "__main__":
    main()
