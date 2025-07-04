import sys
sys.path.insert(0, '../getting_started')
import numpy as np
import casadi as ca
from typing import Union

from acados_template import AcadosOcp, AcadosOcpSolver, AcadosCasadiOcpSolver
from pendulum_model import export_pendulum_ode_model


def get_x_u_traj(ocp_solver: Union[AcadosOcpSolver, AcadosCasadiOcpSolver], N_horizon: int):
    ocp = ocp_solver.acados_ocp if isinstance(ocp_solver, AcadosOcpSolver) else ocp_solver.ocp
    simX = np.zeros((N_horizon+1, ocp.dims.nx))
    simU = np.zeros((N_horizon, ocp.dims.nu))
    for i in range(N_horizon):
        simX[i,:] = ocp_solver.get(i, "x")
        simU[i,:] = ocp_solver.get(i, "u")
    simX[N_horizon,:] = ocp_solver.get(N_horizon, "x")

    return simX, simU

def formulate_ocp(Tf: float = 1.0, N: int = 20)-> AcadosOcp:
    # create ocp object to formulate the OCP
    ocp = AcadosOcp()

    # set model
    model = export_pendulum_ode_model()
    # add dummy control
    model.u = ca.vertcat(model.u, ca.SX.sym('dummy_u'))
    ocp.model = model

    # set h constraints
    ocp.model.con_h_expr_0 = ca.norm_2(model.x)
    ocp.constraints.lh_0 = np.array([0])
    ocp.constraints.uh_0 = np.array([3.16])

    nx = model.x.rows()
    nu = model.u.rows()
    
    # set prediction horizon
    ocp.solver_options.N_horizon = N
    ocp.solver_options.tf = Tf

    # cost matrices
    Q_mat = 2*np.diag([1e3, 1e3, 1e-2, 1e-2])
    R_mat = 2*np.diag([1e-2, 1e-2])

    # path cost
    ocp.cost.cost_type = 'NONLINEAR_LS'
    ocp.model.cost_y_expr = ca.vertcat(model.x, model.u)
    ocp.cost.yref = np.zeros((nx+nu,))
    ocp.cost.W = ca.diagcat(Q_mat, R_mat).full()

    # terminal cost
    ocp.cost.cost_type_e = 'NONLINEAR_LS'
    ocp.cost.yref_e = np.zeros((nx,))
    ocp.model.cost_y_expr_e = model.x
    ocp.cost.W_e = Q_mat

    # set constraints
    Fmax = 80
    ocp.constraints.lbu = np.array([-Fmax])
    ocp.constraints.ubu = np.array([+Fmax])
    ocp.constraints.idxbu = np.array([0])

    ocp.constraints.x0 = np.array([0, np.pi, 0, 0])  # initial state
    ocp.constraints.idxbx_0 = np.array([0, 1, 2, 3])

    # set partial bounds for state
    ocp.constraints.lbx = np.array([-10,-10])
    ocp.constraints.ubx = np.array([10,10])
    ocp.constraints.idxbx = np.array([0,3])

    # set options
    ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM' # FULL_CONDENSING_QPOASES
    ocp.solver_options.hessian_approx = 'GAUSS_NEWTON' # 'GAUSS_NEWTON', 'EXACT'
    ocp.solver_options.integrator_type = 'ERK'
    ocp.solver_options.nlp_solver_type = 'SQP' # SQP_RTI, SQP
    ocp.solver_options.globalization = 'MERIT_BACKTRACKING' # turns on globalization
    
    return ocp

def main():
    N_horizon = 20
    Tf = 1.0
    ocp = formulate_ocp(Tf, N_horizon)

    initial_iterate = ocp.create_default_initial_iterate()

    ## solve using acados
    # create acados solver
    ocp_solver = AcadosOcpSolver(ocp,verbose=False)
    ocp_solver.load_iterate_from_obj(initial_iterate)
    # solve with acados
    status = ocp_solver.solve()
    # get solution
    simX, simU = get_x_u_traj(ocp_solver, N_horizon)
    result = ocp_solver.store_iterate_to_obj()
    lam = ocp_solver.get_flat("lam")
    pi = ocp_solver.get_flat("pi")

    # ## solve using casadi
    casadi_ocp_solver = AcadosCasadiOcpSolver(ocp=ocp,solver="ipopt",verbose=False)
    casadi_ocp_solver.load_iterate_from_obj(result)
    casadi_ocp_solver.solve()
    x_casadi_sol, u_casadi_sol = get_x_u_traj(casadi_ocp_solver, N_horizon)
    lam_casadi = casadi_ocp_solver.get_flat("lam")
    pi_casadi = casadi_ocp_solver.get_flat("pi")

    # evaluate difference
    diff_x = np.linalg.norm(x_casadi_sol - simX)
    print(f"Difference between casadi and acados solution in x: {diff_x}")
    diff_u = np.linalg.norm(u_casadi_sol - simU)
    print(f"Difference between casadi and acados solution in u: {diff_u}")
    diff_lam = np.linalg.norm(lam_casadi - lam)
    print(f"Difference between casadi and acados solution in lam: {diff_lam}")
    diff_pi = np.linalg.norm(pi_casadi - pi)
    print(f"Difference between casadi and acados solution in pi: {diff_pi}")

    test_tol = 1e-4
    if diff_x > test_tol or diff_u > test_tol or diff_lam > test_tol or diff_pi > test_tol:
        raise ValueError(f"Test failed: difference between casadi and acados solution should be smaller than {test_tol}, but got {diff_x} and {diff_u} and {diff_lam} and {diff_pi}.")

if __name__ == "__main__":
    main()