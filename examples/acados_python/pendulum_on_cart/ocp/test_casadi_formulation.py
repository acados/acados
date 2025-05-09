
import sys
sys.path.insert(0, '../common')

import numpy as np
import casadi as ca

from acados_template import AcadosOcp, AcadosOcpSolver
from ocp_example_cost_formulations import formulate_ocp, T_HORIZON

from utils import plot_pendulum

def main():
    ocp = formulate_ocp("CONL")
    ocp.solver_options.tf = T_HORIZON
    N_horizon = ocp.solver_options.N_horizon

    # create casadi formulation
    casadi_nlp, bounds = ocp.create_casadi_nlp_formulation()

    casadi_solver = ca.nlpsol("nlp_solver", 'ipopt', casadi_nlp)
                            #   {'ipopt': {'print_level': 1}, 'print_time': False})
    nlp_sol = casadi_solver(lbx=bounds['lbx'], ubx=bounds['ubx'], lbg=bounds['lbg'], ubg=bounds['ubg'])

    sol_w = nlp_sol['x']
    xtraj = sol_w[:ocp.dims.nx*(N_horizon+1)].reshape((ocp.dims.nx, N_horizon+1)).full()
    utraj = sol_w[ocp.dims.nx*(N_horizon+1):].reshape((ocp.dims.nu, N_horizon)).full()

    xtraj = xtraj.T
    utraj = utraj.T

    # solve with acados
    ocp_solver = AcadosOcpSolver(ocp, verbose=False)
    status = ocp_solver.solve()
    # get solution
    simX = np.zeros((N_horizon+1, ocp.dims.nx))
    simU = np.zeros((N_horizon, ocp.dims.nu))
    for i in range(N_horizon):
        simX[i,:] = ocp_solver.get(i, "x")
        simU[i,:] = ocp_solver.get(i, "u")
    simX[N_horizon,:] = ocp_solver.get(N_horizon, "x")

    # evaluate difference
    diff_x = np.linalg.norm(xtraj - simX)
    print(f"Difference between casadi and acados solution: {diff_x}")
    diff_u = np.linalg.norm(utraj - simU)
    print(f"Difference between casadi and acados solution: {diff_u}")

    test_tol = 1e-5
    if diff_x > test_tol or diff_u > test_tol:
        raise ValueError(f"Test failed: difference between casadi and acados solution should be smaller than {test_tol}, but got {diff_x} and {diff_u}.")

    plot_pendulum(ocp.solver_options.shooting_nodes, ocp.constraints.ubu, utraj, xtraj, latexify=False)


if __name__ == "__main__":
    main()