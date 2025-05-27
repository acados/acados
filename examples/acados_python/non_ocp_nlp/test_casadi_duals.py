
import sys
sys.path.insert(0, '../common')

import numpy as np
import casadi as ca
from typing import Union

from acados_template import AcadosOcp, AcadosOcpSolver, AcadosCasadiOcpSolver

from utils import plot_pendulum


from acados_template import AcadosModel

def formulate_ocp() -> AcadosOcp:

    model = AcadosModel()
    model.x = ca.SX.sym("x", 1)
    model.p_global = ca.SX.sym("p_global", 1)
    model.cost_expr_ext_cost_e = (model.x - model.p_global**2)**2
    model.name = "non_ocp"
    ocp = AcadosOcp()
    ocp.model = model

    ocp.constraints.lbx_e = np.array([-1.0])
    ocp.constraints.ubx_e = np.array([1.0])
    ocp.constraints.idxbx_e = np.array([0])

    ocp.cost.cost_type_e = "EXTERNAL"
    ocp.solver_options.qp_solver = "FULL_CONDENSING_HPIPM"
    ocp.solver_options.N_horizon = 0

    ocp.p_global_values = np.zeros((1,))

    return ocp

def get_x_u_traj(ocp_solver: Union[AcadosOcpSolver, AcadosCasadiOcpSolver], N_horizon: int):
    ocp = ocp_solver.acados_ocp
    simX = np.zeros((N_horizon+1, ocp.dims.nx))
    simU = np.zeros((N_horizon, ocp.dims.nu))
    for i in range(N_horizon):
        simX[i,:] = ocp_solver.get(i, "x")
        simU[i,:] = ocp_solver.get(i, "u")
    simX[N_horizon,:] = ocp_solver.get(N_horizon, "x")

    return simX, simU


def main():
    ocp = formulate_ocp()
    N_horizon = 0

    ## solve using casadi
    casadi_ocp_solver = AcadosCasadiOcpSolver(ocp, verbose=False)
    casadi_ocp_solver.solve()
    x_casadi_sol, u_casadi_sol = get_x_u_traj(casadi_ocp_solver, N_horizon)

    initial_iterate = ocp.create_default_initial_iterate()

    ## solve using acados
    # create acados solver
    ocp_solver = AcadosOcpSolver(ocp, verbose=False)
    # initialize solver
    ocp_solver.load_iterate_from_obj(initial_iterate)
    # solve with acados
    status = ocp_solver.solve()
    # get solution
    simX, simU = get_x_u_traj(ocp_solver, N_horizon)


    # evaluate difference
    diff_x = np.linalg.norm(x_casadi_sol - simX)
    print(f"Difference between casadi and acados solution: {diff_x}")
    diff_u = np.linalg.norm(u_casadi_sol - simU)
    print(f"Difference between casadi and acados solution: {diff_u}")

    # TODO: set solver tolerance and reduce it here.
    test_tol = 5e-5
    if diff_x > test_tol or diff_u > test_tol:
        raise ValueError(f"Test failed: difference between casadi and acados solution should be smaller than {test_tol}, but got {diff_x} and {diff_u}.")

    plot_pendulum(ocp.solver_options.shooting_nodes, ocp.constraints.ubu, u_casadi_sol, x_casadi_sol, latexify=False)


if __name__ == "__main__":
    main()