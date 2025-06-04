
import sys
sys.path.insert(0, '../common')

import numpy as np
import casadi as ca
from typing import Union

from acados_template import AcadosOcp, AcadosOcpSolver, AcadosCasadiOcpSolver
from ocp_example_cost_formulations import formulate_ocp, T_HORIZON

from utils import plot_pendulum

def main():
    ocp = formulate_ocp("CONL")
    ocp.solver_options.tf = T_HORIZON

    initial_iterate = ocp.create_default_initial_iterate()

    ## solve using acados
    # create acados solver
    ocp_solver = AcadosOcpSolver(ocp, verbose=False)
    # initialize solver
    ocp_solver.load_iterate_from_obj(initial_iterate)
    # solve with acados
    status = ocp_solver.solve()
    # get solution
    result_acados = ocp_solver.store_iterate_to_obj()
    result_acados_flat = ocp_solver.store_iterate_to_flat_obj()

    ## solve using casadi
    casadi_ocp_solver = AcadosCasadiOcpSolver(ocp, verbose=False)
    casadi_ocp_solver.load_iterate_from_obj(result_acados)
    casadi_ocp_solver.solve()
    result_casadi = casadi_ocp_solver.store_iterate_to_obj()
    result_casadi_flat = casadi_ocp_solver.store_iterate_to_flat_obj()

    # evaluate difference
    # TODO: set solver tolerance and reduce it here.
    test_tol = 5e-5
    for field in ["x", "u", "pi", "lam"]:
        diff = np.linalg.norm(getattr(result_acados_flat, field) - getattr(result_casadi_flat, field))
        print(f"Difference between acados and casadi solution for {field:4}: {diff:.3e}")
        if diff > test_tol:
            raise ValueError(f"Test failed: difference between acados and casadi solution for {field} should be smaller than {test_tol}, but got {diff}.")

    # plot_pendulum(ocp.solver_options.shooting_nodes, ocp.constraints.ubu, u_casadi_sol, x_casadi_sol, latexify=False)


if __name__ == "__main__":
    main()