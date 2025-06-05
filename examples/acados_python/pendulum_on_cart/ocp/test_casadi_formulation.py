
import sys
sys.path.insert(0, '../common')

import numpy as np
import casadi as ca
from typing import Union

from acados_template import AcadosOcp, AcadosOcpSolver, AcadosCasadiOcpSolver
from ocp_example_cost_formulations import formulate_ocp, T_HORIZON, FMAX

from utils import plot_pendulum

def raise_test_failure_message(msg: str):
    # print(f"ERROR: {msg}")
    raise Exception(msg)

def main(cost_version="CONL", constraint_version='h', casadi_solver_name="ipopt"):
    ocp = formulate_ocp(cost_version=cost_version, constraint_version=constraint_version)
    ocp.solver_options.tf = T_HORIZON

    initial_iterate = ocp.create_default_initial_iterate()

    ## solve using acados
    # create acados solver
    ocp_solver = AcadosOcpSolver(ocp, verbose=False)
    # initialize solver
    ocp_solver.load_iterate_from_obj(initial_iterate)
    # solve with acados
    status = ocp_solver.solve()
    ocp_solver.print_statistics()
    if status != 0:
        raise_test_failure_message(f"acados solver returned status {status}.")
    # get solution
    result_acados = ocp_solver.store_iterate_to_obj()
    result_acados_flat = ocp_solver.store_iterate_to_flat_obj()

    ## solve using casadi
    casadi_nlp_opts = {}
    if casadi_solver_name == "fatrop":
        casadi_nlp_opts["expand"] = True
        casadi_nlp_opts["fatrop"] = {"mu_init": 0.1}
        casadi_nlp_opts["structure_detection"] = "auto"
        casadi_nlp_opts["debug"] = True
    casadi_ocp_solver = AcadosCasadiOcpSolver(ocp, verbose=False, solver=casadi_solver_name, casadi_nlp_opts=casadi_nlp_opts)
    # casadi_ocp_solver.load_iterate_from_obj(result_acados)
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
            raise_test_failure_message(f"Test failed: difference between acados and casadi solution for {field} should be smaller than {test_tol}, but got {diff}.")

    # plot_pendulum(ocp.solver_options.shooting_nodes, FMAX,
    #               np.array(result_casadi.u_traj), np.array(result_casadi.x_traj), latexify=False)
    del ocp_solver, casadi_ocp_solver


if __name__ == "__main__":
    main(cost_version="CONL", constraint_version='bu')
    main(cost_version="CONL", constraint_version='h')
    main(cost_version="CONL", constraint_version='bu', casadi_solver_name="fatrop")
    # main(cost_version="CONL", constraint_version='h', casadi_solver_name="fatrop")
