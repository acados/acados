# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# Licensed under the 2-Clause BSD License.


import numpy as np
from acados_template import AcadosOcpSolver, AcadosFatropOcpSolver
from ocp_example_cost_formulations import formulate_ocp, T_HORIZON

TEST_TOL = 5e-5


def raise_test_failure_message(msg: str):
    raise Exception(msg)


def main(cost_version="LS", constraint_version="bu", hessian_approx="EXACT", integrator_type="ERK"):
    ocp = formulate_ocp(cost_version=cost_version, constraint_version=constraint_version)
    ocp.solver_options.tf = T_HORIZON
    ocp.solver_options.integrator_type = integrator_type
    ocp.solver_options.hessian_approx = hessian_approx
    ocp.solver_options.fixed_hess = False
    ocp.solver_options.nlp_solver_type = "SQP"
    ocp.solver_options.globalization = "MERIT_BACKTRACKING"
    ocp.solver_options.regularize_method = "MIRROR"
    ocp.solver_options.nlp_solver_max_iter = 200
    for tol in ["nlp_solver_tol_stat", "nlp_solver_tol_eq", "nlp_solver_tol_ineq", "nlp_solver_tol_comp"]:
        setattr(ocp.solver_options, tol, 1e-9)
    case = f"{cost_version}_{constraint_version}_{hessian_approx}_{integrator_type}"
    ocp.code_gen_options.code_export_directory = f"c_generated_code_fatrop_test_{case}"
    ocp.code_gen_options.json_file = f"acados_ocp_fatrop_test_{case}.json"
    N = ocp.solver_options.N_horizon

    ## solve using acados
    ocp_solver = AcadosOcpSolver(ocp, verbose=False)
    for k in range(N + 1):
        ocp_solver.set(k, "x", ocp.constraints.x0)
    status = ocp_solver.solve()
    if status != 0:
        raise_test_failure_message(f"acados solver returned status {status}.")

    ## solve using fatrop, from the same initial guess
    ocp.code_gen_options.code_export_directory += "_fatrop"
    ocp.code_gen_options.json_file = ocp.code_gen_options.json_file.replace(".json", "_fatrop.json")
    fatrop_solver = AcadosFatropOcpSolver(ocp, verbose=False)
    for k in range(N + 1):
        fatrop_solver.set(k, "x", ocp.constraints.x0)
    status = fatrop_solver.solve()
    fatrop_solver.print_statistics()
    if status != 0:
        raise_test_failure_message(f"fatrop solver returned status {status}.")

    # evaluate difference
    for field in ["x", "u", "pi", "lam"]:
        stages = range(N) if field in ["u", "pi"] else range(N + 1)
        ref = np.concatenate([ocp_solver.get(k, field) for k in stages])
        sol = np.concatenate([fatrop_solver.get(k, field) for k in stages])
        diff = np.linalg.norm(ref - sol)
        print(f"Difference between acados and fatrop solution for {field:4}: {diff:.3e}")
        if diff > TEST_TOL:
            raise_test_failure_message(f"Test failed: difference for {field} should be smaller than {TEST_TOL}, but got {diff}.")
    if abs(ocp_solver.get_cost() - fatrop_solver.get_cost()) > TEST_TOL:
        raise_test_failure_message("Test failed: cost of acados and fatrop solution differ.")

    del ocp_solver, fatrop_solver


if __name__ == "__main__":
    for cost_version in ["LS", "NLS", "EXTERNAL"]:
        for constraint_version in ["bu", "h"]:
            main(cost_version=cost_version, constraint_version=constraint_version)
    main(cost_version="NLS", constraint_version="h", hessian_approx="GAUSS_NEWTON")
    main(cost_version="NLS", constraint_version="h", integrator_type="IRK")
