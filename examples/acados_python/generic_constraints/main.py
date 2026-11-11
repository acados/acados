#
# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# Licensed under the 2-Clause BSD License.

import argparse

import casadi as ca
import numpy as np

from acados_template import AcadosOcp, AcadosOcpSolver, ocp_get_default_cmake_builder


def create_ocp(use_generic_constraints=True, use_cmake=False):
    nx = 2
    nu = 1
    horizon = 20

    x = ca.SX.sym("x", nx)
    u = ca.SX.sym("u", nu)

    ocp = AcadosOcp()
    ocp.model.name = "generic_constraints"
    ocp.model.x = x
    ocp.model.u = u
    ocp.model.disc_dyn_expr = ca.vertcat(
        0.9 * x[0] + 0.1 * x[1],
        0.1 * x[0] + 0.9 * x[1] + u[0],
    )

    ocp.cost.cost_type = "LINEAR_LS"
    ocp.cost.cost_type_e = "LINEAR_LS"
    ocp.cost.Vx = np.eye(nx + nu)[:, :nx]
    ocp.cost.Vu = np.eye(nx + nu)[:, nx:]
    ocp.cost.W = np.eye(nx + nu)
    ocp.cost.yref = np.zeros(nx + nu)
    ocp.cost.Vx_e = np.eye(nx)
    ocp.cost.W_e = np.eye(nx)
    ocp.cost.yref_e = np.zeros(nx)

    if use_generic_constraints:
        constraint_source = "generic_constraints.c"
        ocp.constraints.constr_h_ext_fun_type_0 = "generic"
        ocp.constraints.constr_h_generic_source_0 = constraint_source
        ocp.constraints.constr_h_fun_0 = "generic_constraints_fun"
        ocp.constraints.constr_h_fun_jac_0 = "generic_constraints_fun_jac"
        ocp.constraints.constr_h_fun_jac_hess_0 = "generic_constraints_fun_jac_hess"
        ocp.constraints.constr_h_ext_fun_type = "generic"
        ocp.constraints.constr_h_generic_source = constraint_source
        ocp.constraints.constr_h_fun = "generic_constraints_fun"
        ocp.constraints.constr_h_fun_jac = "generic_constraints_fun_jac"
        ocp.constraints.constr_h_fun_jac_hess = "generic_constraints_fun_jac_hess"
        ocp.constraints.constr_h_ext_fun_type_e = "generic"
        ocp.constraints.constr_h_generic_source_e = constraint_source
        ocp.constraints.constr_h_fun_e = "generic_constraints_fun_e"
        ocp.constraints.constr_h_fun_jac_e = "generic_constraints_fun_jac_e"
        ocp.constraints.constr_h_fun_jac_hess_e = "generic_constraints_fun_jac_hess_e"
    else:
        ocp.model.con_h_expr_0 = ca.vertcat(u, x)
        ocp.model.con_h_expr = ca.vertcat(u, x)
        ocp.model.con_h_expr_e = x

    ocp.constraints.x0 = np.array([1.0, 1.0])
    ocp.constraints.lh_0 = np.array([-0.5, -2.0, -2.0])
    ocp.constraints.uh_0 = np.array([0.5, 2.0, 2.0])
    ocp.constraints.lh = np.array([-0.5, -2.0, -2.0])
    ocp.constraints.uh = np.array([0.5, 2.0, 2.0])
    ocp.constraints.lh_e = np.array([-2.0, -2.0])
    ocp.constraints.uh_e = np.array([2.0, 2.0])

    ocp.solver_options.N_horizon = horizon
    ocp.solver_options.tf = 2.0
    ocp.solver_options.integrator_type = "DISCRETE"
    ocp.solver_options.hessian_approx = "EXACT"
    ocp.solver_options.nlp_solver_type = "SQP"
    ocp.solver_options.qp_solver = "PARTIAL_CONDENSING_HPIPM"
    ocp.solver_options.print_level = 0

    cmake_builder = ocp_get_default_cmake_builder() if use_cmake else None
    return ocp, cmake_builder


def solve(use_generic_constraints, use_cmake=False):
    ocp, cmake_builder = create_ocp(use_generic_constraints, use_cmake)
    solver = AcadosOcpSolver(ocp, cmake_builder=cmake_builder)
    status = solver.solve()

    if status != 0:
        raise RuntimeError(f"acados solver returned status {status}")
    return (
        np.array([solver.get(stage, "x") for stage in range(ocp.solver_options.N_horizon + 1)]),
        np.array([solver.get(stage, "u") for stage in range(ocp.solver_options.N_horizon)]),
    )


def main(use_cmake=False, compare_casadi=False):
    generic_solution = solve(use_generic_constraints=True, use_cmake=use_cmake)

    if compare_casadi:
        casadi_solution = solve(use_generic_constraints=False, use_cmake=use_cmake)
        np.testing.assert_allclose(generic_solution[0], casadi_solution[0], rtol=1e-10, atol=1e-10)
        np.testing.assert_allclose(generic_solution[1], casadi_solution[1], rtol=1e-10, atol=1e-10)
        print("generic and CasADi nonlinear constraints solutions are identical")
    else:
        print("generic nonlinear constraints example: success")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--cmake", action="store_true", help="build the generated solver with CMake")
    parser.add_argument("--compare-casadi", action="store_true",
                        help="also solve with CasADi constraints and compare the solutions")
    args = parser.parse_args()
    main(use_cmake=args.cmake, compare_casadi=args.compare_casadi)
