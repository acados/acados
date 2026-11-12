#
# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# Licensed under the 2-Clause BSD License.

from acados_template import AcadosModel, AcadosOcp, AcadosOcpSolver

import casadi as ca
import numpy as np


def test_external_cost_hessian_with_algebraic_variables():
    model = AcadosModel()
    x = ca.SX.sym("x")
    xdot = ca.SX.sym("xdot")
    u = ca.SX.sym("u")
    z = ca.SX.sym("z")
    model.x = x
    model.xdot = xdot
    model.u = u
    model.z = z
    model.f_impl_expr = ca.vertcat(xdot, z - x - u**2)
    model.name = "external_cost_hessian_dae"

    ocp = AcadosOcp()
    ocp.model = model
    ocp.solver_options.N_horizon = 1
    ocp.solver_options.tf = 0.1
    ocp.solver_options.integrator_type = "IRK"
    ocp.solver_options.sim_method_num_stages = 1
    ocp.solver_options.sim_method_num_steps = 1
    ocp.solver_options.sim_method_newton_iter = 5
    ocp.solver_options.nlp_solver_type = "SQP_RTI"
    ocp.solver_options.qp_solver = "FULL_CONDENSING_HPIPM"
    ocp.solver_options.hessian_approx = "EXACT"
    ocp.solver_options.print_level = 0

    cost = (
        0.5 * (x**2 + u**2)
        + 0.7 * x * z
        - 0.4 * u * z
        + 0.5 * z**2
    )
    ocp.cost.cost_type = "EXTERNAL"
    ocp.cost.cost_type_e = "EXTERNAL"
    ocp.model.cost_expr_ext_cost = cost
    ocp.model.cost_expr_ext_cost_e = 0.5 * x**2
    ocp.constraints.x0 = np.array([0.7])

    solver = AcadosOcpSolver(ocp, verbose=False)
    x_linearization = np.array([0.7])
    u_linearization = np.array([0.4])
    solver.set(0, "x", x_linearization)
    solver.set(0, "u", u_linearization)
    solver.set(1, "x", x_linearization)

    status = solver.solve()
    if status != 0:
        raise RuntimeError(f"acados returned status {status}.")

    # acados stores the stage Hessian in HPIPM ordering [u, x].
    hessian_acados = solver.get_from_qp_in(0, "R")
    cross_hessian = solver.get_from_qp_in(0, "S")
    hessian_acados = np.block([
        [hessian_acados, cross_hessian],
        [cross_hessian.T, solver.get_from_qp_in(0, "Q")],
    ])

    ux = ca.vertcat(u, x)
    z_reduced = x + u**2
    hess_uxz = ca.hessian(cost, ca.vertcat(u, x, z))[0]
    hess_ux = hess_uxz[:2, :2]
    hess_z = hess_uxz[2:, 2:]
    hess_z_ux = hess_uxz[2:, :2]
    dzdux = ca.jacobian(z_reduced, ux)

    # This is the approximation implemented by acados. The second-order
    # sensitivity term sum_i l_z[i] * d2 z[i] / d(ux)^2 is omitted.
    hessian_reference = (
        hess_ux
        + dzdux.T @ hess_z_ux
        + hess_z_ux.T @ dzdux
        + dzdux.T @ hess_z @ dzdux
    )
    hessian_reference_fun = ca.Function(
        "hessian_reference",
        [x, u, z],
        [hessian_reference],
    )
    hessian_reference = hessian_reference_fun(
        x_linearization,
        u_linearization,
        np.array([x_linearization[0] + u_linearization[0]**2]),
    ).full()
    # External stage costs are scaled by the integration step.
    hessian_reference *= ocp.solver_options.tf / ocp.solver_options.N_horizon

    np.testing.assert_allclose(
        hessian_acados,
        hessian_reference,
        rtol=1e-10,
        atol=1e-10,
    )
    print("test_external_cost_hessian_with_algebraic_variables: passed.")


if __name__ == "__main__":
    test_external_cost_hessian_with_algebraic_variables()
