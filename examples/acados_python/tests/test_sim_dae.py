#
# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# Licensed under the 2-Clause BSD License.

import sys
sys.path.insert(0, '../pendulum_on_cart/common')

from acados_template import AcadosSim, AcadosSimSolver, AcadosModel
from pendulum_model import export_augmented_pendulum_model
from utils import plot_pendulum
import numpy as np
import casadi as ca


def setup_model():
    model = AcadosModel()
    model.x = ca.SX.sym("x")
    model.xdot = ca.SX.sym("xdot")
    model.u = ca.SX.sym("u")
    model.z = ca.SX.sym("z")
    model.f_impl_expr = ca.vertcat(model.xdot - model.u, model.z - model.xdot**2)
    model.name = "dae_nonlinear_z"
    return model

def test_dae(model, x0):
    sim = AcadosSim()

    # set model_name
    sim.model = model

    Tf = 0.1
    nx = model.x.rows()
    nu = model.u.rows()
    N = 200

    # set simulation time
    sim.solver_options.T = Tf
    # set options
    sim.solver_options.integrator_type = 'IRK'
    sim.solver_options.num_stages = 4
    sim.solver_options.num_steps = 3
    sim.solver_options.newton_iter = 3 # for implicit integrator

    sim.solver_options.sens_forw = True
    sim.solver_options.sens_adj = True
    sim.solver_options.sens_algebraic = True
    sim.solver_options.sens_hess = True
    sim.solver_options.output_z = True
    sim.solver_options.sens_algebraic = True
    sim.solver_options.sim_method_jac_reuse = True
    sim.name = model.name

    # create
    integrator = AcadosSimSolver(sim)

    simX = np.zeros((N+1, nx))

    u0_val = 3.0
    u0 = np.array([u0_val])

    # test setter
    integrator.set("u", 3)
    integrator.set("u", 3.0)
    integrator.set("u", u0)

    simX[0,:] = x0

    for i in range(N):
        # set initial state
        integrator.set("x", simX[i,:])
        # solve
        status = integrator.solve()
        # get solution
        simX[i+1,:] = integrator.get("x")

    if status != 0:
        raise Exception(f'acados returned status {status}.')

    S_algebraic = integrator.get("S_algebraic")
    print("S_algebraic (dz_dxu) = ", S_algebraic)

    z = integrator.get("z")
    print("z = ", z)
    last_x0 = simX[N-1,:]
    print(f"{last_x0 = }")

    S_forw = integrator.get("S_forw")
    print("S_forw, sensitivities of simulaition result wrt x,u:\n", S_forw)

    return last_x0, u0_val, z, S_algebraic


def test_external_cost_hessian_with_algebraic_variables():
    from acados_template import AcadosOcp, AcadosOcpSolver

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


if __name__ == "__main__":
    model_1 = export_augmented_pendulum_model()
    x0 = np.array([0.0, np.pi+1, 0.0, 0.0])
    last_x0, u0, z, S_algebraic = test_dae(model_1, x0)

    z_analytic = np.array([last_x0[0], u0**2])
    err_z = np.abs(z - z_analytic)
    if np.any(err_z > 1e-6):
        raise Exception(f'z and z_analytic should match! Difference is {err_z}')


    model_2 = setup_model()
    x0 = np.array([1.])
    _, u0, _, S_algebraic = test_dae(model_2, x0)

    S_algebraic_analytic = np.array([[0, 6]])

    err_S_algebraic = np.abs(S_algebraic - S_algebraic_analytic)
    if np.any(err_S_algebraic > 1e-6):
        raise Exception(f'z and z_analytic should match! Difference is {err_S_algebraic}')

    test_external_cost_hessian_with_algebraic_variables()

    print("Success: z and z_analytic match!")
