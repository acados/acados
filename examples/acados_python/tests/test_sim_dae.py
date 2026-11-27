#
# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# Licensed under the 2-Clause BSD License.


from test_paths import add_example_paths

add_example_paths()

import pytest

from acados_template import AcadosSim, AcadosSimSolver, AcadosModel
from pendulum_model import export_augmented_pendulum_model
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

def run_dae(model, x0, exact_z_output=False):
    sim = AcadosSim()

    # set model_name
    sim.model = model

    Tf = 0.1
    nx = model.x.rows()
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
    sim.solver_options.exact_z_output = exact_z_output
    sim.solver_options.sens_algebraic = True
    sim.solver_options.sim_method_jac_reuse = True
    sim.name = f"{model.name}_exact_z_{int(exact_z_output)}"

    # create
    integrator = AcadosSimSolver(sim, verbose=False)

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

    assert status == 0, f'acados returned status {status}.'

    S_algebraic = integrator.get("S_algebraic")
    print("S_algebraic (dz_dxu) = ", S_algebraic)

    z = integrator.get("z")
    print("z = ", z)
    last_x0 = simX[N-1,:]
    print(f"{last_x0 = }")

    S_forw = integrator.get("S_forw")
    print("S_forw, sensitivities of simulation result wrt x,u:\n", S_forw)

    return last_x0, u0_val, z, S_algebraic


@pytest.mark.parametrize("exact_z_output", [False, True])
def test_augmented_pendulum_dae(exact_z_output):
    model = export_augmented_pendulum_model()
    x0 = np.array([0.0, np.pi+1, 0.0, 0.0])
    last_x0, u0, z, _ = run_dae(model, x0, exact_z_output)
    z_analytic = np.array([last_x0[0], u0**2])
    tol = 1e-12 if exact_z_output else 1e-7
    np.testing.assert_allclose(z, z_analytic, rtol=0.0, atol=tol)


@pytest.mark.parametrize("exact_z_output", [False, True])
def test_nonlinear_algebraic_output(exact_z_output):
    model = setup_model()
    x0 = np.array([1.])
    _, u0, z, S_algebraic = run_dae(model, x0, exact_z_output)
    np.testing.assert_allclose(z, np.array([u0**2]), rtol=0.0, atol=1e-12)
    np.testing.assert_allclose(
        S_algebraic,
        np.array([[0, 2 * u0]]),
        rtol=0.0,
        atol=1e-12,
    )


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__]))
