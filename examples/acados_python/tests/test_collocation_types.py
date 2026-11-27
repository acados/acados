#
# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# Licensed under the 2-Clause BSD License.


import pytest

from acados_template import AcadosSim, AcadosSimSolver
from export_chain_mass_model import export_chain_mass_model
from utils import get_chain_params
from pendulum_model import export_lifted_pendulum_model
import numpy as np

TEST_CONFIGS = [
    ("GAUSS_RADAU_IIA", 2),
    ("GAUSS_RADAU_IIA", 3),
    ("GAUSS_RADAU_IIA", 4),
    ("GAUSS_RADAU_IIA", 5),
    ("GAUSS_RADAU_IIA", 6),
    ("GAUSS_RADAU_IIA", 7),
    ("GAUSS_RADAU_IIA", 8),
    ("GAUSS_LEGENDRE", 2),
    ("GAUSS_LEGENDRE", 3),
    ("GAUSS_LEGENDRE", 4),
    ("GAUSS_LEGENDRE", 5),
    ("GAUSS_LEGENDRE", 6),
    ("GAUSS_LEGENDRE", 7),
    ("GAUSS_LEGENDRE", 8),
    ("SDIRK", 2),
    ("SDIRK", 3),
    ("SDIRK", 5),
    ("SDIRK", 9),
    ("SDIRK", 11),
    ("DIRK", 6),
    ("DIRK", 10),
]

def run_model(model, collocation_type, hessian, ns, num_steps, x0, jac_reuse):
    sim = AcadosSim()

    # set model_name
    sim.model = model

    Tf = 0.01
    nx = model.x.rows()
    nu = model.u.rows()

    # set simulation time
    sim.solver_options.T = Tf
    # set options
    sim.solver_options.integrator_type = 'IRK'
    sim.solver_options.collocation_type = collocation_type
    sim.solver_options.num_stages = ns
    sim.solver_options.num_steps = num_steps
    sim.solver_options.newton_iter = 10 # for implicit integrator
    sim.solver_options.newton_tol = 1e-10 # for implicit integrator

    sim.solver_options.sens_forw = True
    sim.solver_options.sens_adj = True
    sim.solver_options.sens_hess = hessian
    sim.solver_options.jac_reuse = jac_reuse
    sim.name = f"{model.name}_{collocation_type}_{ns}"

    # create
    acados_integrator = AcadosSimSolver(sim, verbose=False)

    u0 = np.zeros((nu, 1))

    S_adj = np.ones(nx)
    acados_integrator.set("seed_adj", S_adj)


    # set initial state
    acados_integrator.set("x", x0)
    # solve
    status = acados_integrator.solve()
    assert status == 0, f'acados returned status {status}.'
    # get solution
    x_end = acados_integrator.get("x")
    sens_forw = acados_integrator.get("S_forw")
    sens_adj = acados_integrator.get("S_adj")
    hess = acados_integrator.get("S_hess")

    return x_end, sens_forw, sens_adj, hess


@pytest.fixture(scope="module")
def chain_mass_case():
    n_mass = 3
    chain_params = get_chain_params()
    m = chain_params["m"]
    D = chain_params["D"]
    L = chain_params["L"]
    model = export_chain_mass_model(n_mass, m, D, L)
    nx = model.x.rows()
    nu = model.u.rows()
    M = n_mass - 2
    # position of last mass
    xPosFirstMass = np.zeros((3,1))
    xEndRef = np.zeros((3,1))
    xEndRef[0] = L * (M+1) * 6
    # initial state
    pos0_x = np.linspace(xPosFirstMass[0], xEndRef[0], n_mass)
    x0 = np.zeros((nx, 1))
    x0[:3*(M+1):3] = pos0_x[1:].reshape((M+1,1))

    reference = run_model(model, "GAUSS_RADAU_IIA", True, 9, 1000, x0, False)
    return model, x0, reference


@pytest.mark.parametrize("hessian", [True, False])
@pytest.mark.parametrize("jac_reuse", [True, False])
@pytest.mark.parametrize(("collocation_type", "ns"), TEST_CONFIGS)
def test_chain_mass_model(chain_mass_case, hessian, jac_reuse, collocation_type, ns):
    model, x0, reference = chain_mass_case
    x_true, sens_forw_true, sens_adj_true, hess_true = reference
    x, sens_forw, sens_adj, hess = run_model(
        model, collocation_type, hessian, ns, 40, x0, jac_reuse
    )
    np.testing.assert_allclose(x, x_true, rtol=0.0, atol=1e-8)
    np.testing.assert_allclose(sens_forw, sens_forw_true, rtol=0.0, atol=1e-8)
    np.testing.assert_allclose(sens_adj, sens_adj_true, rtol=0.0, atol=1e-8)
    if hessian:
        np.testing.assert_allclose(hess, hess_true, rtol=0.0, atol=1e-8)


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__]))


@pytest.fixture(scope="module")
def pendulum_case():
    model = export_lifted_pendulum_model()
    x0 = np.array([0.0, np.pi+1, 0.0, 0.0])
    reference = run_model(model, "GAUSS_RADAU_IIA", True, 9, 1000, x0, False)
    return model, x0, reference


@pytest.mark.parametrize("hessian", [True, False])
@pytest.mark.parametrize("jac_reuse", [True, False])
@pytest.mark.parametrize(("collocation_type", "ns"), TEST_CONFIGS)
def test_pendulum_model(pendulum_case, hessian, jac_reuse, collocation_type, ns):
    model, x0, reference = pendulum_case
    x_true, sens_forw_true, sens_adj_true, hess_true = reference
    x, sens_forw, sens_adj, hess = run_model(
        model, collocation_type, hessian, ns, 40, x0, jac_reuse
    )
    np.testing.assert_allclose(x, x_true, rtol=0.0, atol=1e-8)
    np.testing.assert_allclose(sens_forw, sens_forw_true, rtol=0.0, atol=1e-8)
    np.testing.assert_allclose(sens_adj, sens_adj_true, rtol=0.0, atol=1e-8)
    if hessian:
        np.testing.assert_allclose(hess, hess_true, rtol=0.0, atol=1e-8)
