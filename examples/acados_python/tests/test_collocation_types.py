#
# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# Licensed under the 2-Clause BSD License.

import sys
sys.path.insert(0, '../pendulum_on_cart/common')
sys.path.insert(0, '../chain_mass')

from acados_template import AcadosSim, AcadosSimSolver, AcadosModel
from export_chain_mass_model import export_chain_mass_model
from utils import get_chain_params
from pendulum_model import export_lifted_pendulum_model
import numpy as np
import casadi as ca

def run_model(model, collocation_type, hessian, ns, num_steps, x0):
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
    # get solution
    x_end = acados_integrator.get("x")
    sens_forw = acados_integrator.get("S_forw")
    sens_adj = acados_integrator.get("S_adj")
    hess = acados_integrator.get("S_hess")

    return x_end, sens_forw, sens_adj, hess


def test_chain_mass_model():
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

    x_true, sens_forw_true, sens_adj_true, hess_true = run_model(model, "GAUSS_RADAU_IIA", True, 9, 1000, x0)

    test_configs = [
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
    errors = False
    for hessian in (True,False):
        for (collocation_type, ns) in test_configs:
            x, sens_forw, sens_adj, hess = run_model(model, collocation_type, hessian, ns, 40, x0)
            if np.any(np.abs(x-x_true) >= 1e-8):
                errors = True
                print(f"ERROR: ODE integrator: {ns} stage {collocation_type} x error = {np.max(np.abs(x-x_true))}")
            if np.any(np.abs(sens_forw-sens_forw_true) >= 1e-8):
                errors = True
                print(f"ERROR: ODE integrator: {ns} stage {collocation_type} sens_forw error = {np.max(np.abs(sens_forw-sens_forw_true))}")
            if np.any(np.abs(sens_adj-sens_adj_true) >= 1e-8):
                errors = True
                print(f"ERROR: ODE integrator: {ns} stage {collocation_type} sens_adj error = {np.max(np.abs(sens_adj-sens_adj_true))}")
            if hessian and np.any(np.abs(hess-hess_true) >= 1e-8):
                errors = True
                print(f"ERROR: ODE integrator: {ns} stage {collocation_type} hess error = {np.max(np.abs(hess-hess_true))}")

    return errors

def test_pendulum_model():
    model = export_lifted_pendulum_model()
    x0 = np.array([0.0, np.pi+1, 0.0, 0.0])

    x_true, sens_forw_true, sens_adj_true, hess_true = run_model(model, "GAUSS_RADAU_IIA", True, 9, 1000, x0)

    test_configs = [
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
    errors = False
    for hessian in (True,False):
        for (collocation_type, ns) in test_configs:
            x, sens_forw, sens_adj, hess = run_model(model, collocation_type, hessian, ns, 40, x0)
            if np.any(np.abs(x-x_true) >= 1e-8):
                errors = True
                print(f"ERROR: DAE integrator: {ns} stage {collocation_type} x error = {np.max(np.abs(x-x_true))}")
            if np.any(np.abs(sens_forw-sens_forw_true) >= 1e-8):
                errors = True
                print(f"ERROR: DAE integrator: {ns} stage {collocation_type} sens_forw error = {np.max(np.abs(sens_forw-sens_forw_true))}")
            if np.any(np.abs(sens_adj-sens_adj_true) >= 1e-8):
                errors = True
                print(f"ERROR: DAE integrator: {ns} stage {collocation_type} sens_adj error = {np.max(np.abs(sens_adj-sens_adj_true))}")
            if hessian and np.any(np.abs(hess-hess_true) >= 1e-8):
                errors = True
                print(f"ERROR: DAE integrator: {ns} stage {collocation_type} hess error = {np.max(np.abs(hess-hess_true))}")

    return errors

if __name__ == "__main__":
    errors = test_pendulum_model()
    #errors |= test_chain_mass_model()

    if errors:
        raise Exception("Some errors reported above!")
