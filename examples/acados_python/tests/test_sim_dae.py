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


    print("Success: z and z_analytic match!")
