#
# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# Licensed under the 2-Clause BSD License.

import os, sys
sys.path.insert(0, '../pendulum_on_cart/common')


from acados_template import AcadosModel, AcadosSim, AcadosSimSolver, casadi_length
import numpy as np
import casadi as ca
from utils import plot_pendulum

A_1 = 1
A_2 = 2
N = 2

def export_linear_time_var_model() -> AcadosModel:

    model_name = "linear_time_var"

    # set up states & controls
    x = ca.SX.sym("x")

    nx = casadi_length(x)

    xdot = ca.SX.sym('xdot', nx)

    t = ca.SX.sym("t")

    f_expl = ca.if_else(t < 1.0, A_1@x, A_2@x)

    f_impl = xdot - f_expl

    model = AcadosModel()

    model.f_impl_expr = f_impl
    model.f_expl_expr = f_expl
    model.x = x
    model.xdot = xdot
    model.u = ca.SX.sym('u', 0)
    # model.z = z
    # model.p = p
    model.name = model_name
    model.t = t

    return model


def main():
    sim = AcadosSim()

    # export model
    model = export_linear_time_var_model()

    # set model_name
    sim.model = model

    deltaT = 1.0
    nx = model.x.rows()

    # set simulation time
    sim.solver_options.T = deltaT
    # set options
    sim.solver_options.integrator_type = 'IRK'
    sim.solver_options.num_stages = 4
    sim.solver_options.num_steps = 3
    sim.solver_options.newton_iter = 20 # for implicit integrator
    sim.solver_options.newton_tol = 1e-14 # for implicit integrator
    sim.solver_options.collocation_type = "GAUSS_LEGENDRE"

    # create
    acados_integrator = AcadosSimSolver(sim)

    simX = np.zeros((N+1, nx))
    x0 = np.array([1.0])

    simX[0,:] = x0

    for i in range(N):
        # set initial state
        acados_integrator.set("x", simX[i,:])
        # update initial time
        acados_integrator.set("t0", i * deltaT)
        # initialize IRK
        if sim.solver_options.integrator_type == 'IRK':
            acados_integrator.set("xdot", np.zeros((nx,)))

        # solve
        status = acados_integrator.solve()
        # get solution
        simX[i+1,:] = acados_integrator.get("x")

        if status != 0:
            raise Exception(f'acados returned status {status}.')

        S_forw = acados_integrator.get("S_forw")

        if i == 0:
            true_lin_dyn = A_1
        else:
            true_lin_dyn = A_2

        # print("S_forw", S_forw, np.exp(true_lin_dyn))
        err_solution_sens = np.exp(true_lin_dyn) - S_forw
        assert(err_solution_sens < 1e-7)
    print(f"simulated time varying system accurately")


if __name__ == "__main__":
    main()

