#
# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# Licensed under the 2-Clause BSD License.

import sys
sys.path.insert(0, '../common')

from acados_template import AcadosSim, AcadosSimSolver, sim_get_default_cmake_builder
from pendulum_model import export_pendulum_ode_model
# from utils import plot_pendulum
import numpy as np


def main(build=True, generate=True, use_cmake=True, use_cython=False):
    sim = AcadosSim()
    sim.model = export_pendulum_ode_model()

    Tf = 0.1
    nx = sim.model.x.rows()
    N = 200

    # set simulation time
    sim.solver_options.T = Tf
    # set options
    sim.solver_options.integrator_type = 'IRK'
    sim.solver_options.num_stages = 3
    sim.solver_options.num_steps = 3
    sim.solver_options.newton_iter = 3 # for implicit integrator
    sim.solver_options.collocation_type = "GAUSS_RADAU_IIA"

    cmake_builder = sim_get_default_cmake_builder() if use_cmake else None
    sim.name = "pendulum"

    # create
    if use_cython:
        acados_integrator = AcadosSimSolver.create_cython_solver(sim)
    else:
        sim.simulink_opts = dict()
        acados_integrator = AcadosSimSolver(sim, cmake_builder=cmake_builder, generate=generate, build=build)

    x0 = np.array([0.0, np.pi+1, 0.0, 0.0])
    u0 = np.array([0.0])

    simX = np.zeros((N+1, nx))
    simX[0,:] = x0

    for i in range(N):
        # initialize IRK
        if sim.solver_options.integrator_type == 'IRK':
            acados_integrator.set("xdot", np.zeros((nx,)))

        simX[i+1,:] = acados_integrator.simulate(x=simX[i, :], u=u0)

    S_forw = acados_integrator.get("S_forw")
    print("S_forw, sensitivities of simulaition result wrt x,u:\n", S_forw)

    # plot results
    # plot_pendulum(np.linspace(0, N*Tf, N+1), 10, np.zeros((N, nu)), simX, latexify=False)


if __name__ == "__main__":
    # full cmake build
    main(build=True, generate=True)
    # reuse code
    main(build=False, generate=False)
    # reuse code, but build with make
    main(build=True, generate=False, use_cmake=False)
    # test cython wrapper
    main(build=True, generate=True, use_cmake=False, use_cython=True)
