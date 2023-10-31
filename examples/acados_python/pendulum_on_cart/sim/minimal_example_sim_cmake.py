# -*- coding: future_fstrings -*-
#
# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# The 2-Clause BSD License
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice,
# this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation
# and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.;
#

import sys
sys.path.insert(0, '../common')

from acados_template import AcadosSim, AcadosSimSolver, sim_get_default_cmake_builder
from pendulum_model import export_pendulum_ode_model
# from utils import plot_pendulum
import numpy as np


def main(build=True, generate=True, use_cmake=True, use_cython=False):
    sim = AcadosSim()

    model = export_pendulum_ode_model()
    sim.model = model

    Tf = 0.1
    nx = model.x.size()[0]
    nu = model.u.size()[0]
    N = 200

    # set simulation time
    sim.solver_options.T = Tf
    # set options
    sim.solver_options.integrator_type = 'IRK'
    sim.solver_options.num_stages = 3
    sim.solver_options.num_steps = 3
    sim.solver_options.newton_iter = 3 # for implicit integrator
    sim.solver_options.collocation_type = "GAUSS_RADAU_IIA"

    if use_cmake:
        # use the CMake build pipeline
        cmake_builder = sim_get_default_cmake_builder()
    else:
        cmake_builder = None

    # create
    # acados_integrator = AcadosSimSolver(sim, cmake_builder=cmake_builder)
    if use_cython:
        AcadosSimSolver.generate(sim, json_file='acados_sim.json')
        AcadosSimSolver.build(sim.code_export_directory, with_cython=True)
        acados_integrator = AcadosSimSolver.create_cython_solver('acados_sim.json')
    else:
        acados_integrator = AcadosSimSolver(sim, cmake_builder=cmake_builder, generate=generate, build=build)

    simX = np.ndarray((N+1, nx))
    x0 = np.array([0.0, np.pi+1, 0.0, 0.0])
    u0 = np.array([0.0])
    acados_integrator.set("u", u0)

    simX[0,:] = x0

    for i in range(N):
        # set initial state
        acados_integrator.set("x", simX[i,:])
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
