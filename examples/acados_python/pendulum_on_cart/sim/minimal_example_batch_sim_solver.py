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

from acados_template import AcadosSim, AcadosSimSolver, AcadosSimBatchSolver
from pendulum_model import export_pendulum_ode_model
import numpy as np
import time


"""
This example shows how the AcadosSimBatchSolver can be used to parallelize mulitple integrators.

If you want to use the batch solver, make sure to compile acados with openmp and num_threads set to 1,
i.e. with the flags -DACADOS_WITH_OPENMP=ON -DACADOS_NUM_THREADS=1
The number of threads for the batch solver is then set via the option `num_threads_in_batch_solve`, see below.
"""


def setup_integrator(num_threads_in_batch_solve=1):

    sim = AcadosSim()
    sim.model = export_pendulum_ode_model()
    sim.solver_options.T = 0.2
    sim.solver_options.integrator_type = 'IRK'
    sim.solver_options.num_stages = 5
    sim.solver_options.num_steps = 10
    sim.solver_options.newton_iter = 10 # for implicit integrator
    sim.solver_options.collocation_type = "GAUSS_RADAU_IIA"
    sim.solver_options.num_threads_in_batch_solve = num_threads_in_batch_solve

    return sim


def main_sequential(x0, u0, N_sim):

    nx = x0.shape[0]
    simX = np.zeros((N_sim+1, nx))

    sim = setup_integrator()
    integrator = AcadosSimSolver(sim, verbose=False)

    simX[0,:] = x0

    t0 = time.time()
    for i in range(N_sim):
        simX[i+1,:] = integrator.simulate(x=simX[i, :], u=u0, xdot=np.zeros((nx,)))

    t_elapsed = 1e3 * (time.time() - t0)
    print("main_sequential:", f"{t_elapsed:.3f}ms")

    return simX


def main_batch(Xinit, u0, num_threads_in_batch_solve=1):

    N_batch = Xinit.shape[0] - 1
    sim = setup_integrator(num_threads_in_batch_solve)
    batch_integrator = AcadosSimBatchSolver(sim, N_batch, verbose=False)

    for n in range(N_batch):
        batch_integrator.sim_solvers[n].set("u", u0)
        batch_integrator.sim_solvers[n].set("x", Xinit[n])

    t0 = time.time()
    batch_integrator.solve()
    t_elapsed = 1e3 * (time.time() - t0)

    print(f"main_batch: with {num_threads_in_batch_solve} threads, timing: {t_elapsed:.3f}ms")

    for n in range(N_batch):
        x = batch_integrator.sim_solvers[n].get("x")
        assert np.linalg.norm(x-Xinit[n+1]) < 1e-10


if __name__ == "__main__":

    N_batch = 256
    x0 = np.array([0.0, np.pi+1, 0.0, 0.0])
    u0 = np.array([0.0])

    simX = main_sequential(x0=x0, u0=u0, N_sim=N_batch)

    main_batch(Xinit=simX, u0=u0, num_threads_in_batch_solve=1)
    main_batch(Xinit=simX, u0=u0, num_threads_in_batch_solve=4)

