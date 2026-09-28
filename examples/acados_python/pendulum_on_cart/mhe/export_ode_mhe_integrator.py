#
# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# Licensed under the 2-Clause BSD License.

from acados_template import *
import numpy as np


def export_ode_mhe_integrator(model, h, use_cython=True):

    sim = AcadosSim()

    # set model
    sim.model = model
    sim.parameter_values = np.zeros(np.prod(model.p.shape))

    # set simulation time
    sim.solver_options.T = h
    sim.solver_options.num_stages = 4
    sim.solver_options.num_steps = 3
    sim.solver_options.newton_iter = 3 # for implicit integrator

    if use_cython:
        acados_integrator = AcadosSimSolver.create_cython_solver(sim)
    else:
        acados_integrator = AcadosSimSolver(sim)

    return acados_integrator
