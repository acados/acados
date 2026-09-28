#
# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# Licensed under the 2-Clause BSD License.

# authors: Katrin Baumgaertner, Jonathan Frey

from acados_template import AcadosSim, AcadosSimSolver
import numpy as np


def setup_acados_integrator(
    model,
    dt,
    cstr_param,
    sensitivity_propagation=False,
    num_stages=4,
    num_steps=100,
    integrator_type="ERK",
):

    sim = AcadosSim()

    # set model
    sim.model = model

    # set simulation time
    sim.solver_options.T = dt

    ## set options
    sim.solver_options.integrator_type = integrator_type
    sim.solver_options.num_stages = num_stages
    sim.solver_options.num_steps = num_steps
    # for implicit integrator
    sim.solver_options.newton_iter = 10
    sim.solver_options.newton_tol = 1e-8
    sim.solver_options.collocation_type = "GAUSS_LEGENDRE"
    # sensitivity_propagation
    sim.solver_options.sens_adj = sensitivity_propagation
    sim.solver_options.sens_forw = sensitivity_propagation
    sim.solver_options.sens_hess = sensitivity_propagation

    # nominal parameter values
    sim.parameter_values = np.array([cstr_param.F0])

    # create
    acados_integrator = AcadosSimSolver(sim)

    return acados_integrator


if __name__ == "__main__":
    setup_acados_integrator()
