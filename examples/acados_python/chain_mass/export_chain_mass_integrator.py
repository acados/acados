#
# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# Licensed under the 2-Clause BSD License.

import sys
import numpy as np
import scipy.linalg

from acados_template import AcadosSim, AcadosOcpSolver, AcadosSimSolver

from export_disturbed_chain_mass_model import export_disturbed_chain_mass_model

from plot_utils import *
from utils import *
import matplotlib.pyplot as plt

# create ocp object to formulate the simulation problem
sim = AcadosSim()


def export_chain_mass_integrator(n_mass, m, D, L):

    # simulation options
    Ts = 0.2

    # export model
    M = n_mass - 2 # number of intermediate masses
    model = export_disturbed_chain_mass_model(n_mass, m, D, L)

    # set model
    sim.model = model

    nx = model.x.rows()
    nu = model.u.rows()
    ny = nx + nu
    ny_e = nx

    # disturbances
    nparam = 3*M
    sim.parameter_values = np.zeros((nparam,))

    # solver options
    sim.solver_options.integrator_type = 'IRK'

    sim.solver_options.num_stages = 2
    sim.solver_options.num_steps = 2
    # sim.solver_options.nlp_solver_tol_eq = 1e-9

    # set prediction horizon
    sim.solver_options.T = Ts

    # acados_ocp_solver = AcadosOcpSolver(ocp)
    acados_integrator = AcadosSimSolver(sim)

    return acados_integrator