#
# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# Licensed under the 2-Clause BSD License.

from acados_template import AcadosSim, AcadosSimSolver
from pendulum_model import export_pendulum_ode_model
from utils import plot_pendulum
import numpy as np


def main():

    sim = AcadosSim()
    sim.model = export_pendulum_ode_model()

    Tf = 0.1
    nx = sim.model.x.rows()
    N_sim = 200

    # set simulation time
    sim.solver_options.T = Tf
    # set options
    sim.solver_options.integrator_type = 'IRK'
    sim.solver_options.num_stages = 3
    sim.solver_options.num_steps = 3
    sim.solver_options.newton_iter = 3 # for implicit integrator
    sim.solver_options.collocation_type = "GAUSS_RADAU_IIA"

    # create
    acados_integrator = AcadosSimSolver(sim)

    x0 = np.array([0.0, np.pi+1, 0.0, 0.0])
    u0 = np.array([0.0])
    xdot_init = np.zeros((nx,))

    simX = np.zeros((N_sim+1, nx))
    simX[0,:] = x0

    for i in range(N_sim):
        # Note that xdot is only used if an IRK integrator is used
        simX[i+1,:] = acados_integrator.simulate(x=simX[i,:], u=u0, xdot=xdot_init)

    S_forw = acados_integrator.get("S_forw")
    print("S_forw, sensitivities of simulation result wrt x,u:\n", S_forw)

    plot_pendulum(np.linspace(0, N_sim*Tf, N_sim+1), 10, np.repeat(u0, N_sim), simX,
                  latexify=False, time_label=sim.model.t_label, x_labels=sim.model.x_labels, u_labels=sim.model.u_labels)


if __name__ == "__main__":
    main()
