import numpy as np
from acados_template import AcadosSim, AcadosSimSolver, AcadosSimRosOptions

import sys
import os
script_dir = os.path.dirname(os.path.realpath(__file__))
common_path = os.path.join(script_dir, '..', 'common')
sys.path.insert(0, os.path.abspath(common_path))
from pendulum_model import export_pendulum_ode_model


def main():
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

    sim.ros_opts = AcadosSimRosOptions()
    sim.ros_opts.package_name = "pendulum_on_cart_sim"

    export_code = os.path.join(script_dir, 'generated_sim')
    sim.code_export_directory = str( os.path.join(export_code, "c_generated_code"))
    acados_integrator = AcadosSimSolver(sim, json_file=str(os.path.join(export_code, 'acados_sim.json')))

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


if __name__ == "__main__":
    main()
