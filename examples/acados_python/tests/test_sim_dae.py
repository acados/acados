#
# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# Licensed under the 2-Clause BSD License.

import sys
sys.path.insert(0, '../pendulum_on_cart/common')

from acados_template import AcadosSim, AcadosSimSolver
from pendulum_model import export_augmented_pendulum_model
from utils import plot_pendulum
import numpy as np

sim = AcadosSim()

# export model
model = export_augmented_pendulum_model()

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


# create
acados_integrator = AcadosSimSolver(sim)

simX = np.zeros((N+1, nx))
x0 = np.array([0.0, np.pi+1, 0.0, 0.0])

u0_val = 2.0
u0 = np.array([u0_val])

# test setter
acados_integrator.set("u", 2)
acados_integrator.set("u", 2.0)
acados_integrator.set("u", u0)

simX[0,:] = x0

for i in range(N):
    # set initial state
    acados_integrator.set("x", simX[i,:])
    # solve
    status = acados_integrator.solve()
    # get solution
    simX[i+1,:] = acados_integrator.get("x")

if status != 0:
    raise Exception(f'acados returned status {status}.')

S_algebraic = acados_integrator.get("S_algebraic")
print("S_algebraic (dz_dxu) = ", S_algebraic)

z = acados_integrator.get("z")
print("z = ", z)
last_x0 = simX[N-1,:]
print(f"{last_x0 = }")

z_analytic = np.array([last_x0[0], u0_val**2])
err_z = np.abs(z - z_analytic)
if np.any(err_z > 1e-6):
    raise Exception(f'z and z_analytic should match! Difference is {err_z}')
print("Success: z and z_analytic match!")

S_forw = acados_integrator.get("S_forw")
print("S_forw, sensitivities of simulaition result wrt x,u:\n", S_forw)

# plot results
plot_pendulum(np.linspace(0, Tf, N+1), 10, u0_val * np.ones((N, nu)), simX, latexify=False)
