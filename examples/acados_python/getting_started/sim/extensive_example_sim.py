#
# Copyright 2019 Gianluca Frison, Dimitris Kouzoupis, Robin Verschueren,
# Andrea Zanelli, Niels van Duijkeren, Jonathan Frey, Tommaso Sartor,
# Branimir Novoselnik, Rien Quirynen, Rezart Qelibari, Dang Doan,
# Jonas Koenemann, Yutao Chen, Tobias Schöls, Jonas Schlagenhauf, Moritz Diehl
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

import os, sys, json
sys.path.insert(0, '../common')

from acados_template import AcadosSim, AcadosSimSolver, acados_dae_model_json_dump
from export_pendulum_ode_model import export_pendulum_ode_model
from utils import plot_pendulum
import numpy as np
import matplotlib.pyplot as plt

sim = AcadosSim()

# export model 
model = export_pendulum_ode_model()

# set model_name 
sim.model = model

Tf = 0.1
nx = model.x.size()[0]
nu = model.u.size()[0]
N = 200

# set simulation time
sim.solver_options.T = Tf
# set options
sim.solver_options.num_stages = 4
sim.solver_options.num_steps = 3
sim.solver_options.newton_iter = 3 # for implicit integrator
sim.solver_options.integrator_type = "IRK" # ERK, IRK, GNSF
sim.solver_options.sens_forw = True
sim.solver_options.sens_adj = True
sim.solver_options.sens_hess = True
sim.solver_options.sens_algebraic = False
sim.solver_options.output_z = False
sim.solver_options.jac_reuse = False


if sim.solver_options.integrator_type == "GNSF":
    # Perform GNSF structure detection in Octave
    # export OCTAVE_PATH=$OCTAVE_PATH:$ACADOS_INSTALL_DIR/external/casadi-octave
    # export OCTAVE_PATH=$OCTAVE_PATH:$ACADOS_INSTALL_DIR/interfaces/acados_matlab_octave/
    # export OCTAVE_PATH=$OCTAVE_PATH:$ACADOS_INSTALL_DIR/interfaces/acados_matlab_octave/acados_template_mex/

    # acados_dae_model_json_dump(model)
    # status = os.system('octave convert_dae2gnsf.m')
    with open(model.name + '_gnsf_functions.json', 'r') as f:
        gnsf_dict = json.load(f)
    sim.gnsf_model = gnsf_dict

# create
acados_integrator = AcadosSimSolver(sim)

simX = np.ndarray((N+1, nx))
x0 = np.array([0.0, np.pi+1, 0.0, 0.0])
u0 = np.array([0.0])
acados_integrator.set("u", u0)

simX[0,:] = x0

acados_integrator.set("S_adj", np.ones((nx+nu, 1)))

for i in range(N):
    # set initial state
    acados_integrator.set("x", simX[i,:])
    # solve
    status = acados_integrator.solve()
    # get solution
    simX[i+1,:] = acados_integrator.get("x")

if status != 0:
    raise Exception('acados returned status {}. Exiting.'.format(status))

S_forw = acados_integrator.get("S_forw")
S_hess = acados_integrator.get("S_hess")
S_adj = acados_integrator.get("S_adj")

print("S_forw, sensitivities of simulaition result wrt x,u:\n", S_forw)
print("S_adj, adjoint sensitivities:\n", S_adj)
print("S_hess, second order sensitivities:\n", S_hess)

# plot results
plot_pendulum(np.linspace(0, Tf, N+1), 10, np.zeros((N, nu)), simX, latexify=False)
