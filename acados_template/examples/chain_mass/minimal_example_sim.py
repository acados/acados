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

# authors: Florian Messerer, Jonathan Frey

import sys, os
import numpy as np

from acados_template import AcadosSim, AcadosSimSolver
from export_chain_mass_model import export_chain_mass_model
from utils import *
from plot_utils import *


sim = AcadosSim()

# chain parameters
chain_params = get_chain_params()

n_mass = chain_params["n_mass"]
M = chain_params["n_mass"] - 2 # number of intermediate masses
# Ts = chain_params["Ts"]
# Tsim = chain_params["Tsim"]
N = chain_params["N"]
# u_init = chain_params["u_init"]
with_wall = chain_params["with_wall"]
yPosWall = chain_params["yPosWall"]
m = chain_params["m"]
D = chain_params["D"]
L = chain_params["L"]

# export model 
model = export_chain_mass_model(n_mass, m, D, L)


#%% simulate model
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
sim.solver_options.integrator_type = 'GNSF'
sim.solver_options.newton_iter = 3 # for implicit integrator

# create
acados_integrator = AcadosSimSolver(sim)

simX = np.ndarray((N+1, nx))

# position of last mass
xPosFirstMass = np.zeros((3,1))
xEndRef = np.zeros((3,1))
xEndRef[0] = L * (M+1) * 6
# initial state
pos0_x = np.linspace(xPosFirstMass[0], xEndRef[0], n_mass)
x0 = np.zeros((nx, 1))
x0[:3*(M+1):3] = pos0_x[1:].reshape((M+1,1))
    

u0 = np.zeros((nu, 1))
acados_integrator.set("u", u0)

simX[0,:] = x0.flatten()

for i in range(N):
    # set initial state
    acados_integrator.set("x", simX[i,:])
    # solve
    status = acados_integrator.solve()
    # get solution
    simX[i+1,:] = acados_integrator.get("x")

if status != 0:
    raise Exception(f'acados returned status {status}.')

plot_chain_position_traj(simX)

#%% find steady state

xrest = compute_steady_state(n_mass, m, D, L, xPosFirstMass, xEndRef)

if os.environ.get('ACADOS_ON_CI') is None:

    plot_chain_position(xrest, xPosFirstMass)

    plot_chain_position_3D((xrest, simX[-1,:]), xPosFirstMass)

    ani = animate_chain_position_3D(simX, xPosFirstMass)
    ani = animate_chain_position(simX, xPosFirstMass)
