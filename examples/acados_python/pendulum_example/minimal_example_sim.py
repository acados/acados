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

from acados_template import *
from export_pendulum_ode_model import export_pendulum_ode_model
import numpy as np
import matplotlib.pyplot as plt

sim = acados_sim()

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
sim.solver_options.num_stages = 4
sim.solver_options.num_steps = 3
sim.solver_options.newton_iter = 3 # for implicit integrator


# create
acados_integrator = generate_sim_solver(sim)

simX = np.ndarray((N+1, nx))
x0 = np.array([0.0, np.pi+1, 0.0, 0.0])
u0 = np.array([0.0])
acados_integrator.set("u", u0)

simX[0,:] = x0

for i in range(N-1):
    # set initial state
    acados_integrator.set("x", simX[i,:])
    # solve
    status = acados_integrator.solve()
    # get solution
    simX[i+1,:] = acados_integrator.get("x")

if status != 0:
    raise Exception('acados returned status {}. Exiting.'.format(status))

# plot results
t = np.linspace(0.0, Tf/N, N)


# plot results
t = np.linspace(0.0, Tf/N, N)

plt.subplot(4, 1, 1)
plt.plot(t, simX[:-1,0])
plt.ylabel('p')
plt.xlabel('t')
plt.grid(True)

plt.subplot(4, 1, 2)
plt.plot(t, simX[:-1,1])
plt.ylabel('theta')
plt.xlabel('t')
plt.grid(True)

plt.subplot(4, 1, 3)
plt.plot(t, simX[:-1,2])
plt.ylabel('v')
plt.xlabel('t')
plt.grid(True)

plt.subplot(4, 1, 4)
plt.plot(t, simX[:-1,3])
plt.ylabel('dtheta')
plt.xlabel('t')
plt.grid(True)

plt.title('closed-loop simulation')
plt.subplots_adjust(left=None, bottom=None, right=None, top=None, hspace=0.4)

# avoid plotting when running on Travis
if os.environ.get('ACADOS_ON_TRAVIS') is None: 
    plt.show()

