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
import acados_template as at
from export_ode_model import *
import numpy as np
import scipy.linalg
from ctypes import *

model_name = "example"
sim_shared_lib = "libacados_sim_{}".format(model_name)

acados_solver = generate_solver(sim_shared_lib)
acados_simulator = acados_sim(shared_lib)

Nsim = 100

simX = np.ndarray((Nsim, nx))
simU = np.ndarray((Nsim, nu))

for i in range(Nsim):
    status = acados_solver.solve()

    # get solution
    x0 = acados_solver.get(0, "x")
    u0 = acados_solver.get(0, "u")
    
    for j in range(nx):
        simX[i,j] = x0[j]

    for j in range(nu):
        simU[i,j] = u0[j]
    
    # update initial condition
    x0 = acados_solver.get(1, "x")

    acados_solver.set(0, "lbx", x0)
    acados_solver.set(0, "ubx", x0)

    # update reference
    for j in range(N):
        acados_solver.set(j, "yref", np.array([0, 0, 0, 0, 0]))
    acados_solver.set(N, "yref", np.array([0, 0, 0, 0]))

# plot results
import matplotlib
import matplotlib.pyplot as plt
t = np.linspace(0.0, Tf/N, Nsim)
plt.subplot(2, 1, 1)
plt.step(t, simU, color='r')
plt.title('closed-loop simulation')
plt.ylabel('u')
plt.xlabel('t')
plt.grid(True)
plt.subplot(2, 1, 2)
plt.plot(t, simX[:,1])
plt.ylabel('theta')
plt.xlabel('t')
plt.grid(True)
plt.show()

