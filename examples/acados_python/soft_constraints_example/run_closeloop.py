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

from define_ocp import define_ocp
from export_pendulum_ode_model import export_pendulum_ode_model
from plot_utils import plot_ocp
from acados_template import *

# import the model
model = export_pendulum_ode_model()

# Define the ocp problem
ocp = define_ocp(model)
#  sim = define_sim(model)

#  store_ocp_solver(ocp, "test.json")
ocp = load_ocp_solver("test.json")

# generate c code
# compile the src
# wrap the solver with ctypes interface
acados_ocp_solver = generate_ocp_solver(ocp, model)
acados_sim_solver = generate_sim_solver(sim, model)

nx = ocp.dims.nx
nu = ocp.dims.nu
N = ocp.dims.N
x0 = np.array([0.0, 3.14, 0.0, 0.0])

# Close loop simulation
Nsim = 100
simX = np.ndarray((Nsim, nx))
simU = np.ndarray((Nsim, nu))

#  acados_sim_solver.set("S_forw", forward_seed)
acados_sim_solver.set("x", x0)

for i in range(Nsim):
    acados_ocp_solver.set(0, "lbx", x0)
    acados_ocp_solver.set(0, "ubx", x0)

    status = acados_ocp_solver.solve()

    # get solution
    u0 = acados_ocp_solver.get(0, "u")

    for j in range(nx):
        simX[i,j] = x0[j]

    for j in range(nu):
        simU[i,j] = u0[j]

    # simulate the systemd
    acados_sim_solver.set("u", u0)
    acados_sim_solver.set("x", x0)
    status = acados_sim_solver.solve()
    x0 = acados_sim_solver.get("xn")

    #  #  update initial condition
    #  x0 = acados_ocp_solver.get(1, "x")


    # update reference
    for j in range(N):
        acados_ocp_solver.set(j, "yref", np.array([0, 0, 0, 0, 0]))
    acados_ocp_solver.set(N, "yref", np.array([0, 0, 0, 0]))


# save solution
np.save("simX.npy", simX)
np.save("simU.npy", simU)

# plot solution
plot_ocp(dt=Tf/N, Nsim=Nsim)
