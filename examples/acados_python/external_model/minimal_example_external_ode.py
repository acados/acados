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

# NOTE: requires CasADi system installation and compilation
#    of the test external library. To do so:
# cd test_external_lib
# mkdir build
# cd build
# cmake ..
# make
#

from acados_template import AcadosOcp, AcadosOcpSolver, AcadosSimSolver
from export_external_ode_model import export_external_ode_model
import numpy as np
import os
import scipy.linalg
import matplotlib.pyplot as plt

# os.environ["ACADOS_SOURCE_DIR"] =

# create ocp object to formulate the OCP
ocp = AcadosOcp()

# set model
model = export_external_ode_model()
ocp.model = model
ocp.solver_options.model_external_shared_lib_dir = os.getcwd()+"/test_external_lib/build"
ocp.solver_options.model_external_shared_lib_name = "external_ode_casadi"

Tf = 1.0
nx = model.x.size()[0]
nu = model.u.size()[0]
ny = nx + nu
ny_e = nx
N = 30
x0 = np.array([0, 0])
xT = np.array([1/2, 1])

# set dimensions
ocp.dims.N = N

# set cost module
ocp.cost.cost_type = 'LINEAR_LS'
ocp.cost.cost_type_e = 'LINEAR_LS'

Q = 2*np.diag([1e-3, 1e-3])
R = 2*np.diag([1e-4])
ocp.cost.W = scipy.linalg.block_diag(Q, R)
ocp.cost.W_e = Q

ocp.cost.Vx = np.zeros((ny, nx))

Vu = np.zeros((ny, nu))
Vu[2,0] = 1.0
ocp.cost.Vu = Vu

ocp.cost.Vx_e = np.zeros((nx,nx))
ocp.cost.yref  = np.zeros((ny, ))
ocp.cost.yref_e = np.zeros((ny_e, ))

# set constraints
Fmax = 10
ocp.constraints.x0 = x0
ocp.dims.nbx_0 = nx
# ocp.constraints.lbu = np.array([-Fmax])
# ocp.constraints.ubu = np.array([+Fmax])
# ocp.constraints.idxbu = np.array([0])

# terminal constraints
ocp.constraints.Jbx_e = np.eye(nx)
ocp.constraints.ubx_e = xT
ocp.constraints.lbx_e = xT
ocp.constraints.idxbx_e = np.array(range(nx))
ocp.dims.nbx_e = nx

ocp.solver_options.tf = Tf
ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM' # FULL_CONDENSING_QPOASES
ocp.solver_options.hessian_approx = 'GAUSS_NEWTON'
ocp.solver_options.integrator_type = 'ERK'
ocp.solver_options.nlp_solver_type = 'SQP' # SQP_RTI
ocp_solver = AcadosOcpSolver(ocp, json_file = 'acados_ocp.json')

# initial guess
t_traj = np.linspace(0, Tf, N+1)
x_traj = np.linspace(x0,xT,N+1)
u_traj = np.ones((N,1))+np.random.rand(N,1)*1e-6
for n in range(N+1):
    ocp_solver.set(n, 'x', x_traj[n,:])
for n in range(N):
    ocp_solver.set(n, 'u', u_traj[n])


# solve
status = ocp_solver.solve()

if status != 0:
    raise Exception(f'acados returned status {status}.')

# get solution
stat_fields = ['time_tot', 'time_lin', 'time_qp', 'time_qp_solver_call', 'time_reg', 'sqp_iter']
for field in stat_fields:
    print(f"{field} : {ocp_solver.get_stats(field)}")
simX = np.ndarray((N + 1, nx))
simU = np.ndarray((N, nu))
for i in range(N):
    simX[i,:] = ocp_solver.get(i, "x")
    simU[i,:] = ocp_solver.get(i, "u")
simX[N,:] = ocp_solver.get(N, "x")

print(simX)
plt.plot(simX[:,:nu],'o',label='opt_sol')
plt.plot(x_traj[:,0],'x',label='init_sol')
plt.legend()
plt.title('position')
plt.figure()
plt.plot(simX[:,nu:],'o',label='opt_sol')
plt.plot(x_traj[:,1],'x',label='init_sol')
plt.legend()
plt.title('velocity')
plt.figure()
plt.plot(simU,'o',label='opt_sol')
plt.plot(u_traj,'x',label='init_traj')
plt.legend()
plt.title('control')
plt.show(block=True)