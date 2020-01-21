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
from export_pendulum_ode_model import *
import numpy as np
import scipy.linalg
from ctypes import *
import matplotlib.pyplot as plt

# FORMULATION = 'NLS' # 'LS'
FORMULATION = 'LS' # 'LS'

# create ocp object to formulate the OCP
ocp = acados_ocp_nlp()

# export model 
model = export_pendulum_ode_model()

# set model_name 
ocp.model = model

Tf = 1.0
nx = model.x.size()[0]
nu = model.u.size()[0]
ny = nx + nu
ny_e = nx
N = 20

# set ocp_nlp_dimensions
nlp_dims     = ocp.dims
nlp_dims.nx  = nx 
nlp_dims.ny  = ny 
nlp_dims.ny_e = ny_e 
nlp_dims.nbx = 0
nlp_dims.nbu = nu 
nlp_dims.nu  = nu
nlp_dims.N   = N

# set cost module
nlp_cost = ocp.cost

if FORMULATION == 'LS':
    nlp_cost.cost_type = 'LINEAR_LS'
    nlp_cost.cost_type_e = 'LINEAR_LS'
elif FORMULATION == 'NLS':
    nlp_cost.cost_type = 'NONLINEAR_LS'
    nlp_cost.cost_type_e = 'NONLINEAR_LS'
else:
    raise Exception('Unknown FORMULATION. Possible values are \'LS\' and \'NLS\'.')

Q = 2*np.diag([1e3, 1e3, 1e-2, 1e-2])
R = 2*np.diag([1e-2])

if FORMULATION == 'NLS':
    nlp_cost.W = scipy.linalg.block_diag(R, Q) 
else:
    nlp_cost.W = scipy.linalg.block_diag(Q, R) 

nlp_cost.W_e = Q 

if FORMULATION == 'LS':
    nlp_cost.Vx = np.zeros((ny, nx))
    nlp_cost.Vx[:nx,:nx] = np.eye(nx)

    Vu = np.zeros((ny, nu))
    Vu[4,0] = 1.0
    nlp_cost.Vu = Vu

    nlp_cost.Vx_e = np.eye(nx)

elif FORMULATION == 'NLS':
    x = SX.sym('x', 4, 1)
    u = SX.sym('u', 1, 1)
    ocp.cost_r.expr = vertcat(u, x) 
    ocp.cost_r.x = x
    ocp.cost_r.u = u 
    ocp.cost_r.name = 'lin_res' 
    ocp.cost_r.ny = nx + nu

    ocp.cost_r_e.expr = x
    ocp.cost_r_e.x = x 
    ocp.cost_r_e.name = 'lin_res' 
    ocp.cost_r_e.ny = nx 
else:
    raise Exception("Invalid cost formulation. Exiting.")


nlp_cost.yref  = np.zeros((ny, ))
nlp_cost.yref_e = np.zeros((ny_e, ))

# set constraints
Fmax = 80
ocp.constraints.constr_type = 'BGH'
ocp.constraints.lbu = np.array([-Fmax])
ocp.constraints.ubu = np.array([+Fmax])
ocp.constraints.x0 = np.array([0.0, np.pi, 0.0, 0.0])
# ocp.constraints.x0 = np.array([0.0, 0.5, 0.0, 0.0])
ocp.constraints.idxbu = np.array([0])

ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM'
# ocp.solver_options.qp_solver = 'FULL_CONDENSING_QPOASES'
ocp.solver_options.hessian_approx = 'GAUSS_NEWTON'
ocp.solver_options.integrator_type = 'ERK'

ocp.solver_options.qp_solver_cond_N = 5

# set prediction horizon
ocp.solver_options.tf = Tf
ocp.solver_options.nlp_solver_type = 'SQP'
# ocp.solver_options.nlp_solver_type = 'SQP_RTI'

# set header path
ocp.acados_include_path  = '../../../../include'
ocp.acados_lib_path      = '../../../../lib'

acados_solver = generate_solver(ocp, json_file = 'acados_ocp.json')


simX = np.ndarray((N+1, nx))
simU = np.ndarray((N, nu))

status = acados_solver.solve()

if status != 0:
    raise Exception('acados returned status {}. Exiting.'.format(status))

# get solution
for i in range(N):
    simX[i,:] = acados_solver.get(i, "x")
    simU[i,:] = acados_solver.get(i, "u")
simX[N,:] = acados_solver.get(N, "x")


# plot results
t = np.linspace(0.0, Tf/N, N)

plt.subplot(5, 1, 1)
plt.step(t, simU, color='r')
plt.ylabel('u')
plt.xlabel('t')
plt.ylim([-Fmax, Fmax])
plt.grid(True)

plt.subplot(5, 1, 2)
plt.plot(t, simX[:-1,0])
plt.ylabel('p')
plt.xlabel('t')
plt.grid(True)

plt.subplot(5, 1, 3)
plt.plot(t, simX[:-1,1])
plt.ylabel('theta')
plt.xlabel('t')
plt.grid(True)

plt.subplot(5, 1, 4)
plt.plot(t, simX[:-1,2])
plt.ylabel('v')
plt.xlabel('t')
plt.grid(True)

plt.subplot(5, 1, 5)
plt.plot(t, simX[:-1,3])
plt.ylabel('dtheta')
plt.xlabel('t')
plt.grid(True)

plt.title('closed-loop simulation')
plt.subplots_adjust(left=None, bottom=None, right=None, top=None, hspace=0.4)

# avoid plotting when running on Travis
if os.environ.get('ACADOS_ON_TRAVIS') is None: 
    plt.show()

