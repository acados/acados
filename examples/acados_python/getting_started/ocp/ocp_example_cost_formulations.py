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

import sys
sys.path.insert(0, '../common')

from acados_template import AcadosOcp, AcadosOcpSolver
from export_pendulum_ode_model import export_pendulum_ode_model
import numpy as np
import scipy.linalg
from utils import plot_pendulum
from casadi import vertcat

COST_MODULE = 'EXTERNAL' # 'LS', 'EXTERNAL'
HESSIAN_APPROXIMATION = 'EXACT' # 'GAUSS_NEWTON
EXTERNAL_COST_USE_NUM_HESS = 1

# create ocp object to formulate the OCP
ocp = AcadosOcp()

# set model
model = export_pendulum_ode_model()
ocp.model = model

Tf = 1.0
nx = model.x.size()[0]
nu = model.u.size()[0]
ny = nx + nu
ny_e = nx
N = 20

# set dimensions
ocp.dims.nx = nx
ocp.dims.ny = ny
ocp.dims.ny_e = ny_e
ocp.dims.nbu = nu 
ocp.dims.nu = nu
ocp.dims.N = N

# set cost
Q = 2*np.diag([1e3, 1e3, 1e-2, 1e-2])
R = 2*np.diag([1e-2])

x = ocp.model.x
u = ocp.model.u

ocp.cost.W_e = Q
ocp.cost.W = scipy.linalg.block_diag(Q, R)

if COST_MODULE == 'LS':
    ocp.cost.cost_type = 'LINEAR_LS'
    ocp.cost.cost_type_e = 'LINEAR_LS'

    ocp.cost.Vx = np.zeros((ny, nx))
    ocp.cost.Vx[:nx,:nx] = np.eye(nx)

    Vu = np.zeros((ny, nu))
    Vu[4,0] = 1.0
    ocp.cost.Vu = Vu

    ocp.cost.Vx_e = np.eye(nx)

elif COST_MODULE == 'NLS':
    ocp.cost.cost_type = 'NONLINEAR_LS'
    ocp.cost.cost_type_e = 'NONLINEAR_LS'

    ocp.model.cost_y_expr = vertcat(x, u)
    ocp.model.cost_y_expr_e = x

elif COST_MODULE == 'EXTERNAL':
    ocp.cost.cost_type = 'EXTERNAL'
    ocp.cost.cost_type_e = 'EXTERNAL'

    ocp.model.cost_expr_ext_cost = vertcat(x, u).T @ ocp.cost.W @ vertcat(x, u)
    ocp.model.cost_expr_ext_cost_e = x.T @ Q @ x

else:
    raise Exception('Unknown COST_MODULE. Possible values are \'LS\' and \'NLS\'.')

ocp.cost.yref = np.zeros((ny, ))
ocp.cost.yref_e = np.zeros((ny_e, ))

# set constraints
Fmax = 80
ocp.constraints.constr_type = 'BGH'
ocp.constraints.lbu = np.array([-Fmax])
ocp.constraints.ubu = np.array([+Fmax])
ocp.constraints.x0 = np.array([0.0, np.pi, 0.0, 0.0])
ocp.constraints.idxbu = np.array([0])

ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM' # FULL_CONDENSING_QPOASES
ocp.solver_options.hessian_approx = HESSIAN_APPROXIMATION
ocp.solver_options.regularize_method = 'CONVEXIFY'
ocp.solver_options.integrator_type = 'ERK'

ocp.solver_options.qp_solver_cond_N = 5

# set prediction horizon
ocp.solver_options.tf = Tf
ocp.solver_options.nlp_solver_type = 'SQP' # SQP_RTI
ocp.solver_options.ext_cost_num_hess = EXTERNAL_COST_USE_NUM_HESS

ocp_solver = AcadosOcpSolver(ocp, json_file = 'acados_ocp.json')

# from casadi import jacobian
# ux = vertcat(ocp.model.u, ocp.model.x)
# jacobian(jacobian(ocp.model.cost_expr_ext_cost, ux), ux)
# SX(@1=0.04, @2=4000,
# [[@1, 00, 00, 00, 00],
#  [00, @2, 00, 00, 00],
#  [00, 00, @2, 00, 00],
#  [00, 00, 00, @1, 00],
#  [00, 00, 00, 00, @1]])

# NOTE: hessian is wrt [u,x]
if EXTERNAL_COST_USE_NUM_HESS:
    for i in range(N):
        ocp_solver.cost_set(i, "ext_cost_num_hess", np.diag([0.04, 4000, 4000, 0.04, 0.04, ]))
    ocp_solver.cost_set(N, "ext_cost_num_hess", np.diag([4000, 4000, 0.04, 0.04, ]))


simX = np.ndarray((N+1, nx))
simU = np.ndarray((N, nu))

status = ocp_solver.solve()

ocp_solver.print_statistics()

if status != 0:
    raise Exception('acados returned status {}. Exiting.'.format(status))

# get solution
for i in range(N):
    simX[i,:] = ocp_solver.get(i, "x")
    simU[i,:] = ocp_solver.get(i, "u")
simX[N,:] = ocp_solver.get(N, "x")


plot_pendulum(np.linspace(0, Tf, N+1), Fmax, simU, simX, latexify=False)

