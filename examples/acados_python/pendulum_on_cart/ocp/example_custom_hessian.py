# -*- coding: future_fstrings -*-
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

import sys
sys.path.insert(0, '../common')

from acados_template import AcadosOcp, AcadosOcpSolver
from pendulum_model import export_pendulum_ode_model
import numpy as np
from utils import plot_pendulum
from casadi import tanh, SX, Sparsity, hessian, if_else, horzcat, DM, blockcat

# create ocp object to formulate the OCP
ocp = AcadosOcp()

# set model
model = export_pendulum_ode_model()
ocp.model = model

Tf = 2.
nx = model.x.rows()
nu = model.u.rows()
ny = nx + nu
ny_e = nx
N = 40

# set dimensions
ocp.solver_options.N_horizon = N

# set cost
ocp.cost.cost_type = 'EXTERNAL'
ocp.cost.cost_type_e = 'EXTERNAL'
W_u = 1e-3
theta = model.x[1]

tanh_theta_squared = tanh(theta)**2

ocp.model.cost_expr_ext_cost = tanh_theta_squared + .5*(model.x[0]**2 + W_u*model.u**2)
ocp.model.cost_expr_ext_cost_e = tanh_theta_squared + .5*model.x[0]**2

custom_hess_u = W_u

J = horzcat(SX.eye(2), SX(2,2))

print(DM(J.sparsity()))

# diagonal matrix with second order terms of outer loss function.
D = SX.sym('D', Sparsity.diag(2))
D[0, 0] = 1
[hess_tan, grad_tan] = hessian(tanh_theta_squared, theta)
D[1, 1] = if_else(theta == 0, hess_tan, grad_tan/theta)

custom_hess_x = J.T @ D @ J

zeros = SX(1, nx)
cost_expr_ext_cost_custom_hess = blockcat(custom_hess_u, zeros, zeros.T, custom_hess_x)
cost_expr_ext_cost_custom_hess_e = custom_hess_x

ocp.model.cost_expr_ext_cost_custom_hess = cost_expr_ext_cost_custom_hess
ocp.model.cost_expr_ext_cost_custom_hess_e = cost_expr_ext_cost_custom_hess_e

# set constraints
Fmax = 35
ocp.constraints.lbu = np.array([-Fmax])
ocp.constraints.ubu = np.array([+Fmax])
ocp.constraints.idxbu = np.array([0])

x0 = np.array([0.0, np.pi, 0.0, 0.0])
xf = np.array([0.0, 0.0, 0.0, 0.0])
ocp.constraints.x0 = x0

# set options
ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM' # FULL_CONDENSING_QPOASES
ocp.solver_options.integrator_type = 'ERK'
ocp.solver_options.nlp_solver_type = 'SQP' # SQP_RTI, SQP
ocp.solver_options.globalization = 'MERIT_BACKTRACKING'
ocp.solver_options.nlp_solver_max_iter = 150
ocp.solver_options.store_iterates = True # store intermediate SQP iterates to access them after the solve
ocp.solver_options.tol = 1e-4

# set prediction horizon
ocp.solver_options.tf = Tf

ocp_solver = AcadosOcpSolver(ocp, json_file = 'acados_ocp.json')
ocp_solver.options_set("globalization_line_search_use_sufficient_descent", 1)
ocp_solver.options_set("globalization_full_step_dual", 1)

simX = np.zeros((N+1, nx))
simU = np.zeros((N, nu))

# initialization
for i, tau in enumerate(np.linspace(0, 1, N+1)):
    ocp_solver.set(i, 'x', x0*(1-tau) + tau*xf)

status = ocp_solver.solve()

if status != 0:
    print(f'acados returned status {status}.')

# get solution
for i in range(N):
    simX[i,:] = ocp_solver.get(i, "x")
    simU[i,:] = ocp_solver.get(i, "u")
simX[N,:] = ocp_solver.get(N, "x")

ocp_solver.print_statistics()

iteration = 2
print(f"primal iterate at iteration {iteration}:")
iterate = ocp_solver.get_iterate(iteration)
print(iterate.x_traj)
# get all iterates
iterates = ocp_solver.get_iterates()
# print(iterates.iterate_traj[-1])

plot_pendulum(np.linspace(0, Tf, N+1), Fmax, simU, simX, latexify=False)
