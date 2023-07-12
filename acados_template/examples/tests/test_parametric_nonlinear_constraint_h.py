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
sys.path.insert(0, '../pendulum_on_cart/common')

from acados_template import AcadosOcp, AcadosOcpSolver
from pendulum_model import export_pendulum_ode_model
import numpy as np
import scipy.linalg
from utils import plot_pendulum
from casadi import SX, vertcat


COST_MODULE = 'EXTERNAL'

# create ocp object to formulate the OCP
ocp = AcadosOcp()

# set model
model = export_pendulum_ode_model()
ocp.model = model
x = model.x
u = model.u

Tf = 1.0
nx = x.size()[0]
nu = u.size()[0]
ny = nx + nu
ny_e = nx
N = 20

# set dimensions
ocp.dims.N = N

# set cost
Q = 2*np.diag([1e3, 1e3, 1e-2, 1e-2])
R = 2*np.diag([1e-2])
cost_W = scipy.linalg.block_diag(Q, R)


#
n_param = 42
p = SX.sym('p', n_param)
constraint_quotient = p[0]
y_param = p[1:nx+nu+1]
ocp.model.p = p

if COST_MODULE == 'LS':
    ocp.cost.W_e = Q
    ocp.cost.W = cost_W

    ocp.cost.cost_type = 'LINEAR_LS'
    ocp.cost.cost_type_e = 'LINEAR_LS'

    ocp.cost.Vx = np.zeros((ny, nx))
    ocp.cost.Vx[:nx,:nx] = np.eye(nx)

    Vu = np.zeros((ny, nu))
    Vu[4,0] = 1.0
    ocp.cost.Vu = Vu

    ocp.cost.Vx_e = np.eye(nx)

    ocp.cost.yref = np.zeros((ny, ))
    ocp.cost.yref_e = np.zeros((ny_e, ))

# elif COST_MODULE == 'NLS':
#     ocp.cost.cost_type = 'NONLINEAR_LS'
#     ocp.cost.cost_type_e = 'NONLINEAR_LS'

#     ocp.model.cost_y_expr = vertcat(x, u)
#     ocp.model.cost_y_expr_e = x

elif COST_MODULE == 'EXTERNAL':
    ocp.cost.cost_type = 'EXTERNAL'
    ocp.cost.cost_type_e = 'EXTERNAL'

    residual = y_param - vertcat(x, u)
    ocp.model.cost_expr_ext_cost = residual.T @ cost_W @ residual
    res_e = y_param[0:nx] - x
    ocp.model.cost_expr_ext_cost_e = res_e.T @ Q @ res_e

# set constraints

Fmax = 80
# use equivalent formulation with h constraint
# ocp.constraints.lbu = np.array([-Fmax])
# ocp.constraints.ubu = np.array([+Fmax])
# ocp.constraints.idxbu = np.array([0])

ocp.constraints.lh = np.array([-Fmax])
ocp.constraints.uh = np.array([+Fmax])
ocp.model.con_h_expr = model.u / constraint_quotient

p_0 = np.zeros(n_param)
p_0[0] = 1.0
p_0[1] = 1.0
ocp.parameter_values = p_0

ocp.constraints.x0 = np.array([0.0, np.pi, 0.0, 0.0])

ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM' # FULL_CONDENSING_QPOASES
ocp.solver_options.hessian_approx = 'EXACT' # GAUSS_NEWTON, EXACT
ocp.solver_options.regularize_method = 'CONVEXIFY' # GAUSS_NEWTON, EXACT
ocp.solver_options.integrator_type = 'ERK'
ocp.solver_options.print_level = 0
ocp.solver_options.nlp_solver_type = 'SQP' # SQP_RTI, SQP

# set prediction horizon
ocp.solver_options.tf = Tf

# Cython
if 0:
    AcadosOcpSolver.generate(ocp, json_file='acados_ocp.json')
    AcadosOcpSolver.build(ocp.code_export_directory, with_cython=True)
    ocp_solver = AcadosOcpSolver.create_cython_solver('acados_ocp.json')
else:
    ocp_solver = AcadosOcpSolver(ocp, json_file = 'acados_ocp.json')

for i in range(N):
    ## Two equivalent ways to set parameters
    if i < N/2:
        # set all parameters
        ocp_solver.set(i, "p", p_0)
    else:
        # set subset of parameters
        ocp_solver.set_params_sparse(i, np.array(range(n_param)), np.zeros(n_param))
        ocp_solver.set_params_sparse(i, np.ascontiguousarray([0, 1]), np.array([1.0, 1.0]))

simX = np.ndarray((N+1, nx))
simU = np.ndarray((N, nu))

status = ocp_solver.solve()

if status != 0:
    raise Exception(f'acados returned status {status}.')

# get solution
for i in range(N):
    simX[i,:] = ocp_solver.get(i, "x")
    simU[i,:] = ocp_solver.get(i, "u")
simX[N,:] = ocp_solver.get(N, "x")

ocp_solver.print_statistics() # encapsulates: stat = ocp_solver.get_stats("statistics")

plot_pendulum(np.linspace(0, Tf, N+1), Fmax, simU, simX, latexify=False)
