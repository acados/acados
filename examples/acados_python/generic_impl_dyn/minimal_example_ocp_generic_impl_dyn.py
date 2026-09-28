#
# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# Licensed under the 2-Clause BSD License.

import sys
sys.path.insert(0, '../pendulum_on_cart/common')

from acados_template import AcadosOcp, AcadosOcpSolver
from pendulum_model import export_pendulum_ode_model
import numpy as np
import scipy.linalg
from utils import plot_pendulum

# create ocp object to formulate the OCP
ocp = AcadosOcp()

# set model
model = export_pendulum_ode_model()
ocp.model = model

Tf = 1.0
nx = model.x.rows()
nu = model.u.rows()
ny = nx + nu
ny_e = nx
N = 20


# GENERIC
ocp.model.dyn_ext_fun_type = 'generic'
ocp.model.dyn_generic_source = 'generic_impl_dyn.c'
ocp.model.dyn_impl_dae_fun = 'generic_impl_dyn_fun'
ocp.model.dyn_impl_dae_fun_jac = 'generic_impl_dyn_fun_jac'
ocp.model.dyn_impl_dae_jac = 'generic_impl_dyn_jac'

# set dimensions
ocp.solver_options.N_horizon = N

# set cost
Q = 2*np.diag([1e3, 1e3, 1e-2, 1e-2])
R = 2*np.diag([1e-2])

ocp.cost.W_e = Q
ocp.cost.W = scipy.linalg.block_diag(Q, R)

ocp.cost.cost_type = 'LINEAR_LS'
ocp.cost.cost_type_e = 'LINEAR_LS'

ocp.cost.Vx = np.zeros((ny, nx))
ocp.cost.Vx[:nx,:nx] = np.eye(nx)

Vu = np.zeros((ny, nu))
Vu[4,0] = 1.0
ocp.cost.Vu = Vu

ocp.cost.Vx_e = np.eye(nx)

ocp.cost.yref  = np.zeros((ny, ))
ocp.cost.yref_e = np.zeros((ny_e, ))

# set constraints
Fmax = 80
ocp.constraints.lbu = np.array([-Fmax])
ocp.constraints.ubu = np.array([+Fmax])
ocp.constraints.idxbu = np.array([0])

ocp.constraints.x0 = np.array([0.0, np.pi, 0.0, 0.0])

# set options
ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM' # FULL_CONDENSING_QPOASES
# PARTIAL_CONDENSING_HPIPM, FULL_CONDENSING_QPOASES, FULL_CONDENSING_HPIPM,
# PARTIAL_CONDENSING_QPDUNES, PARTIAL_CONDENSING_OSQP
ocp.solver_options.hessian_approx = 'GAUSS_NEWTON'
ocp.solver_options.integrator_type = 'IRK'
ocp.solver_options.print_level = 0
ocp.solver_options.nlp_solver_type = 'SQP' # SQP_RTI, SQP

# set prediction horizon
ocp.solver_options.tf = Tf
ocp.solver_options.nlp_solver_max_iter = 200
ocp.solver_options.sim_method_num_stages = 1
ocp.solver_options.sim_method_num_steps = 1



ocp_solver = AcadosOcpSolver(ocp, json_file = 'acados_ocp.json')

simX = np.zeros((N+1, nx))
simU = np.zeros((N, nu))

status = ocp_solver.solve()

if status != 0:
    ocp_solver.print_statistics() # encapsulates: stat = ocp_solver.get_stats("statistics")
    raise Exception(f'acados returned status {status}.')

# get solution
for i in range(N):
    simX[i,:] = ocp_solver.get(i, "x")
    simU[i,:] = ocp_solver.get(i, "u")
simX[N,:] = ocp_solver.get(N, "x")

ocp_solver.print_statistics() # encapsulates: stat = ocp_solver.get_stats("statistics")

plot_pendulum(np.linspace(0, Tf, N+1), Fmax, simU, simX, latexify=False)
