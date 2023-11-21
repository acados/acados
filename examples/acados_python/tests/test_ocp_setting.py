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

from acados_template import *
from pendulum_model import export_pendulum_ode_model
from casadi import vertcat
import json
import numpy as np
import scipy.linalg
import argparse

# set to 'True' to generate test data
GENERATE_DATA = False

LOCAL_TEST = False
TEST_TOL = 1e-8

if LOCAL_TEST is True:
    COST_MODULE = 'LS'
    COST_MODULE_N = 'LS'
    SOLVER_TYPE = 'SQP_RTI'
    QP_SOLVER = 'FULL_CONDENSING_QPOASES'
    INTEGRATOR_TYPE = 'IRK'
    HESS_APPROX = 'GAUSS_NEWTON'
else:
    parser = argparse.ArgumentParser(description='test Python interface on pendulum example.')
    parser.add_argument('--COST_MODULE', dest='COST_MODULE',
                        default='LS',
                        help='COST_MODULE: linear least-squares (LS) or nonlinear \
                                least-squares (NLS) (default: LS), external (EXTERNAL)')

    parser.add_argument('--COST_MODULE_N', dest='COST_MODULE_N',
                        default='LS',
                        help='COST_MODULE_N: linear least-squares (LS) or nonlinear \
                                least-squares (NLS) (default: LS), external (EXTERNAL)')

    parser.add_argument('--QP_SOLVER', dest='QP_SOLVER',
                        default='PARTIAL_CONDENSING_HPIPM',
                        help='QP_SOLVER: PARTIAL_CONDENSING_HPIPM, FULL_CONDENSING_HPIPM, ' \
                                'FULL_CONDENSING_HPIPM (default: PARTIAL_CONDENSING_HPIPM)')

    parser.add_argument('--INTEGRATOR_TYPE', dest='INTEGRATOR_TYPE',
                        default='ERK',
                        help='INTEGRATOR_TYPE: explicit (ERK) or implicit (IRK) or GNSF-IRK (GNSF) ' \
                                ' Runge-Kutta (default: ERK)')

    parser.add_argument('--SOLVER_TYPE', dest='SOLVER_TYPE',
                        default='SQP_RTI',
                        help='SOLVER_TYPE: (full step) sequential quadratic programming (SQP) or ' \
                                ' real-time iteration (SQP_RTI) (default: SQP_RTI)')

    parser.add_argument('--HESS_APPROX', dest='HESS_APPROX',
                        default='GAUSS_NEWTON',
                        help='HESS_APPROX: GAUSS_NEWTON or ' \
                                ' EXACT (default: GAUSS_NEWTON)')

    parser.add_argument('--REGULARIZATION', dest='REGULARIZATION',
                        default='NO_REGULARIZE',
                        help='REGULARIZATION: NO_REGULARIZE or MIRROR or PROJECT or CONVEXIFY' \
                                ' or PROJECT_REDUC_HESS (default: NO_REGULARIZE)')

    args = parser.parse_args()

    COST_MODULE = args.COST_MODULE
    COST_MODULE_values = ['LS', 'NLS', 'EXTERNAL']
    if COST_MODULE not in COST_MODULE_values:
        raise Exception('Invalid unit test value {} for parameter COST_MODULE. Possible values are' \
                ' {}. Exiting.'.format(COST_MODULE, COST_MODULE_values))

    COST_MODULE_N = args.COST_MODULE_N
    COST_MODULE_N_values = ['LS', 'NLS', 'EXTERNAL']
    if COST_MODULE_N not in COST_MODULE_N_values:
        raise Exception('Invalid unit test value {} for parameter COST_MODULE_N. Possible values are' \
                ' {}. Exiting.'.format(COST_MODULE_N, COST_MODULE_N_values))

    QP_SOLVER = args.QP_SOLVER
    QP_SOLVER_values = ['PARTIAL_CONDENSING_HPIPM', 'FULL_CONDENSING_HPIPM', 'FULL_CONDENSING_QPOASES']
    if QP_SOLVER not in QP_SOLVER_values:
        raise Exception('Invalid unit test value {} for parameter QP_SOLVER. Possible values are' \
                ' {}. Exiting.'.format(QP_SOLVER, QP_SOLVER_values))

    INTEGRATOR_TYPE = args.INTEGRATOR_TYPE
    INTEGRATOR_TYPE_values = ['ERK', 'IRK', 'GNSF']
    if INTEGRATOR_TYPE not in INTEGRATOR_TYPE_values:
        raise Exception('Invalid unit test value {} for parameter INTEGRATOR_TYPE. Possible values are' \
                ' {}. Exiting.'.format(INTEGRATOR_TYPE, INTEGRATOR_TYPE_values))

    SOLVER_TYPE = args.SOLVER_TYPE
    SOLVER_TYPE_values = ['SQP', 'SQP_RTI']
    if SOLVER_TYPE not in SOLVER_TYPE_values:
        raise Exception('Invalid unit test value {} for parameter SOLVER_TYPE. Possible values are' \
                ' {}. Exiting.'.format(SOLVER_TYPE, SOLVER_TYPE_values))

    HESS_APPROX = args.HESS_APPROX
    HESS_APPROX_values = ['GAUSS_NEWTON', 'EXACT']
    if HESS_APPROX not in HESS_APPROX_values:
        raise Exception('Invalid unit test value {} for parameter HESS_APPROX. Possible values are' \
                ' {}. Exiting.'.format(HESS_APPROX, HESS_APPROX_values))

    REGULARIZATION = args.REGULARIZATION
    REGULARIZATION_values = ['NO_REGULARIZE', 'MIRROR', 'PROJECT', 'PROJECT_REDUC_HESS', 'CONVEXIFY']
    if REGULARIZATION not in REGULARIZATION_values:
        raise Exception('Invalid unit test value {} for parameter REGULARIZATION. Possible values are' \
                ' {}. Exiting.'.format(REGULARIZATION, REGULARIZATION_values))

# print test setting
print("Running test with:", \
      "\n\tcost module terminal:", COST_MODULE_N,\
      "\n\tcost module:", COST_MODULE, \
      "\n\thessian approximation:", HESS_APPROX, \
      "\n\tintergrator:", INTEGRATOR_TYPE, \
      "\n\tqp solver:", QP_SOLVER,\
      "\n\tregularization:", REGULARIZATION, \
      "\n\tsolver:", SOLVER_TYPE)

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
ocp.dims.N = N

# set cost
Q = 2*np.diag([1e3, 1e3, 1e-2, 1e-2])
R = 2*np.diag([1e-2])

x = ocp.model.x
u = ocp.model.u

Vx = np.zeros((ny, nx))
Vx[:nx,:nx] = np.eye(nx)

Vu = np.zeros((ny, nu))
Vu[4,0] = 1.0

if COST_MODULE == 'LS':
    ocp.cost.cost_type = 'LINEAR_LS'

    ocp.cost.Vx = Vx
    ocp.cost.Vu = Vu

    ocp.cost.W = scipy.linalg.block_diag(Q, R)
    ocp.cost.yref  = np.zeros((ny, ))

elif COST_MODULE == 'NLS':
    ocp.cost.cost_type = 'NONLINEAR_LS'
    ocp.model.cost_y_expr = vertcat(x, u)

    ocp.cost.W = scipy.linalg.block_diag(Q, R)
    ocp.cost.yref  = np.zeros((ny, ))

elif COST_MODULE == 'EXTERNAL':
    ocp.cost.cost_type = 'EXTERNAL'
    ocp.model.cost_expr_ext_cost = 0.5 * vertcat(x, u).T @ scipy.linalg.block_diag(Q, R) @ vertcat(x, u)

else:
    raise Exception('Unknown COST_MODULE. Possible values are \'LS\', \'NLS\', \'EXTERNAL\'.')


if COST_MODULE_N == 'LS':
    ocp.cost.cost_type_e = 'LINEAR_LS'
    ocp.cost.Vx_e = np.eye(nx)
    ocp.cost.W_e = Q
    ocp.cost.yref_e = np.zeros((ny_e, ))

elif COST_MODULE_N == 'NLS':
    ocp.cost.cost_type_e = 'NONLINEAR_LS'
    ocp.model.cost_y_expr_e = x
    ocp.cost.W_e = Q
    ocp.cost.yref_e = np.zeros((ny_e, ))

elif COST_MODULE_N == 'EXTERNAL':
    ocp.cost.cost_type_e = 'EXTERNAL'
    ocp.model.cost_expr_ext_cost_e = 0.5 * x.T @ Q @ x

else:
    raise Exception('Unknown COST_MODULE_N. Possible values are \'LS\', \'NLS\', \'EXTERNAL\'.')



# set constraints
Fmax = 80
ocp.constraints.lbu = np.array([-Fmax])
ocp.constraints.ubu = np.array([+Fmax])
ocp.constraints.x0 = np.array([0.0, np.pi, 0.0, 0.0])
ocp.constraints.idxbu = np.array([0])

# set options
ocp.solver_options.qp_solver = QP_SOLVER
ocp.solver_options.hessian_approx = HESS_APPROX
ocp.solver_options.exact_hess_constr = 1
ocp.solver_options.exact_hess_cost = 1
ocp.solver_options.exact_hess_dyn = 1
ocp.solver_options.regularize_method = REGULARIZATION

ocp.solver_options.integrator_type = INTEGRATOR_TYPE
ocp.solver_options.sim_method_num_stages = 2
ocp.solver_options.sim_method_num_steps = 5
ocp.solver_options.sim_method_newton_iter = 3

ocp.solver_options.nlp_solver_tol_stat = TEST_TOL
ocp.solver_options.nlp_solver_tol_eq = TEST_TOL
ocp.solver_options.nlp_solver_tol_ineq = TEST_TOL
ocp.solver_options.nlp_solver_tol_comp = TEST_TOL

ocp.solver_options.qp_solver_cond_N = int(N/2)

ocp.solver_options.nlp_solver_max_iter = 200
ocp.solver_options.qp_solver_iter_max = 50
ocp.solver_options.print_level = 0

ocp.solver_options.ext_fun_compile_flags = ''

# set prediction horizon
ocp.solver_options.tf = Tf
ocp.solver_options.nlp_solver_type = SOLVER_TYPE

if ocp.solver_options.integrator_type == 'GNSF':
    with open('../pendulum_on_cart/common/' + model.name + '_gnsf_functions.json', 'r') as f:
        gnsf_dict = json.load(f)
    ocp.gnsf_model = gnsf_dict

ocp_solver = AcadosOcpSolver(ocp, json_file = 'acados_ocp.json')

# initialize solver
x_traj_init = np.transpose( np.vstack( [np.zeros((N+1,)), \
     np.arange(np.pi, -np.pi/N,- np.pi/N), np.zeros((N+1,)), np.zeros((N+1,))]) )
for i in range(N+1):
    ocp_solver.set(i, "x", x_traj_init[i])

pi_init = np.ones((N, nx))
for i in range(N):
    ocp_solver.set(i, "pi", pi_init[i])

u_init = np.zeros((N, nu))
for i in range(N):
    ocp_solver.set(i, "u", u_init[i])

# solve ocp
simX = np.ndarray((N+1, nx))
simU = np.ndarray((N, nu))

status = ocp_solver.solve()

ocp_solver.print_statistics()

if status != 0:
    # import pdb; pdb.set_trace()
    raise Exception(f'acados returned status {status}.')

sqp_iter = ocp_solver.get_stats('sqp_iter')
if SOLVER_TYPE in {'SQP'}:
    print("Problem solved: SQP iterations ", sqp_iter, "\n")

# get solution
for i in range(N):
    simX[i,:] = ocp_solver.get(i, "x")
    simU[i,:] = ocp_solver.get(i, "u")
simX[N,:] = ocp_solver.get(N, "x")

if COST_MODULE in {'LINEAR_LS', 'NONLINEAR_LS'}:
    # update reference
    for j in range(N):
        ocp_solver.cost_set(j, "yref", np.array([0, 0, 0, 0, 0]))
    ocp_solver.cost_set(N, "yref", np.array([0, 0, 0, 0]))

# dump result to JSON file for unit testing
test_file_name = 'test_data/pendulum_ocp_formulations/test_ocp_' + COST_MODULE + \
                    '_' + COST_MODULE_N + '_' + QP_SOLVER + \
                    '_' + INTEGRATOR_TYPE + \
                    '_' + SOLVER_TYPE + \
                    '_' + HESS_APPROX + \
                    '_' + REGULARIZATION + \
                    '.json'

if GENERATE_DATA:
    with open(test_file_name, 'w') as f:
        json.dump({"simX": simX.tolist(), "simU": simU.tolist()}, f, indent=4, sort_keys=True)
else:
    ## Perform checks
    # check trajectory against known results
    with open(test_file_name, 'r') as f:
        test_data = json.load(f)
    simX_error = np.linalg.norm(test_data['simX'] - simX)
    simU_error = np.linalg.norm(test_data['simU'] - simU)

    CHECK_TOL = 50 * TEST_TOL
    if simX_error > CHECK_TOL or simU_error > CHECK_TOL:
        raise Exception("Python acados test failure with accuracies" +
                        " {:.2E} and {:.2E} ({:.2E} required)".format(simX_error, simU_error, CHECK_TOL) +
                        " on pendulum example! Exiting.\n")
    else:
        print('Python test passed with accuracy {:.2E}'.format(max(simU_error, simX_error)))

    # check cost function computation
    acados_cost = ocp_solver.get_cost()

    total_cost = 0
    W = scipy.linalg.block_diag(Q, R)

    for i in range(N):
        x = ocp_solver.get(i, "x")
        u = ocp_solver.get(i, "u")

        y = Vx @ x + Vu @ u
        cost =  0.5 * (y.T @ W @ y)
        # print("stage ", i, " y ", y, " unscaled cost acados style ", cost, "scaled ", cost * Tf/N)
        total_cost += Tf/N * cost

    i = N
    x = ocp_solver.get(i, "x")
    cost = .5 * (x.T @ Q @ x)
    # print("stage ", i, " x ", x, " unscaled cost acados style ", cost, "scaled ", cost * Tf/N)
    total_cost += cost

    cost_err = np.abs(acados_cost - total_cost)/ np.max((acados_cost, total_cost))

    if cost_err < 1e-10:
        print("Passed cost computation test")
        print("cost acados: ", acados_cost, " manually computed cost ", total_cost, "relative error", cost_err)
    else:
        raise Exception("Cost function computation test failed \n" +
            "cost acados: {:.2E}, manually computed cost {:.2E} relative error {:.2E}.".format(acados_cost, total_cost, cost_err))