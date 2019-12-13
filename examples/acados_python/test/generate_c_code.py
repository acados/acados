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
import json
import argparse

# set to 'True' to generate test data
GENERATE_DATA = False

LOCAL_TEST = False
TEST_TOL = 1e-8

if LOCAL_TEST is True:
    FORMULATION = 'LS'
    SOLVER_TYPE = 'SQP_RTI'
    QP_SOLVER = 'FULL_CONDENSING_QPOASES'
    INTEGRATOR_TYPE = 'IRK'
else:
    parser = argparse.ArgumentParser(description='test Python interface on pendulum example.')
    parser.add_argument('--FORMULATION', dest='FORMULATION',
                        default='LS',
                        help='FORMULATION: linear least-squares (LS) or nonlinear \
                                least-squares (NLS) (default: LS)')

    parser.add_argument('--QP_SOLVER', dest='QP_SOLVER',
                        default='PARTIAL_CONDENSING_HPIPM',
                        help='QP_SOLVER: PARTIAL_CONDENSING_HPIPM, FULL_CONDENSING_HPIPM, ' \
                                'FULL_CONDENSING_HPIPM (default: PARTIAL_CONDENSING_HPIPM)')

    parser.add_argument('--INTEGRATOR_TYPE', dest='INTEGRATOR_TYPE',
                        default='ERK',
                        help='INTEGRATOR_TYPE: explicit (ERK) or implicit (IRK) ' \
                                ' Runge-Kutta (default: ERK)')

    parser.add_argument('--SOLVER_TYPE', dest='SOLVER_TYPE',
                        default='SQP_RTI',
                        help='SOLVER_TYPE: (full step) sequential quadratic programming (SQP) or ' \
                                ' real-time iteration (SQP-RTI) (default: SQP-RTI)')


    args = parser.parse_args()

    FORMULATION = args.FORMULATION
    FORMULATION_values = ['LS', 'NLS']
    if FORMULATION not in FORMULATION_values:
        raise Exception('Invalid unit test value {} for parameter FORMULATION. Possible values are' \
                ' {}. Exiting.'.format(FORMULATION, FORMULATION_values))

    QP_SOLVER = args.QP_SOLVER
    QP_SOLVER_values = ['PARTIAL_CONDENSING_HPIPM', 'FULL_CONDENSING_HPIPM', 'FULL_CONDENSING_QPOASES']
    if QP_SOLVER not in QP_SOLVER_values:
        raise Exception('Invalid unit test value {} for parameter QP_SOLVER. Possible values are' \
                ' {}. Exiting.'.format(QP_SOLVER, QP_SOLVER_values))

    INTEGRATOR_TYPE = args.INTEGRATOR_TYPE
    INTEGRATOR_TYPE_values = ['ERK', 'IRK']
    if INTEGRATOR_TYPE not in INTEGRATOR_TYPE:
        raise Exception('Invalid unit test value {} for parameter INTEGRATOR_TYPE. Possible values are' \
                ' {}. Exiting.'.format(INTEGRATOR_TYPE, INTEGRATOR_TYPE_values))

    SOLVER_TYPE = args.SOLVER_TYPE
    SOLVER_TYPE_values = ['SQP', 'SQP-RTI']
    if SOLVER_TYPE not in SOLVER_TYPE:
        raise Exception('Invalid unit test value {} for parameter SOLVER_TYPE. Possible values are' \
                ' {}. Exiting.'.format(SOLVER_TYPE, SOLVER_TYPE_values))


# print test setting
print("Running test with:\n\tformulation:", FORMULATION, "\n\tqp solver: ", QP_SOLVER,\
      "\n\tintergrator: ", INTEGRATOR_TYPE, "\n\tsolver: ", SOLVER_TYPE)

# create render arguments
ocp = acados_ocp_nlp()

# export model
model = export_ode_model()

# set model_name
ocp.model = model

Tf = 2.0
nx = model.x.size()[0]
nu = model.u.size()[0]
ny = nx + nu
ny_e = nx
N = 50

# set ocp_nlp_dimensions
nlp_dims     = ocp.dims
nlp_dims.nx  = nx
nlp_dims.ny  = ny
nlp_dims.ny_e = ny_e
nlp_dims.nbx = 0
nlp_dims.nbu = nu
nlp_dims.nu  = model.u.size()[0]
nlp_dims.N   = N

# set weighting matrices
nlp_cost = ocp.cost

if FORMULATION == 'LS':
    nlp_cost.cost_type = 'LINEAR_LS'
    nlp_cost.cost_type_e = 'LINEAR_LS'
elif FORMULATION == 'NLS':
    nlp_cost.cost_type = 'NONLINEAR_LS'
    nlp_cost.cost_type_e = 'NONLINEAR_LS'
else:
    raise Exception('Unknown FORMULATION. Possible values are \'LS\' and \'NLS\'.')

Q = np.eye(4)
Q[0,0] = 1e0
Q[1,1] = 1e2
Q[2,2] = 1e-3
Q[3,3] = 1e-2

R = np.eye(1)
R[0,0] = 1e0

unscale = N/Tf
Q = Q * unscale
R = R * unscale

if FORMULATION == 'NLS':
    nlp_cost.W = scipy.linalg.block_diag(R, Q)
else:
    nlp_cost.W = scipy.linalg.block_diag(Q, R)

nlp_cost.W_e = Q/unscale

Vx = np.zeros((ny, nx))
Vx[0,0] = 1.0
Vx[1,1] = 1.0
Vx[2,2] = 1.0
Vx[3,3] = 1.0

nlp_cost.Vx = Vx

Vu = np.zeros((ny, nu))
Vu[4,0] = 1.0
nlp_cost.Vu = Vu


Vx_e = np.zeros((ny_e, nx))
Vx_e[0,0] = 1.0
Vx_e[1,1] = 1.0
Vx_e[2,2] = 1.0
Vx_e[3,3] = 1.0

nlp_cost.Vx_e = Vx_e
if FORMULATION == 'NLS':
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


nlp_cost.yref  = np.zeros((ny, ))
nlp_cost.yref_e = np.zeros((ny_e, ))

# setting bounds
Fmax = 2.0
nlp_con = ocp.constraints
nlp_con.lbu = np.array([-Fmax])
nlp_con.ubu = np.array([+Fmax])
nlp_con.x0 = np.array([0.0, 3.14, 0.0, 0.0])
nlp_con.idxbu = np.array([0])

# set QP solver
ocp.solver_options.qp_solver = QP_SOLVER
ocp.solver_options.hessian_approx = 'GAUSS_NEWTON'
ocp.solver_options.integrator_type = INTEGRATOR_TYPE
ocp.solver_options.sim_method_num_stages = 2
ocp.solver_options.sim_method_num_steps = 5

# set prediction horizon
ocp.solver_options.tf = Tf
ocp.solver_options.nlp_solver_type = SOLVER_TYPE

# set header path
ocp.acados_include_path  = '../../../../include'
ocp.acados_lib_path      = '../../../../lib'

acados_solver = generate_solver(ocp, json_file = 'acados_ocp.json')

Nsim = 100

simX = np.ndarray((Nsim, nx))
simU = np.ndarray((Nsim, nu))

for i in range(Nsim):
    status = acados_solver.solve()

    if status !=0:
        print("acados failure! Exiting. \n")
        sys.exit(status)

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

# dump result to JSON file for unit testing
test_file_name = 'test_data/generate_c_code_out_' + FORMULATION + '_' + QP_SOLVER + '_' + \
            INTEGRATOR_TYPE + '_' + SOLVER_TYPE + '.json'

if GENERATE_DATA:
    with open(test_file_name, 'w') as f:
        json.dump({"simX": simX.tolist(), "simU": simU.tolist()}, f, indent=4, sort_keys=True)
else:
    with open(test_file_name, 'r') as f:
        test_data = json.load(f)
    simX_error = np.linalg.norm(test_data['simX'] - simX)
    simU_error = np.linalg.norm(test_data['simU'] - simU)
    if  simX_error > TEST_TOL or  simU_error > TEST_TOL:
        raise Exception("Python acados test failure with accuracies {:.2E} and {:.2E} ({:.2E} required) on pendulum example! Exiting.\n".format(simX_error, simU_error, TEST_TOL))
    else:
        print('Python test passed with accuracy {:.2E}'.format(max(simU_error, simX_error)))
