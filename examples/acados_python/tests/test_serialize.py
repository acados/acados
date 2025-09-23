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

from acados_template import AcadosOcp, AcadosOcpSolver, AcadosModel
from pendulum_model import export_pendulum_ode_model
import numpy as np
import scipy.linalg
from utils import plot_pendulum
import casadi as ca

def main(cost_type='NONLINEAR_LS', hessian_approximation='EXACT', ext_cost_use_num_hess=0,
         integrator_type='ERK'):
    print(f"using: cost_type {cost_type}, integrator_type {integrator_type}")

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
    N_horizon = 20

    ocp.solver_options.N_horizon = N_horizon

    # set cost
    Q_mat = 2*np.diag([1e3, 1e3, 1e-2, 1e-2])
    R_mat = 2*np.diag([1e-2])

    x = ocp.model.x
    u = ocp.model.u

    cost_W = scipy.linalg.block_diag(Q_mat, R_mat)

    if cost_type == 'LS':
        ocp.cost.cost_type = 'LINEAR_LS'
        ocp.cost.cost_type_e = 'LINEAR_LS'

        ocp.cost.Vx = np.zeros((ny, nx))
        ocp.cost.Vx[:nx,:nx] = np.eye(nx)

        Vu = np.zeros((ny, nu))
        Vu[4,0] = 1.0
        ocp.cost.Vu = Vu

        ocp.cost.Vx_e = np.eye(nx)

    elif cost_type == 'NONLINEAR_LS':
        ocp.cost.cost_type = 'NONLINEAR_LS'
        ocp.cost.cost_type_e = 'NONLINEAR_LS'

        ocp.model.cost_y_expr = ca.vertcat(x, u)
        ocp.model.cost_y_expr_e = x

    elif cost_type == 'EXTERNAL':
        ocp.cost.cost_type = 'EXTERNAL'
        ocp.cost.cost_type_e = 'EXTERNAL'

        ocp.model.cost_expr_ext_cost = ca.vertcat(x, u).T @ cost_W @ ca.vertcat(x, u)
        ocp.model.cost_expr_ext_cost_e = x.T @ Q_mat @ x

    else:
        raise Exception('Unknown cost_type. Possible values are \'LS\' and \'NONLINEAR_LS\'.')

    if cost_type in ['LS', 'NONLINEAR_LS']:
        ocp.cost.yref = np.zeros((ny, ))
        ocp.cost.yref_e = np.zeros((ny_e, ))
        ocp.cost.W_e = Q_mat
        ocp.cost.W = cost_W

    # set constraints
    Fmax = 80
    ocp.constraints.lbu = np.array([-Fmax])
    ocp.constraints.ubu = np.array([+Fmax])
    x0 = np.array([0.0, np.pi, 0.0, 0.0])
    ocp.constraints.x0 = x0
    ocp.constraints.idxbu = np.array([0])

    ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM'
    ocp.solver_options.hessian_approx = hessian_approximation
    ocp.solver_options.regularize_method = 'CONVEXIFY'
    ocp.solver_options.integrator_type = integrator_type

    if ocp.solver_options.integrator_type == 'GNSF':
        import json
        with open('../pendulum_on_cart/common/' + model.name + '_gnsf_functions.json', 'r') as f:
            gnsf_dict = json.load(f)
        ocp.gnsf_model = gnsf_dict

    # set prediction horizon
    ocp.solver_options.tf = Tf
    ocp.solver_options.nlp_solver_type = 'SQP'
    ocp.solver_options.ext_cost_num_hess = ext_cost_use_num_hess

    model_dict = ocp.model.to_dict()

    model_from_dict = AcadosModel.from_dict(model_dict)
    ocp.model = model_from_dict

    solver = AcadosOcpSolver(ocp, json_file=f"acados_ocp_{model.name}.json")



if __name__ == '__main__':
    for integrator_type in ['GNSF', 'ERK', 'IRK']:
        for cost_type in ['EXTERNAL', 'LS', 'NONLINEAR_LS']:
            hessian_approximation = 'GAUSS_NEWTON'
            ext_cost_use_num_hess = 1
            main(cost_type=cost_type, hessian_approximation=hessian_approximation,
                ext_cost_use_num_hess=ext_cost_use_num_hess, integrator_type=integrator_type)
