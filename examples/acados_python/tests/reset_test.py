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
from casadi import vertcat


def main(cost_type='NONLINEAR_LS', hessian_approximation='EXACT', ext_cost_use_num_hess=0,
         integrator_type='ERK'):
    print(f"using: cost_type {cost_type}, integrator_type {integrator_type}")
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

    ocp.dims.N = N

    # set cost
    Q = 2*np.diag([1e3, 1e3, 1e-2, 1e-2])
    R = 2*np.diag([1e-2])

    x = ocp.model.x
    u = ocp.model.u

    cost_W = scipy.linalg.block_diag(Q, R)

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

        ocp.model.cost_y_expr = vertcat(x, u)
        ocp.model.cost_y_expr_e = x

    elif cost_type == 'EXTERNAL':
        ocp.cost.cost_type = 'EXTERNAL'
        ocp.cost.cost_type_e = 'EXTERNAL'

        ocp.model.cost_expr_ext_cost = vertcat(x, u).T @ cost_W @ vertcat(x, u)
        ocp.model.cost_expr_ext_cost_e = x.T @ Q @ x

    else:
        raise Exception('Unknown cost_type. Possible values are \'LS\' and \'NONLINEAR_LS\'.')

    if cost_type in ['LS', 'NONLINEAR_LS']:
        ocp.cost.yref = np.zeros((ny, ))
        ocp.cost.yref_e = np.zeros((ny_e, ))
        ocp.cost.W_e = Q
        ocp.cost.W = cost_W

    # set constraints
    Fmax = 80
    ocp.constraints.lbu = np.array([-Fmax])
    ocp.constraints.ubu = np.array([+Fmax])
    x0 = np.array([0.0, np.pi, 0.0, 0.0])
    ocp.constraints.x0 = x0
    ocp.constraints.idxbu = np.array([0])

    ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM' # FULL_CONDENSING_QPOASES
    ocp.solver_options.hessian_approx = hessian_approximation
    ocp.solver_options.regularize_method = 'CONVEXIFY'
    ocp.solver_options.integrator_type = integrator_type
    if ocp.solver_options.integrator_type == 'GNSF':
        import json
        with open('../pendulum_on_cart/common/' + model.name + '_gnsf_functions.json', 'r') as f:
            gnsf_dict = json.load(f)
        ocp.gnsf_model = gnsf_dict

    ocp.solver_options.qp_solver_cond_N = 5

    # set prediction horizon
    ocp.solver_options.tf = Tf
    ocp.solver_options.nlp_solver_type = 'SQP' # SQP_RTI
    ocp.solver_options.ext_cost_num_hess = ext_cost_use_num_hess

    ocp_solver = AcadosOcpSolver(ocp, json_file = 'acados_ocp.json')

    # set NaNs as input to test reset() -> NOT RECOMMENDED!!!
    # ocp_solver.options_set('print_level', 2)
    for i in range(N):
        ocp_solver.set(i, 'x', np.NaN * np.ones((nx,)))
        ocp_solver.set(i, 'u', np.NaN * np.ones((nu,)))
    status = ocp_solver.solve()
    ocp_solver.print_statistics() # encapsulates: stat = ocp_solver.get_stats("statistics")
    if status == 0:
        raise Exception(f'acados returned status {status}, although NaNs were given.')
    else:
        print(f'acados returned status {status}, which is expected, since NaNs were given.')

    # RESET
    ocp_solver.reset()
    for i in range(N):
        ocp_solver.set(i, 'x', x0)

    if cost_type == 'EXTERNAL':
        # NOTE: hessian is wrt [u,x]
        if ext_cost_use_num_hess:
            for i in range(N):
                ocp_solver.cost_set(i, "ext_cost_num_hess", np.diag([0.04, 4000, 4000, 0.04, 0.04, ]))
            ocp_solver.cost_set(N, "ext_cost_num_hess", np.diag([4000, 4000, 0.04, 0.04, ]))

    simX = np.ndarray((N+1, nx))
    simU = np.ndarray((N, nu))

    status = ocp_solver.solve()

    ocp_solver.print_statistics()
    if status != 0:
        raise Exception(f'acados returned status {status} for cost_type {cost_type}\n'
                        f'integrator_type = {integrator_type}.')

    # get solution
    for i in range(N):
        simX[i,:] = ocp_solver.get(i, "x")
        simU[i,:] = ocp_solver.get(i, "u")
    simX[N,:] = ocp_solver.get(N, "x")


if __name__ == '__main__':
    # for integrator_type in ['ERK', 'IRK']:
    # ['LIFTED_IRK']
    for integrator_type in ['GNSF', 'ERK', 'IRK']:
        for cost_type in ['EXTERNAL', 'LS', 'NONLINEAR_LS']:
            # for cost_type in ['EXTERNAL', 'LS', 'NONLINEAR_LS']:
            hessian_approximation = 'GAUSS_NEWTON' # 'GAUSS_NEWTON, EXACT
            ext_cost_use_num_hess = 1
            main(cost_type=cost_type, hessian_approximation=hessian_approximation,
                ext_cost_use_num_hess=ext_cost_use_num_hess, integrator_type=integrator_type)
