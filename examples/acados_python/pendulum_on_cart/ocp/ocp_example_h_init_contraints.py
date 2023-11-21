# -*- coding: future_fstrings -*-
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

from acados_template import AcadosOcp, AcadosOcpSolver, AcadosModel
from pendulum_model import export_pendulum_ode_model, export_augmented_pendulum_model
import numpy as np
import scipy.linalg
from utils import plot_pendulum
from casadi import vertcat, SX

# contraint to test and compare
# * non-linear constraint expression
# * non-linear constraint expression + relaxing the initial state constraints
# * state bounds constraint
CONSTRAINT_VERSIONS = ['nl', 'nl_relxd','bound']

def test_initial_h_constraints(constraint_version: str):
    print(f'#################################################################################### {constraint_version} constraint ####################################################################################')
    ocp = AcadosOcp()

    model = export_pendulum_ode_model()
    model.name = f'{model.name}_{constraint_version}_LS'
    # set model
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
    ocp.cost.W_e = Q
    ocp.cost.W = cost_W

    # set constraints
    Fmax = 50
    ocp.constraints.lbu = np.array([-Fmax])
    ocp.constraints.ubu = np.array([+Fmax])
    ocp.constraints.idxbu = np.array([0])

    ocp.constraints.x0 = np.array([-2, np.pi, 0.0, 0.0])

    lbx = np.array([-2., -np.pi, -4, -5])
    ubx = -lbx

    if constraint_version == 'bound':
        ocp.constraints.lbx = lbx
        ocp.constraints.ubx = ubx
        ocp.constraints.idxbx = np.arange(nx)
        ocp.constraints.lbx_e = lbx
        ocp.constraints.ubx_e = ubx
        ocp.constraints.idxbx_e = np.arange(nx)

    elif constraint_version == 'nl':
        # relaxed initial state h constraints
        lbh_0 = np.array([-2, -np.pi, -4, -5])
        ubh_0 = -lbh_0
        ocp.model.con_h_expr = ocp.model.x
        ocp.constraints.lh = lbx
        ocp.constraints.uh = ubx
        ocp.model.con_h_expr_0 = ocp.model.x
        ocp.constraints.lh_0 = lbh_0
        ocp.constraints.uh_0 = ubh_0
        ocp.model.con_h_expr_e = ocp.model.x
        ocp.constraints.lh_e = lbx
        ocp.constraints.uh_e = ubx

    elif constraint_version == 'nl_relxd':
        # relaxed initial state h constraints
        lbh_0 = 10 * np.array([-2, -np.pi, -4, -5])
        ubh_0 = -lbh_0
        ocp.model.con_h_expr = ocp.model.x
        ocp.constraints.lh = lbx
        ocp.constraints.uh = ubx
        ocp.model.con_h_expr_0 = ocp.model.x
        ocp.constraints.lh_0 = lbh_0
        ocp.constraints.uh_0 = ubh_0
        ocp.model.con_h_expr_e = ocp.model.x
        ocp.constraints.lh_e = lbx
        ocp.constraints.uh_e = ubx
    else:
        raise NotImplementedError(f'constraint {constraint_version} is not implemented in this example!')

    ocp.solver_options.qp_solver = 'FULL_CONDENSING_HPIPM'
    ocp.solver_options.hessian_approx = 'GAUSS_NEWTON'
    ocp.solver_options.integrator_type = 'IRK'

    ocp.solver_options.qp_solver_cond_N = 5

    # set prediction horizon
    ocp.solver_options.tf = Tf
    ocp.solver_options.nlp_solver_type = 'SQP'

    ocp_solver = AcadosOcpSolver(ocp, json_file = 'acados_ocp.json')

    simX = np.ndarray((N+1, nx))
    simU = np.ndarray((N, nu))

    status = ocp_solver.solve()
    if status != 0:
        print(f'acados returned status {status}.')

    # for debugging
    # ocp_solver.print_statistics()
    # ocp_solver.dump_last_qp_to_json(f'{constraint_version}_last_qp.json', overwrite=True)

    # get solution
    for i in range(N):
        simX[i,:] = ocp_solver.get(i, "x")
        simU[i,:] = ocp_solver.get(i, "u")
    simX[N,:] = ocp_solver.get(N, "x")
    # plot results
    plot_pendulum(np.linspace(0, Tf, N+1), Fmax, simU, simX, latexify=False)

    return status

if __name__ == "__main__":
    # we expect that the OCP with nl constraints will fail because it contains two active constraints at the initial node.
    expected_status = [2, 0, 0]
    for constraint_version in CONSTRAINT_VERSIONS:
        if test_initial_h_constraints(constraint_version=constraint_version) != expected_status[CONSTRAINT_VERSIONS.index(constraint_version)]:
            raise Exception(f'constraint {constraint_version} failed!')
        else:
            print(f'constraint {constraint_version} passed!')
