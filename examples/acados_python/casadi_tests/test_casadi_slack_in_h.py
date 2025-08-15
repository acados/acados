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
sys.path.insert(0, '../getting_started')
from acados_template import AcadosOcp, AcadosOcpSolver, AcadosCasadiOcpSolver
from pendulum_model import export_pendulum_ode_model
import numpy as np
import casadi as ca


N_horizon = 20
Tf = 1.0

def formulate_ocp(using_soft_constraints=True):
    # create ocp object to formulate the OCP
    ocp = AcadosOcp()

    # set model
    model = export_pendulum_ode_model()
    ocp.model = model

    nx = model.x.rows()
    nu = model.u.rows()

    # set prediction horizon
    ocp.solver_options.N_horizon = N_horizon
    ocp.solver_options.tf = Tf

    # cost matrices
    Q_mat = 2*np.diag([1e3, 1e3, 1e-2, 1e-2])
    R_mat = 2*np.diag([1e-2])

    # path cost
    ocp.cost.cost_type = 'NONLINEAR_LS'
    ocp.model.cost_y_expr = ca.vertcat(model.x, model.u)
    ocp.cost.yref = np.zeros((nx+nu,))
    ocp.cost.W = ca.diagcat(Q_mat, R_mat).full()

    # terminal cost
    ocp.cost.cost_type_e = 'NONLINEAR_LS'
    ocp.cost.yref_e = np.zeros((nx,))
    ocp.model.cost_y_expr_e = model.x
    ocp.cost.W_e = Q_mat

    # set constraints
    ocp.constraints.x0 = np.array([0, np.pi, 0, 0])  # initial state
    ocp.constraints.idxbx_0 = np.array([0, 1, 2, 3])

    if using_soft_constraints:
        xmax = 2
        vmax = 3
        Fmax = 50
        # soft bound on x, using constraint h
        v1 = ocp.model.x[2]
        x1 = ocp.model.x[0]
        u = ocp.model.u[0]

        # initial soft constraint on h
        ocp.model.con_h_expr_0 = ca.vertcat(u,v1)
        ocp.constraints.lh_0 = np.array([-Fmax, -vmax])
        ocp.constraints.uh_0 = np.array([+Fmax, +vmax])
        ocp.constraints.idxsh_0 = np.array([1])  # indices of slacked constraints within h
        #set initial penalty weight for slack variables
        ocp.cost.zl_0 = np.ones((1,))
        ocp.cost.Zl_0 = np.ones((1,))
        ocp.cost.zu_0 = np.ones((1,))
        ocp.cost.Zu_0 = np.ones((1,))

        # intermidiate soft constraints on h
        ocp.model.con_h_expr = ca.vertcat(x1, v1, u)
        ocp.constraints.lh = np.array([-xmax, -vmax, -Fmax])
        ocp.constraints.uh = np.array([+xmax, +vmax, +Fmax])
        ocp.constraints.idxsh = np.array([1]) # indices of slacked constraints within h
        # set penalty weight for slack variables
        ocp.cost.zl = np.ones((1,))
        ocp.cost.Zl = np.ones((1,))
        ocp.cost.zu = np.ones((1,))
        ocp.cost.Zu = np.ones((1,))

        # terminal soft constraint on h
        ocp.model.con_h_expr_e = ca.vertcat(v1, x1)
        ocp.constraints.lh_e = np.array([-vmax, -xmax])
        ocp.constraints.uh_e = np.array([+vmax, +xmax])
        ocp.constraints.idxsh_e = np.array([0,1]) # indices of slacked constraints within h
        # set penalty weight for terminal slack variable
        ocp.cost.zl_e = np.ones((2,))
        ocp.cost.Zl_e = np.ones((2,))
        ocp.cost.zu_e = np.ones((2,))
        ocp.cost.Zu_e = np.ones((2,))

    else:
        # hard constraints on u
        Fmax = 80
        ocp.constraints.lbu = np.array([-Fmax])
        ocp.constraints.ubu = np.array([+Fmax])
        ocp.constraints.idxbu = np.array([0])

    # set options
    ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM' # FULL_CONDENSING_QPOASES
    ocp.solver_options.hessian_approx = 'GAUSS_NEWTON' # 'GAUSS_NEWTON', 'EXACT'
    ocp.solver_options.integrator_type = 'ERK'
    ocp.solver_options.nlp_solver_type = 'SQP' # SQP_RTI, SQP
    ocp.solver_options.globalization = 'MERIT_BACKTRACKING' # turns on globalization

    return ocp

def main():
    ocp = formulate_ocp()

    # create solver
    ocp_solver = AcadosOcpSolver(ocp, verbose=False)
    # solve OCP
    status = ocp_solver.solve()
    if status != 0:
        raise Exception(f'acados OCP solver returned status {status}')
    result = ocp_solver.store_iterate_to_obj()
    print('acados_cost:', ocp_solver.get_cost())

    casadi_ocp_solver = AcadosCasadiOcpSolver(ocp)
    casadi_ocp_solver.load_iterate_from_obj(result)
    casadi_ocp_solver.solve()
    result_casadi = casadi_ocp_solver.store_iterate_to_obj()
    print('casadi_cost:', casadi_ocp_solver.get_cost())

    result.flatten().allclose(other=result_casadi.flatten())

    casadi_sqp_ocp_solver = AcadosCasadiOcpSolver(ocp, solver="sqpmethod", verbose=False)
    casadi_sqp_ocp_solver.load_iterate_from_obj(result)
    casadi_sqp_ocp_solver.solve()
    iteration = casadi_sqp_ocp_solver.get_stats('nlp_iter')
    if iteration > 1:
        raise Exception(f'casadi SQP solver returned {iteration} iterations, expected less than 1.')

if __name__ == '__main__':
    main()