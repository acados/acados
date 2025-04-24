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

import casadi as ca


def main(qp_solver: str = 'PARTIAL_CONDENSING_HPIPM'):
    # create ocp object to formulate the OCP
    ocp = AcadosOcp()

    # set model
    model = export_pendulum_ode_model()
    ocp.model = model

    Tf = 1.0
    nx = model.x.rows()
    nu = model.u.rows()
    N = 20

    # set prediction horizon
    ocp.solver_options.N_horizon = N
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
    Fmax = 80
    ocp.constraints.lbu = np.array([-Fmax])
    ocp.constraints.ubu = np.array([+Fmax])
    ocp.constraints.idxbu = np.array([0])

    ocp.constraints.x0 = np.array([0.0, np.pi, 0.0, 0.0])

    # set options
    ocp.solver_options.qp_solver = qp_solver
    ocp.solver_options.hessian_approx = 'GAUSS_NEWTON'
    ocp.solver_options.integrator_type = 'IRK'
    ocp.solver_options.nlp_solver_type = 'SQP' # SQP_RTI, SQP
    ocp.solver_options.globalization = 'MERIT_BACKTRACKING'
    ocp.solver_options.qp_tol = 1e-8

    ocp_solver = AcadosOcpSolver(ocp, verbose=False)

    status = ocp_solver.solve()
    ocp_solver.print_statistics() # encapsulates: stat = ocp_solver.get_stats("statistics")

    if status != 0:
        raise Exception(f'acados returned status {status}.')

    sol = ocp_solver.store_iterate_to_obj()
    qp_iter = ocp_solver.get_stats("qp_iter")
    nlp_iter = ocp_solver.get_stats("nlp_iter")
    assert nlp_iter == 10, f"cold start should require 10 iterations, got {nlp_iter}"

    ocp_solver.load_iterate_from_obj(sol)
    ocp_solver.solve()
    ocp_solver.print_statistics()
    nlp_iter = ocp_solver.get_stats("nlp_iter")
    assert nlp_iter == 0, f"hot start should require 0 iterations, got {nlp_iter}"

    disturbed_sol = sol
    disturbed_sol.x_traj[0] += 0.1 * disturbed_sol.x_traj[0]
    ocp_solver.load_iterate_from_obj(disturbed_sol)
    status = ocp_solver.solve()
    ocp_solver.print_statistics()
    nlp_iter = ocp_solver.get_stats("nlp_iter")
    assert nlp_iter == 6, f"warm start should require 6 iterations, got {nlp_iter}"

    ocp_solver.load_iterate_from_obj(disturbed_sol)
    ocp_solver.options_set("qp_warm_start", 1)
    ocp_solver.options_set("warm_start_first_qp_from_nlp", True)
    ocp_solver.options_set("warm_start_first_qp", True)
    # ocp_solver.options_set("qp_mu0", 1e-3)
    status = ocp_solver.solve()
    ocp_solver.print_statistics()


if __name__ == '__main__':
    main('FULL_CONDENSING_HPIPM')
    main('PARTIAL_CONDENSING_HPIPM')
