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
from pendulum_model import export_pendulum_ode_model

from acados_template.mpc_utils import detect_constraint_structure
import casadi as ca
import numpy as np
from acados_template import AcadosOcp, AcadosOcpSolver


def main(constr_to_test: str):
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
    if constr_to_test == "bounds":
        ocp.model.con_h_expr_0 = ocp.model.u
        ocp.model.con_h_expr = ocp.model.u
        ocp.constraints.lh = np.array([-Fmax])
        ocp.constraints.uh = np.array([+Fmax])
        ocp.constraints.lh_0 = np.array([-Fmax])
        ocp.constraints.uh_0 = np.array([+Fmax])
        detect_constraint_structure(ocp.model, ocp.constraints, "initial")
        detect_constraint_structure(ocp.model, ocp.constraints, "path")
    elif constr_to_test == "linear":
        ocp.model.con_h_expr_0 = ocp.model.u + \
            1e-9 * ca.DM(np.ones(4)).T @ ocp.model.x
        ocp.model.con_h_expr = ocp.model.u + \
            1e-9 * ca.DM(np.ones(4)).T @ ocp.model.x
        ocp.constraints.lh = np.array([-Fmax])
        ocp.constraints.uh = np.array([+Fmax])
        ocp.constraints.lh_0 = np.array([-Fmax])
        ocp.constraints.uh_0 = np.array([+Fmax])
        detect_constraint_structure(ocp.model, ocp.constraints, "initial")
        detect_constraint_structure(ocp.model, ocp.constraints, "path")
    elif constr_to_test == "ground_truth":
        ocp.constraints.lbu = np.array([-Fmax])
        ocp.constraints.ubu = np.array([+Fmax])
        ocp.constraints.idxbu = np.array([0])

    ocp.constraints.x0 = np.array([0.0, np.pi, 0.0, 0.0])

    # set options
    ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM'  # FULL_CONDENSING_QPOASES
    # PARTIAL_CONDENSING_HPIPM, FULL_CONDENSING_QPOASES, FULL_CONDENSING_HPIPM,
    # PARTIAL_CONDENSING_QPDUNES, PARTIAL_CONDENSING_OSQP, FULL_CONDENSING_DAQP
    ocp.solver_options.hessian_approx = 'GAUSS_NEWTON'  # 'GAUSS_NEWTON', 'EXACT'
    ocp.solver_options.integrator_type = 'IRK'
    # ocp.solver_options.print_level = 1
    ocp.solver_options.nlp_solver_type = 'SQP'  # SQP_RTI, SQP
    ocp.solver_options.globalization = 'MERIT_BACKTRACKING'  # turns on globalization

    ocp_solver = AcadosOcpSolver(ocp)

    simX = np.zeros((N+1, nx))
    simU = np.zeros((N, nu))

    status = ocp_solver.solve()
    # encapsulates: stat = ocp_solver.get_stats("statistics")
    ocp_solver.print_statistics()

    if status != 0:
        raise Exception(f'acados returned status {status}.')

    # get solution
    for i in range(N):
        simX[i, :] = ocp_solver.get(i, "x")
        simU[i, :] = ocp_solver.get(i, "u")
    simX[N, :] = ocp_solver.get(N, "x")

    return simX, simU


if __name__ == '__main__':
    simX_1, simU_1 = main(constr_to_test="ground_truth")
    simX_2, simU_2 = main(constr_to_test="bounds")
    simX_3, simU_3 = main(constr_to_test="linear")

    for i in range(simX_1.shape[1]):
        assert (np.allclose(simX_1[i, :], simX_2[i, :], atol=1e-6, rtol=0))
        assert (np.allclose(simX_1[i, :], simX_3[i, :], atol=1e-6, rtol=0))
    assert (np.allclose(simU_1, simU_2, atol=1e-6, rtol=0))
    assert (np.allclose(simU_1, simU_3, atol=1e-6, rtol=0))
