# This test is an extension of the 'minimal_example_ocp_reuse_code.py' example.
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

from acados_template import AcadosOcp, AcadosOcpSolver, ACADOS_INFTY
from pendulum_model import export_pendulum_ode_model
import numpy as np
import scipy.linalg
from utils import plot_pendulum

PLOT = False

def main(constraint_variant='one_sided',
         qp_solver='PARTIAL_CONDENSING_HPIPM'):

    # create ocp object to formulate the OCP
    ocp = AcadosOcp()

    # set model
    model = export_pendulum_ode_model()
    ocp.model = model

    nx = model.x.rows()
    nu = model.u.rows()
    ny = nx + nu
    ny_e = nx

    N_horizon = 20  # number of shooting nodes
    Tf = 1.0

    # set dimensions
    ocp.solver_options.N_horizon = N_horizon

    # set cost
    Q = 2 * np.diag([1e3, 1e3, 1e-2, 1e-2])
    R = 2 * np.diag([1e-2])

    ocp.cost.W_e = Q
    ocp.cost.W = scipy.linalg.block_diag(Q, R)

    ocp.cost.cost_type = 'LINEAR_LS'
    ocp.cost.cost_type_e = 'LINEAR_LS'

    ocp.cost.Vx = np.zeros((ny, nx))
    ocp.cost.Vx[:nx, :nx] = np.eye(nx)

    Vu = np.zeros((ny, nu))
    Vu[4, 0] = 1.0
    ocp.cost.Vu = Vu

    ocp.cost.Vx_e = np.eye(nx)

    ocp.cost.yref = np.zeros((ny,))
    ocp.cost.yref_e = np.zeros((ny_e,))

    # set constraints
    Fmax = 80
    ocp.constraints.lbu = np.array([-Fmax])
    ocp.constraints.ubu = np.array([+Fmax])
    ocp.constraints.idxbu = np.array([0])

    ocp.constraints.x0 = np.array([0.0, np.pi, 0.0, 0.0])

    if constraint_variant == 'one_sided':
        ocp.constraints.lbx = np.array([-ACADOS_INFTY])
        ocp.constraints.ubx = np.array([+5.0])
        ocp.constraints.idxbx = np.array([0])
        expected_status = 0
    elif constraint_variant == 'one_sided_wrong_infty':
        ocp.constraints.lbx = np.array([-0.5*ACADOS_INFTY])
        ocp.constraints.ubx = np.array([+5.0])
        ocp.constraints.idxbx = np.array([0])
        if qp_solver == 'FULL_CONDENSING_DAQP':
            expected_status = 0
        elif qp_solver in ['FULL_CONDENSING_HPIPM', 'PARTIAL_CONDENSING_HPIPM']:
            # complementarity residual does not converge to tolerance if infty is not defined properly
            expected_status = 2

    # set options
    ocp.solver_options.qp_solver = qp_solver
    ocp.solver_options.hessian_approx = 'GAUSS_NEWTON'
    ocp.solver_options.integrator_type = 'ERK'
    ocp.solver_options.nlp_solver_type = 'SQP'
    ocp.solver_options.tol = 1e-7

    # set prediction horizon
    ocp.solver_options.tf = Tf

    # create solver
    ocp_solver = AcadosOcpSolver(ocp)

    # solve the problem defined here (original from code export), analog to 'minimal_example_ocp.py'
    simX0 = np.zeros((N_horizon + 1, nx))
    simU0 = np.zeros((N_horizon, nu))

    print(80*'-')
    print(f'solve original code with N = {N_horizon} and Tf = {Tf} s:')
    status = ocp_solver.solve()
    ocp_solver.print_statistics()

    if status != expected_status:
        raise Exception(f"expected status {expected_status}, got {status} for constraint_variant {constraint_variant}.")

    # get solution
    for i in range(N_horizon):
        simX0[i, :] = ocp_solver.get(i, "x")
        simU0[i, :] = ocp_solver.get(i, "u")
    simX0[N_horizon, :] = ocp_solver.get(N_horizon, "x")

    lambdas = [ocp_solver.get(i, "lam") for i in range(1, N_horizon)]
    for lam in lambdas:
        assert np.all(lam >= 0)

    # if unbounded constraint is defined properly, lambda should be zero
    i_infty = 1
    if constraint_variant == 'one_sided':
        for lam in lambdas:
            assert lam[i_infty] == 0
    elif (constraint_variant == 'one_sided_wrong_infty' and
        qp_solver in ['FULL_CONDENSING_HPIPM', 'PARTIAL_CONDENSING_HPIPM']):
        for lam in lambdas:
            assert lam[i_infty] != 0

    if constraint_variant == 'one_sided':
        stage = 1
        # check setting lambdas
        ocp_solver.set(stage, "lam", np.ones(lambdas[0].shape))
        lam = ocp_solver.get(stage, "lam")
        assert lam[i_infty] == 0

        # check updating bound
        i_infty_new = 2 # lam is ordered as lbu, lbx, ubu, ubx
        ocp_solver.constraints_set(stage, "lbx", -10.)
        ocp_solver.set(stage, "lam", np.ones(lambdas[0].shape))
        ocp_solver.constraints_set(stage, "ubu", ACADOS_INFTY)
        lam = ocp_solver.get(stage, "lam")

        assert lam[i_infty_new] == 0
        assert lam[i_infty] == 1.

    if PLOT:
        plot_pendulum(np.linspace(0, Tf, N_horizon + 1), Fmax, simU0, simX0, latexify=False, plt_show=True, X_true_label=f'N={N_horizon}, Tf={Tf}')

    ocp_solver = None

if __name__ == "__main__":
    for qp_solver in ['FULL_CONDENSING_HPIPM', 'PARTIAL_CONDENSING_HPIPM', 'FULL_CONDENSING_DAQP']:
        for constraint_variant in ['one_sided', 'one_sided_wrong_infty']:
            print(80*'-')
            print(f'constraint_variant = {constraint_variant}, qp_solver = {qp_solver}')
            main(constraint_variant=constraint_variant, qp_solver=qp_solver)

    # main(constraint_variant='one_sided', qp_solver='PARTIAL_CONDENSING_HPIPM')
