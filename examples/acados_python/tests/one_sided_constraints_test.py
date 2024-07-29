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

def main(constraint_variant='one_sided'):

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
    ocp.dims.N = N_horizon

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
        ocp.constraints.lbx = np.array([-2*ACADOS_INFTY])
        ocp.constraints.ubx = np.array([+5.0])
        ocp.constraints.idxbx = np.array([0])
        # complementarity residual does not converge to tolerance if infty is too large
        expected_status = 2

    # set options
    ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM'
    ocp.solver_options.hessian_approx = 'GAUSS_NEWTON'
    ocp.solver_options.integrator_type = 'ERK'
    ocp.solver_options.nlp_solver_type = 'SQP'

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


    if PLOT:
        plot_pendulum(np.linspace(0, Tf, N_horizon + 1), Fmax, simU0, simX0, latexify=False, plt_show=True, X_true_label=f'N={N_horizon}, Tf={Tf}')

    ocp_solver = None

if __name__ == "__main__":
    main(constraint_variant='one_sided')
    main(constraint_variant='one_sided_wrong_infty')
