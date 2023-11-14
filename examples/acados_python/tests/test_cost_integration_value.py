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

from acados_template import AcadosOcp, AcadosOcpSolver
from pendulum_model import export_pendulum_ode_model
import numpy as np
import casadi as ca
import scipy.linalg
from utils import plot_pendulum

COST_VARIANTS = ['FULL_STATE_PENALTY', 'PARTIAL_STATE_PENALTY', 'DOUBLE_STATE_PENALTY', 'CREATIVE_NONLINEAR']
PLOT = False
NUM_STAGES = [1, 3]

TOL = 1e-10


def solve_ocp(cost_variant, num_stages):

    model = export_pendulum_ode_model()

    ocp = AcadosOcp()
    ocp.model = model

    nx = model.x.size()[0]
    nu = model.u.size()[0]
    ny = nx + nu
    ny_e = nx

    Tf = 1.0
    N = 20
    ocp.dims.N = N

    # set cost
    Q = 2 * np.diag([1e3, 1e3, 1e-2, 1e-2])
    R = 2 * np.diag([1e-2])

    ocp.cost.cost_type = 'NONLINEAR_LS'

    if cost_variant == "FULL_STATE_PENALTY":
        ny = nx + nu
        ocp.model.cost_y_expr = ca.vertcat(model.x, model.u)
        ocp.cost.W = scipy.linalg.block_diag(Q, R)
        ocp.cost.yref = np.zeros((ny, ))

    elif cost_variant == 'PARTIAL_STATE_PENALTY':
        nyx = 2
        ny = nyx + nu
        ocp.model.cost_y_expr = ca.vertcat(model.x[:nyx], model.u)
        ocp.cost.W = scipy.linalg.block_diag(Q[:nyx, :nyx], R)
        ocp.cost.yref = np.zeros((ny, ))

    elif cost_variant == 'DOUBLE_STATE_PENALTY':
        ny = 2*nx + nu
        ocp.model.cost_y_expr = ca.vertcat(model.x, model.x, model.u)
        ocp.cost.W = scipy.linalg.block_diag(Q, Q, R)
        ocp.cost.yref = np.zeros((ny, ))

    elif cost_variant == 'CREATIVE_NONLINEAR':
        ocp.model.cost_y_expr = ca.vertcat(model.x[2], model.u, 0.1*(model.x[0]+model.u[0]+1.)**3)
        ny = max(ocp.model.cost_y_expr.shape)
        ocp.cost.W = Q[:ny, :ny]
        ocp.cost.yref = np.zeros((ny, ))
    else:
        raise Exception(f"cost_variant {cost_variant} not supported")

    ny_e = nx
    ocp.cost.cost_type_e = 'NONLINEAR_LS'
    ocp.model.cost_y_expr_e = model.x
    ocp.cost.W_e = Q
    ocp.cost.yref_e = np.zeros((ny_e, ))

    # augment with cost state
    cost_state = ca.SX.sym('cost_state')
    cost_state_dot = ca.SX.sym('cost_state_dot')
    res = ocp.model.cost_y_expr - ocp.cost.yref
    cost = 0.5*res.T @ ocp.cost.W @ res

    ocp.model.f_expl_expr = ca.vertcat(ocp.model.f_expl_expr, cost)
    ocp.model.x = ca.vertcat(ocp.model.x, cost_state)
    ocp.model.xdot = ca.vertcat(ocp.model.xdot, cost_state_dot)
    ocp.model.f_impl_expr = ocp.model.f_expl_expr - ocp.model.xdot

    # set constraints
    Fmax = 80
    ocp.constraints.lbu = np.array([-Fmax])
    ocp.constraints.ubu = np.array([+Fmax])
    ocp.constraints.idxbu = np.array([0])
    ocp.constraints.x0 = np.array([0.0, np.pi, 0.0, 0.0, 0.0])

    # set options
    ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM'  # FULL_CONDENSING_QPOASES
    ocp.solver_options.hessian_approx = 'GAUSS_NEWTON'
    ocp.solver_options.integrator_type = 'IRK'
    ocp.solver_options.sim_method_num_stages = num_stages
    ocp.solver_options.sim_method_num_steps = 1
    ocp.solver_options.nlp_solver_type = 'SQP'  # SQP_RTI, SQP
    ocp.solver_options.cost_discretization = 'INTEGRATOR'

    # set prediction horizon
    ocp.solver_options.tf = Tf
    ocp_solver = AcadosOcpSolver(ocp, json_file='acados_ocp.json')

    # test setting HPIPM options
    ocp_solver.options_set('qp_tol_ineq', 1e-8)
    ocp_solver.options_set('qp_tau_min', 1e-10)
    ocp_solver.options_set('qp_mu0', 1e0)

    simX = np.ndarray((N + 1, nx+1))
    simU = np.ndarray((N, nu))

    print(80*'-')
    print(f'solve original code with N = {N} and Tf = {Tf} s:')
    status = ocp_solver.solve()

    ocp_solver.print_statistics()

    if status != 0:
        raise Exception(f'acados returned status {status}.')

    # get solution
    for i in range(N):
        simX[i, :] = ocp_solver.get(i, "x")
        simU[i, :] = ocp_solver.get(i, "u")
    simX[N, :] = ocp_solver.get(N, "x")

    # compare cost and value of cost state
    cost_solver = ocp_solver.get_cost()

    xN = simX[N, :nx]
    terminal_cost = 0.5*xN @ ocp.cost.W_e @ xN
    cost_state = simX[-1, -1] + terminal_cost

    abs_diff = np.abs(cost_solver - cost_state)

    print(f"\nComparing solver cost and cost state for {cost_variant=}, {num_stages=}:\n  {abs_diff=:.3e}")
    if abs_diff < TOL:
        print('  SUCCESS!\n')
    else:
        raise Exception(f"  ERROR for {cost_variant=}, {num_stages=}:\n  {abs_diff=:.3e}\n")

    if PLOT:# plot but don't halt
        plot_pendulum(np.linspace(0, Tf, N + 1), Fmax, simU, simX, latexify=False, plt_show=False, X_true_label=f'original: N={N}, Tf={Tf}')


if __name__ == "__main__":
    for num_stages in NUM_STAGES:
        for cost_variant in COST_VARIANTS:
            solve_ocp(cost_variant, num_stages)



