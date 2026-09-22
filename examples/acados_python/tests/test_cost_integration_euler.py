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

COST_VARIANTS = ['PARTIAL_STATE_PENALTY', 'FULL_STATE_PENALTY', 'DOUBLE_STATE_PENALTY', 'CREATIVE_NONLINEAR']
PLOT = False
COST_DISCRETIZATIONS = ['EULER', 'INTEGRATOR']

T_HORIZON = 1.0
N_HORIZON = 20
F_MAX = 80

def formulate_ocp(cost_variant):
    ocp = AcadosOcp()
    model = export_pendulum_ode_model()
    ocp.model = model

    nx = model.x.rows()
    nu = model.u.rows()
    ny = nx + nu
    ny_e = nx

    # set dimensions
    ocp.solver_options.N_horizon = N_HORIZON

    # set cost
    Q = 2 * np.diag([1e3, 1e3, 1e-2, 1e-2])
    Q[0, 1] = 10.
    Q[1, 0] = 10.
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

    # set constraints
    ocp.constraints.lbu = np.array([-F_MAX])
    ocp.constraints.ubu = np.array([+F_MAX])
    ocp.constraints.idxbu = np.array([0])

    ocp.constraints.x0 = np.array([0.0, np.pi, 0.0, 0.0])

    return ocp


def set_options(ocp, cost_discretization):
    ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM'
    ocp.solver_options.hessian_approx = 'GAUSS_NEWTON'
    ocp.solver_options.integrator_type = 'IRK'
    ocp.solver_options.collocation_type = 'EXPLICIT_RUNGE_KUTTA'
    ocp.solver_options.sim_method_num_stages = 1
    ocp.solver_options.sim_method_num_steps = 1
    ocp.solver_options.nlp_solver_type = 'SQP'
    ocp.solver_options.cost_discretization = cost_discretization
    ocp.solver_options.tf = 1.0


def solve_ocp(cost_discretization, cost_variant):
    ocp = formulate_ocp(cost_variant)
    set_options(ocp, cost_discretization)
    ocp_solver = AcadosOcpSolver(ocp)

    # test setting HPIPM options
    ocp_solver.options_set('qp_tol_ineq', 1e-8)
    ocp_solver.options_set('qp_tau_min', 1e-10)
    ocp_solver.options_set('qp_mu0', 1e0)

    print(80*'-')
    print(f'solve OCP with cost variant {cost_variant} discretization {cost_discretization} N_HORIZON = {N_HORIZON} and T_HORIZON = {T_HORIZON} s:')
    status = ocp_solver.solve()
    ocp_solver.print_statistics()

    if status != 0:
        raise Exception(f'acados returned status {status}.')

    iterate = ocp_solver.get_iterate()
    simX = np.array(iterate.x)
    simU = np.array(iterate.u)

    if PLOT:
        plot_pendulum(np.linspace(0, T_HORIZON, N_HORIZON + 1), F_MAX, simU, simX, latexify=False, plt_show=False, X_true_label=f'original: N_HORIZON={N_HORIZON}, T_HORIZON={T_HORIZON}')

    return iterate


def compare_iterates(cost_variant, reference_iterate, iterate):
    if not reference_iterate.allclose(iterate, atol=1e-10, rtol=0.0):
        raise Exception(f"comparing {cost_variant=} failed with mismatching iterates")

    print(f"successfuly compared {len(COST_DISCRETIZATIONS)} cost discretizations for {cost_variant}")


if __name__ == "__main__":
    for cost_variant in COST_VARIANTS:
        reference_iterate = None
        for cost_discretization in COST_DISCRETIZATIONS:
            iterate = solve_ocp(cost_discretization, cost_variant)
            if reference_iterate is None:
                reference_iterate = iterate
            else:
                compare_iterates(cost_variant, reference_iterate, iterate)

