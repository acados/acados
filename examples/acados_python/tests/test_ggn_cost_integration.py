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

COST_TYPE = ['NONLINEAR_LS', 'CONVEX_OVER_NONLINEAR']
PLOT = False
COST_DISCRETIZATIONS = ['EULER', 'INTEGRATOR']

TOL = 1e-10


def solve_ocp(cost_discretization, cost_type, num_stages, collocation_type):

    model = export_pendulum_ode_model()

    ocp = AcadosOcp()
    ocp.model = model

    nx = model.x.size()[0]
    nu = model.u.size()[0]
    ny = nx + nu

    Tf = 1.0
    N = 20
    ocp.dims.N = N

    Q = 2 * np.diag([1e3, 1e3, 1e-2, 1e-2, 0.1, 1e-3])
    R = 2 * np.diag([1e-2])
    cost_W = scipy.linalg.block_diag(Q, R)

    ocp.cost.cost_type = cost_type
    ocp.cost.cost_type_e = cost_type

    ny = nx + nu + 2
    ny_e = nx + 2
    ocp.model.cost_y_expr = ca.vertcat(model.x, model.x[-1]**2, model.x[-1]*model.u, model.u)
    ocp.model.cost_y_expr_e = ca.vertcat(model.x, model.x[-1]**2, model.x[-1])
    ocp.cost.yref = np.zeros((ny, ))
    ocp.cost.yref_e = np.zeros((ny_e, ))

    if cost_type == "NONLINEAR_LS":
        ocp.cost.W = cost_W
        ocp.cost.W_e = Q

    elif cost_type == 'CONVEX_OVER_NONLINEAR':
        r = ca.SX.sym('r', ny)
        r_e = ca.SX.sym('r_e', ny_e)

        ocp.model.cost_psi_expr = 0.5 * (r.T @ cost_W @ r)
        ocp.model.cost_psi_expr_e = 0.5 * (r_e.T @ Q @ r_e)

        ocp.model.cost_r_in_psi_expr = r
        ocp.model.cost_r_in_psi_expr_e = r_e

    else:
        raise Exception(f"cost_type {cost_type} not supported")


    # augment with cost state
    cost_state = ca.SX.sym('cost_state')
    cost_state_dot = ca.SX.sym('cost_state_dot')
    res = ocp.model.cost_y_expr - ocp.cost.yref
    cost = 0.5*res.T @ cost_W @ res

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
    ocp.solver_options.cost_discretization = cost_discretization
    ocp.solver_options.nlp_solver_max_iter = 100

    # for debugging:
    # ocp.solver_options.nlp_solver_max_iter = 1
    ocp.solver_options.collocation_type = collocation_type
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
    print(f'solve OCP with {cost_type} {cost_discretization} N = {N} and Tf = {Tf} s:')
    status = ocp_solver.solve()
    # ocp_solver.dump_last_qp_to_json(f'qp_{cost_discretization}.json', overwrite=True)

    ocp_solver.print_statistics()

    if status != 0:
        raise Exception(f'acados returned status {status}.')

    ocp_solver.store_iterate(filename=get_iterate_filename(cost_discretization, cost_type), overwrite=True)

    # get solution
    for i in range(N):
        simX[i, :] = ocp_solver.get(i, "x")
        simU[i, :] = ocp_solver.get(i, "u")
    simX[N, :] = ocp_solver.get(N, "x")

    # compare cost and value of cost state
    cost_solver = ocp_solver.get_cost()

    xN = simX[N, :nx]

    resN = ca.vertcat(xN, xN[-1]**2, xN[-1]).full()
    terminal_cost = 0.5* resN.T @ Q @ resN
    cost_state = simX[-1, -1] + terminal_cost

    abs_diff = np.abs(cost_solver - cost_state).item()

    print(f"\nComparing solver cost and cost state for {cost_type=}, {num_stages=}:\n  {abs_diff=:.3e}")
    if abs_diff < TOL:
        print('  SUCCESS!\n')
    else:
        raise Exception(f"  ERROR for {cost_type=}, {num_stages=}:\n  {abs_diff=:.3e}\n")

    if PLOT:# plot but don't halt
        plot_pendulum(np.linspace(0, Tf, N + 1), Fmax, simU, simX[:, :-1], latexify=False, plt_show=True, X_true_label=f'original: N={N}, Tf={Tf}')


def get_iterate_filename(cost_discretization, cost_type):
    return f'final_iterate_{cost_discretization}_{cost_type}.json'

def compare_iterates(cost_type):
    import json
    ref_cost_discretization = COST_DISCRETIZATIONS[0]

    ref_iterate_filename = get_iterate_filename(ref_cost_discretization, cost_type)
    with open(ref_iterate_filename, 'r') as f:
        ref_iterate = json.load(f)

    tol = 1e-10
    for cost_discretization in COST_DISCRETIZATIONS[1:]:
        iterate_filename = get_iterate_filename(cost_discretization, cost_type)
        with open(iterate_filename, 'r') as f:
            iterate = json.load(f)

        assert iterate.keys() == ref_iterate.keys()

        errors = [np.max(np.abs((np.array(iterate[k]) - np.array(ref_iterate[k])))) for k in iterate]
        max_error = max(errors)
        print(f"max error {max_error:e}")
        if (max_error < tol):
            print(f"successfuly compared {len(COST_DISCRETIZATIONS)} cost discretizations for {cost_type}")
        else:
            raise Exception(f"comparing {cost_type=}, {cost_discretization=} failed with {max_error=}")

if __name__ == "__main__":

    for cost_type in COST_TYPE:
        for cost_discretization in COST_DISCRETIZATIONS:
            solve_ocp(cost_discretization, cost_type, num_stages=1, collocation_type='EXPLICIT_RUNGE_KUTTA')
        compare_iterates(cost_type)

    for cost_type in COST_TYPE:
            solve_ocp('INTEGRATOR', cost_type, num_stages=3, collocation_type='GAUSS_LEGENDRE')
