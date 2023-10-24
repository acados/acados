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

def solve_ocp(cost_discretization, cost_variant):

    # create ocp object to formulate the OCP
    ocp = AcadosOcp()

    # set model
    model = export_pendulum_ode_model()
    ocp.model = model

    nx = model.x.size()[0]
    nu = model.u.size()[0]
    ny = nx + nu
    ny_e = nx

    Tf = 1.0
    N = 20

    # set dimensions
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

    # set constraints
    Fmax = 80
    ocp.constraints.lbu = np.array([-Fmax])
    ocp.constraints.ubu = np.array([+Fmax])
    ocp.constraints.idxbu = np.array([0])

    ocp.constraints.x0 = np.array([0.0, np.pi, 0.0, 0.0])

    # set options
    ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM'
    ocp.solver_options.hessian_approx = 'GAUSS_NEWTON'
    ocp.solver_options.integrator_type = 'IRK'
    ocp.solver_options.collocation_type = 'EXPLICIT_RUNGE_KUTTA'
    ocp.solver_options.sim_method_num_stages = 1
    ocp.solver_options.sim_method_num_steps = 1
    ocp.solver_options.nlp_solver_type = 'SQP'
    ocp.solver_options.cost_discretization = cost_discretization

    # set prediction horizon
    ocp.solver_options.tf = Tf
    ocp_solver = AcadosOcpSolver(ocp, json_file='acados_ocp.json')

    # test setting HPIPM options
    ocp_solver.options_set('qp_tol_ineq', 1e-8)
    ocp_solver.options_set('qp_tau_min', 1e-10)
    ocp_solver.options_set('qp_mu0', 1e0)

    simX = np.ndarray((N + 1, nx))
    simU = np.ndarray((N, nu))

    print(80*'-')
    print(f'solve OCP with cost variant {cost_variant} discretization {cost_discretization} N = {N} and Tf = {Tf} s:')
    status = ocp_solver.solve()
    ocp_solver.print_statistics()

    if status != 0:
        raise Exception(f'acados returned status {status}.')

    # get solution
    for i in range(N):
        simX[i, :] = ocp_solver.get(i, "x")
        simU[i, :] = ocp_solver.get(i, "u")
    simX[N, :] = ocp_solver.get(N, "x")

    ocp_solver.store_iterate(filename=get_iterate_filename(cost_discretization, cost_variant), overwrite=True)

    if PLOT:# plot but don't halt
        plot_pendulum(np.linspace(0, Tf, N + 1), Fmax, simU, simX, latexify=False, plt_show=False, X_true_label=f'original: N={N}, Tf={Tf}')


def get_iterate_filename(cost_discretization, cost_variant):
    return f'final_iterate_{cost_discretization}_{cost_variant}.json'


def compare_iterates(cost_variant):
    import json
    ref_cost_discretization = COST_DISCRETIZATIONS[0]

    ref_iterate_filename = get_iterate_filename(ref_cost_discretization, cost_variant)
    with open(ref_iterate_filename, 'r') as f:
        ref_iterate = json.load(f)

    tol = 1e-10
    for cost_discretization in COST_DISCRETIZATIONS[1:]:
        iterate_filename = get_iterate_filename(cost_discretization, cost_variant)
        with open(iterate_filename, 'r') as f:
            iterate = json.load(f)

        assert iterate.keys() == ref_iterate.keys()

        errors = [np.max(np.abs((np.array(iterate[k]) - np.array(ref_iterate[k])))) for k in iterate]
        max_error = max(errors)
        print(f"max error {max_error:e}")
        if (max_error < tol):
            print(f"successfuly compared {len(COST_DISCRETIZATIONS)} cost discretizations for {cost_variant}")
        else:
            raise Exception(f"comparing {cost_variant=}, {cost_discretization=} failed with {max_error=}")


if __name__ == "__main__":
    for cost_variant in COST_VARIANTS:
        for cost_discretization in COST_DISCRETIZATIONS:
            solve_ocp(cost_discretization, cost_variant)
        compare_iterates(cost_variant)

