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

from acados_template import AcadosOcp, AcadosOcpSolver, AcadosModel, ACADOS_INFTY
import numpy as np
from casadi import *
from matplotlib import pyplot as plt
from itertools import product
# Simplest NLP with Maratos effect
#
# Problem is in one dimension
#
# min -x
#
# s.t.  x <= 1
#       x^2 >= 4

# Settings
PLOT = False
FOR_LOOPING = False # call solver in for loop to get all iterates
TOL = 1e-6

def main():
    # run test cases
    params = {'globalization': ['FUNNEL_L1PEN_LINESEARCH']}

    keys, values = zip(*params.items())
    for combination in product(*values):
        setting = dict(zip(keys, combination))
        solve_infeasible_linearization(setting)


def solve_infeasible_linearization(setting):

    globalization = setting['globalization']

    # create ocp object to formulate the OCP
    ocp = AcadosOcp()

    # set model
    model = AcadosModel()
    x = SX.sym('x')

    # dynamics: identity
    model.disc_dyn_expr = x
    model.x = x
    model.u = SX.sym('u', 0, 0) # [] / None doesnt work
    model.p = []
    model.name = f'inconsistent_qp_linearization'
    ocp.model = model

    # discretization
    Tf = 1
    N = 1
    ocp.solver_options.N_horizon = N
    ocp.solver_options.tf = Tf

    # cost
    ocp.cost.cost_type_e = 'EXTERNAL'
    ocp.model.cost_expr_ext_cost_e = -x

    # constraints
    ocp.model.con_h_expr_0 = x**2
    ocp.constraints.lh_0 = np.array([4.0])
    ocp.constraints.uh_0 = np.array([ACADOS_INFTY])

    # add bounds on x
    nx = 1
    ocp.constraints.idxbx_0 = np.array(range(nx))
    ocp.constraints.lbx_0 = -ACADOS_INFTY * np.ones((nx))
    ocp.constraints.ubx_0 = 1 * np.ones((nx))

    # set options
    ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM'
    ocp.solver_options.qp_solver_cond_N = N
    ocp.solver_options.hessian_approx = 'EXACT'
    ocp.solver_options.integrator_type = 'DISCRETE'
    if globalization == 'FUNNEL_L1PEN_LINESEARCH':
        ocp.solver_options.print_level = 1
    ocp.solver_options.tol = TOL
    ocp.solver_options.nlp_solver_type = 'SQP'
    ocp.solver_options.globalization = globalization
    ocp.solver_options.alpha_min = 1e-15
    SQP_max_iter = 300
    ocp.solver_options.qp_solver_iter_max = 400
    ocp.solver_options.regularize_method = 'MIRROR'
    ocp.solver_options.qp_tol = 5e-7

    ocp.solver_options.nlp_solver_max_iter = SQP_max_iter
    ocp_solver = AcadosOcpSolver(ocp, json_file=f'{model.name}.json')

    # initialize solver
    xinit = np.array([1.0])
    [ocp_solver.set(i, "x", xinit) for i in range(N+1)]

    # solve
    status = ocp_solver.solve()
    if status != 0:
        raise RuntimeError("Solve failed, since QP infeasible!")

    # get solution
    solution = ocp_solver.get(0, "x")

    # compare to analytical solution
    exact_solution = np.array([-2.0])
    sol_err = max(np.abs(solution - exact_solution ))


if __name__ == '__main__':
    main()
