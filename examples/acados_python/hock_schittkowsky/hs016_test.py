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

from acados_template import AcadosOcp, AcadosOcpSolver, AcadosModel, ACADOS_INFTY, AcadosOcpFlattenedIterate
import numpy as np
from casadi import *


def solve_problem(qp_solver: str = 'FULL_CONDENSING_HPIPM', scale_qp_constraints: bool = False):
    print(f"Solving with {qp_solver} and scale_qp_constraints={scale_qp_constraints}")

    # create ocp object to formulate the OCP
    ocp = AcadosOcp()

    # set model
    model = AcadosModel()
    x = SX.sym('x', 2)

    # dynamics: identity
    model.disc_dyn_expr = x
    model.x = x
    model.name = f'hs_016'
    ocp.model = model

    # cost
    ocp.cost.cost_type_e = 'EXTERNAL'
    ocp.model.cost_expr_ext_cost_e = 100*(x[1] - x[0]**2)**2 + (1 - x[0])**2

    # constraints
    ocp.model.con_h_expr_e = vertcat(x[0]**2 + x[1], x[0] + x[1]**2)
    ocp.constraints.lh_e = np.array([0.0, 0.0])
    ocp.constraints.uh_e = np.array([ACADOS_INFTY, ACADOS_INFTY])

    # add bounds on x;
    ocp.constraints.idxbx_e = np.arange(2)
    ocp.constraints.ubx_e = np.array([0.5, 1.0])
    ocp.constraints.lbx_e = np.array([-0.5, -ACADOS_INFTY])

    # set options
    ocp.solver_options.N_horizon = 0
    ocp.solver_options.qp_solver = qp_solver
    ocp.solver_options.qp_solver_mu0 = 1e3
    ocp.solver_options.hessian_approx = 'EXACT'
    ocp.solver_options.regularize_method = 'MIRROR'
    ocp.solver_options.print_level = 1
    ocp.solver_options.nlp_solver_max_iter = 1000
    ocp.solver_options.qp_solver_iter_max = 1000

    # Search direction
    ocp.solver_options.nlp_solver_type = 'SQP_WITH_FEASIBLE_QP'

    # Globalization
    ocp.solver_options.globalization = 'FUNNEL_L1PEN_LINESEARCH'
    ocp.solver_options.globalization_full_step_dual = True
    ocp.solver_options.globalization_funnel_use_merit_fun_only = False

    # Scaling
    if scale_qp_constraints:
        ocp.solver_options.qpscaling_scale_objective = 'OBJECTIVE_GERSHGORIN'
        ocp.solver_options.qpscaling_scale_constraints = 'INF_NORM'
    ocp.code_export_directory = f'c_generated_code_{model.name}'

    ocp_solver = AcadosOcpSolver(ocp, json_file=f'{model.name}.json', verbose=False)

    # initialize solver
    xinit = np.array([-2, 1])
    ocp_solver.set(0, "x", xinit)

    # solve
    status = ocp_solver.solve()

    # get solution
    assert status == 0, f"Solver failed with status {status}"

    sol = ocp_solver.store_iterate_to_flat_obj()

    return sol


def main():
    sol_list = []
    for qp_solver in ['FULL_CONDENSING_HPIPM', 'PARTIAL_CONDENSING_HPIPM']:
        for scaling in [False, True]:
            sol = solve_problem(qp_solver, scaling)
            sol_list.append(sol)

    ref_sol = sol_list[0]
    for i, sol in enumerate(sol_list[1:]):
        if not AcadosOcpFlattenedIterate.allclose(ref_sol, sol):
            raise ValueError(f"Solution does not match reference close enough!")
        else:
            print(f"Solution {i+1} matches reference solution.")

if __name__ == '__main__':
    main()
