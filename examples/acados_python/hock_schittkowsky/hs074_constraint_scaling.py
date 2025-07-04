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
import casadi as ca


def solve_problem_with_constraint_scaling(scale_constraints):

    # create ocp object to formulate the OCP
    ocp = AcadosOcp()

    # set model
    model = AcadosModel()
    x = ca.SX.sym('x', 4)

    # dynamics: identity
    model.disc_dyn_expr = x
    model.x = x
    model.name = f'hs_074'
    ocp.model = model

    a = 0.55

    # cost
    ocp.cost.cost_type_e = 'EXTERNAL'
    ocp.model.cost_expr_ext_cost_e = 3*x[0] + 1.0e-6*x[0]**3 + 2*x[1] + 2.0e-6*x[1]**3/3

    # constraints
    g = ca.SX.zeros(4, 1)
    g[0] =  x[3] - x[2]
    g[1] = x[0]  - 1000*ca.sin(-x[2] - 0.25) - 1000*ca.sin(-x[3] - 0.25)
    g[2] = x[1]  - 1000*ca.sin(x[2] - 0.25) - 1000*ca.sin(x[2]-x[3] - 0.25)
    g[3] = 1000*ca.sin(x[3] - 0.25) + 1000*ca.sin(x[3] - x[2] - 0.25)

    ocp.model.con_h_expr_e = g
    ocp.constraints.lh_e = np.array([-a, 894.8, 894.8, -1294.8])
    ocp.constraints.uh_e = np.array([a, 894.8, 894.8, -1294.8])

    # add bounds on x;
    ocp.constraints.idxbx_e = np.arange(4)
    ocp.constraints.ubx_e = np.array([1200, 1200, a, a])
    ocp.constraints.lbx_e = np.array([0.0, 0.0, -a, -a])

    # set options
    ocp.solver_options.N_horizon = 0
    ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM'
    ocp.solver_options.qp_solver_mu0 = 1e3
    ocp.solver_options.hessian_approx = 'EXACT'
    ocp.solver_options.regularize_method = 'MIRROR'
    ocp.solver_options.print_level = 1
    ocp.solver_options.nlp_solver_max_iter = 1000
    ocp.solver_options.qp_solver_iter_max = 1000
    ocp.solver_options.nlp_solver_type = 'SQP_WITH_FEASIBLE_QP'

    # Globalization
    ocp.solver_options.globalization = 'FUNNEL_L1PEN_LINESEARCH'
    ocp.solver_options.globalization_full_step_dual = True
    ocp.solver_options.globalization_funnel_use_merit_fun_only = False

    # Scaling
    ocp.solver_options.qpscaling_scale_objective = 'OBJECTIVE_GERSHGORIN'
    if scale_constraints:
        ocp.solver_options.qpscaling_scale_constraints = 'INF_NORM'

    ocp.code_export_directory = f'c_generated_code_{model.name}'
    ocp_solver = AcadosOcpSolver(ocp, json_file=f'{model.name}.json', verbose = False)

    # initialize solver
    xinit = np.zeros(4)
    ocp_solver.set(0, "x", xinit)

    # solve
    status = ocp_solver.solve()

    # checks
    obj_scale = ocp_solver.get_qp_scaling_objective()
    print(f"Objective scaling: {obj_scale:.4e}")
    if scale_constraints:
        constr_scale = ocp_solver.get_qp_scaling_constraints(stage=0)
        print(f"Constraints scaling factors: {constr_scale}")

    if scale_constraints:
        assert status == 0, "Scaling of the constraints was not succesful!"
    else:
        assert status == 4, "Problem should not be solvable without scaling!"
    del ocp_solver

def main():
    # run test cases
    print("\nTest standard unscaled version, HPIPM should fail:")
    solve_problem_with_constraint_scaling(scale_constraints=False)
    print("\n\n----------------------------------------------")
    print("\nTest constraint scaling version, HPIPM should fail:")
    solve_problem_with_constraint_scaling(scale_constraints=True)
    print("\n\n----------------------------------------------")

if __name__ == '__main__':
    main()
