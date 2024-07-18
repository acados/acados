from acados_template import AcadosOcp, AcadosOcpSolver, AcadosModel
import numpy as np
from casadi import *
import casadi as cs
from itertools import product
# Simplest NLP with Marathos effect
#
# min log(exp(x) + exp(-x))
#

# Settings
TOL = 1e-6

def main():
    # run test cases
    params = {'globalization': ['FIXED_STEP', 'FUNNEL_L1PEN_LINESEARCH', 'MERIT_BACKTRACKING']}

    keys, values = zip(*params.items())
    for combination in product(*values):
        setting = dict(zip(keys, combination))
        solve_problem_with_setting(setting)


def solve_problem_with_setting(setting):

    globalization = setting['globalization']

    # create ocp object to formulate the OCP
    ocp = AcadosOcp()

    # set model
    model = AcadosModel()
    x = SX.sym('x')

    # dynamics: identity
    model.disc_dyn_expr = x
    model.x = x
    model.u = SX.sym('u', 0, 0)
    model.p = []
    model.name = f'convex_globalization_problem'
    ocp.model = model

    # discretization
    Tf = 1
    N = 1
    ocp.dims.N = N
    ocp.solver_options.tf = Tf

    # cost
    ocp.cost.cost_type_e = 'EXTERNAL'
    ocp.model.cost_expr_ext_cost_e = cs.log(cs.exp(x) + cs.exp(-x))

    # set options
    ocp.solver_options.qp_solver = 'FULL_CONDENSING_HPIPM'
    ocp.solver_options.hessian_approx = 'EXACT'
    ocp.solver_options.integrator_type = 'DISCRETE'
    ocp.solver_options.print_level = 1
    ocp.solver_options.tol = TOL
    ocp.solver_options.nlp_solver_type = 'SQP'
    ocp.solver_options.globalization = globalization
    ocp.solver_options.qp_solver_iter_max = 400
    ocp.solver_options.regularize_method = 'MIRROR'
    ocp.solver_options.qp_tol = 5e-7

    ocp.solver_options.nlp_solver_max_iter = 100
    ocp_solver = AcadosOcpSolver(ocp, json_file=f'{model.name}.json')

    # initialize solver
    xinit = np.array([1.5])
    [ocp_solver.set(i, "x", xinit) for i in range(N+1)]

    print(f"solved convex globalization test problem with settings {setting}")
    # solve
    ocp_solver.solve()

    status = ocp_solver.get_status()

    # get solution
    solution = ocp_solver.get(0, "x")

    # compare to analytical solution
    exact_solution = np.array([0])
    sol_err = max(np.abs(solution - exact_solution ))

    if globalization in ['FUNNEL_L1PEN_LINESEARCH', 'MERIT_BACKTRACKING']:
        assert status == 0, f'{globalization} did not converge. Algorithm should converge!'
        assert sol_err <= TOL*1e1, f"numerical solutions do not match to analytical solution with tolerance {TOL}"
    else:
        assert status != 0, 'Fixed step converged. Theoretically impossible!'

if __name__ == '__main__':
    main()
