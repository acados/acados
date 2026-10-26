import numpy as np
import casadi as ca

from acados_template import AcadosOcp, AcadosModel, AcadosCasadiOcpSolver

N_HORIZON = 2

def formulate_ocp(u_ref: float, lb: float, ub: float, constraint_type: str = 'bu') -> AcadosOcp:
    """
    L(u, lambda_u, lambda_l) = 0.5(u - u_ref)^2 + lambda_u (u - 1) + lambda_l (-u - 1)
    The optimal control is u = clip(u_ref, lb, ub)
    """
    ocp = AcadosOcp()

    model = AcadosModel()
    model.name = 'strict_complementarity_toy'
    model.x = ca.SX.sym('x', 1)
    model.u = ca.SX.sym('u', 1)
    model.disc_dyn_expr = model.x + model.u
    ocp.model = model

    ocp.solver_options.N_horizon = N_HORIZON
    ocp.solver_options.tf = float(N_HORIZON)
    ocp.solver_options.integrator_type = 'DISCRETE'
    ocp.solver_options.hessian_approx = 'EXACT'

    ocp.cost.cost_type = 'EXTERNAL'
    ocp.cost.cost_type_e = 'EXTERNAL'
    model.cost_expr_ext_cost = 0.5 * (model.u - u_ref)**2
    model.cost_expr_ext_cost_e = 0 * model.x

    ocp.constraints.x0 = np.array([0.0])
    if constraint_type == 'bu':
        ocp.constraints.idxbu = np.array([0])
        ocp.constraints.lbu = np.array([lb])
        ocp.constraints.ubu = np.array([ub])
    elif constraint_type == 'g':
        ocp.constraints.C = np.zeros((1, 1))
        ocp.constraints.D = np.ones((1, 1))
        ocp.constraints.lg = np.array([lb])
        ocp.constraints.ug = np.array([ub])
    else:
        raise ValueError(f'Unknown constraint_type {constraint_type}.')

    return ocp

def solve_ocp(u_ref: float, lb: float, ub: float, constraint_type: str) -> AcadosCasadiOcpSolver:
    ocp = formulate_ocp(u_ref, lb, ub, constraint_type)
    solver = AcadosCasadiOcpSolver(ocp, solver="ipopt", verbose=False,
                                   casadi_solver_opts={'ipopt.print_level': 0, 'ipopt.sb': 'yes'})
    status = solver.solve()
    if status != 0:
        raise Exception(f'casadi solver returned status {status}.')
    return solver

def main(test_case):
    description, u_ref, expected = test_case

    for constraint_type in ['bu', 'g']:
        solver = solve_ocp(u_ref, -1.0, 1.0, constraint_type)
        result = solver.satisfies_strict_complementarity()
        result_stages = solver.satisfies_strict_complementarity_stages()
        expected_stages = N_HORIZON * [expected] + [True]
        print(f"{constraint_type}, {description}: strict complementarity = {result}, stage-wise = {result_stages}")

        if result != expected:
            raise Exception(f"{constraint_type}, {description}: expected strict complementarity {expected}, got {result}.")
        if result_stages != expected_stages:
            raise Exception(f"{constraint_type}, {description}: expected stage-wise result {expected_stages}, got {result_stages}.")

if __name__ == "__main__":
    test_cases = [
        ("upper bound strictly active", 2.0, True),
        ("lower bound strictly active", -2.0, True),
        ("inactive", 0.0, True),
        ("upper bound weakly active", 1.0, False),
        ("lower bound weakly active", -1.0, False),
    ]
    for test_case in test_cases:
        main(test_case)
