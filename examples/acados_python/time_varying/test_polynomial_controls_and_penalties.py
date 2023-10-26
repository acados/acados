import itertools
import numpy as np

from acados_template import casadi_length

from piecewiese_polynomial_control_example import create_ocp_solver

COST_TYPES = ['NONLINEAR_LS', 'CONVEX_OVER_NONLINEAR', 'NLS_TO_CONL']
EXPLICIT_SYMMETRIC_PENALTIES = [True, False]


def test_polynomial_controls_and_penalties():

    settings = [s for s in itertools.product(COST_TYPES, EXPLICIT_SYMMETRIC_PENALTIES)]

    N_horizon = 5
    degree_u_polynom = 2

    x_traj_list = []
    u_traj_list = []
    for s in settings:
        print(f"{s=}")
        (cost_type, explicit_symmetric_penalties) = s
        # create solver and extract
        ocp_solver, evaluate_polynomial_u_fun = create_ocp_solver(cost_type, N_horizon, degree_u_polynom, explicit_symmetric_penalties=explicit_symmetric_penalties)
        ocp = ocp_solver.acados_ocp
        model = ocp.model
        nu = casadi_length(model.u)
        nx = casadi_length(model.x)

        # solve OCP
        status = ocp_solver.solve()
        ocp_solver.print_statistics()

        if status != 0:
            raise Exception(f'acados returned status {status}.')

        cost_val = ocp_solver.get_cost()
        print(f"found optimal cost: {cost_val:.4e}")

        # get solution
        simX = np.ndarray((N_horizon+1, nx))
        simU = np.ndarray((N_horizon, nu))
        for i in range(N_horizon):
            simX[i,:] = ocp_solver.get(i, "x")
            simU[i,:] = ocp_solver.get(i, "u")
        simX[N_horizon,:] = ocp_solver.get(N_horizon, "x")

        x_traj_list.append(simX)
        u_traj_list.append(simU)
        del ocp_solver

    # check that the trajectories are the same
    tol_test = 1e-5
    i = 0
    max_diff = np.zeros(len(settings))
    for j in range(i+1, len(settings)):
        max_diff[i] = max(np.max(np.abs(u_traj_list[i] - u_traj_list[j])), np.max(np.abs(x_traj_list[i] - x_traj_list[j])))
        print(f"max_diff between {settings[i]} and {settings[j]}: {max_diff[i]}")

    for j in range(i+1, len(settings)):
        if not max_diff[i] < tol_test:
            raise Exception(f"Trajectories for settings {settings[i]} and {settings[j]} differ with {max_diff[i]} > {tol_test}.")

if __name__ == "__main__":
    test_polynomial_controls_and_penalties()