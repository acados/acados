#
# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# Licensed under the 2-Clause BSD License.

from acados_template import AcadosOcpSolver
import numpy as np
from sensitivity_utils import export_parametric_ocp, plot_cost_gradient_results


def main():
    """
    Evaluate policy and calculate its gradient for the pendulum on a cart with a parametric model.
    """

    p_nominal = 1.0
    x0 = np.array([0.0, np.pi / 2, 0.0, 0.0])
    delta_p = 0.002
    p_test = np.arange(p_nominal - 0.5, p_nominal + 0.5, delta_p)

    np_test = p_test.shape[0]
    N_horizon = 50
    T_horizon = 2.0
    Fmax = 80.0

    ocp = export_parametric_ocp(x0=x0, N_horizon=N_horizon, T_horizon=T_horizon, Fmax=Fmax, qp_solver_ric_alg=1)
    ocp.code_gen_options.with_value_sens_wrt_params = True
    acados_ocp_solver = AcadosOcpSolver(ocp)

    optimal_value_grad = np.zeros(np_test)
    optimal_value = np.zeros(np_test)

    pi = np.zeros(np_test)
    for i, p in enumerate(p_test):
        acados_ocp_solver.set_p_global_and_precompute_dependencies(np.array([p]))
        pi[i] = acados_ocp_solver.solve_for_x0(x0)[0]
        optimal_value[i] = acados_ocp_solver.get_cost()
        optimal_value_grad[i] = acados_ocp_solver.eval_and_get_optimal_value_gradient("p_global").item()

    # evaluate cost gradient
    optimal_value_grad_via_fd = np.gradient(optimal_value, delta_p)
    cost_reconstructed_np_grad = np.cumsum(optimal_value_grad_via_fd) * delta_p + optimal_value[0]

    plot_cost_gradient_results(p_test, optimal_value, optimal_value_grad, optimal_value_grad_via_fd, cost_reconstructed_np_grad)

    # checks
    test_tol = 1e-2
    mean_rel_diff = np.mean(np.abs(optimal_value_grad - optimal_value_grad_via_fd) / np.abs(optimal_value_grad_via_fd))
    print(f"Mean difference between value function gradient obtained by acados and via FD is {mean_rel_diff} should be < {test_tol}.")
    assert mean_rel_diff <= test_tol



if __name__ == "__main__":
    main()
