# -*- coding: future_fstrings -*-
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
sys.path.insert(0, '../pendulum_on_cart/solution_sensitivities')
import numpy as np
from sensitivity_utils import plot_cost_gradient_results
from acados_template import AcadosOcpSolver
from setup_parametric_ocp import PARAM_VALUE_DICT, export_parametric_ocp


def main():

    learnable_params = ["A", "Q", "b"]
    x0 = np.array([0.1, -0.2])

    delta_p = 0.005
    p_nominal = np.concatenate([PARAM_VALUE_DICT[k].flatten() for k in learnable_params]).flatten()

    p_test = np.arange(0.0, 1.0, delta_p)
    np_test = p_test.shape[0]

    dp = np.ones(p_nominal.shape).reshape((-1,))
    # Q
    dp[4] = 10.0
    dp[5] = 0.0
    dp[6] = 0.0
    dp[7] = 10.0
    # b
    dp[8] = -0.2
    dp[9] = -0.2

    ocp = export_parametric_ocp(PARAM_VALUE_DICT, learnable_params = learnable_params)
    ocp.solver_options.with_value_sens_wrt_params = True
    acados_ocp_solver = AcadosOcpSolver(ocp)

    optimal_value_grad = np.zeros((np_test,))
    optimal_value = np.zeros((np_test,))

    pi = np.zeros(np_test)
    for i, p in enumerate(p_test):
        p_val = p_nominal + p * dp
        acados_ocp_solver.set_p_global_and_precompute_dependencies(p_val)
        pi[i] = acados_ocp_solver.solve_for_x0(x0, fail_on_nonzero_status=False)[0]

        status = acados_ocp_solver.get_status()
        if status != 0:
            print(f"Solver failed with status {status} for p_val = {p_val}.")

        optimal_value[i] = acados_ocp_solver.get_cost()
        optimal_value_grad[i] = acados_ocp_solver.eval_and_get_optimal_value_gradient("p_global") @ dp

    # evaluate cost gradient
    optimal_value_grad_via_fd = np.gradient(optimal_value, delta_p)
    cost_reconstructed_np_grad = np.cumsum(optimal_value_grad_via_fd) * delta_p + optimal_value[0]
    cost_reconstructed_acados = np.cumsum(optimal_value_grad) * delta_p + optimal_value[0]

    plot_cost_gradient_results(p_test, optimal_value, optimal_value_grad,
                               optimal_value_grad_via_fd, cost_reconstructed_np_grad,
                               cost_reconstructed_acados, y_scale_log=True,
                               title=f"varying parameters {', '.join(learnable_params)} in direction {dp}",
                               xlabel=r"$\alpha$ in $p+\alpha \Delta p$")

    # checks
    test_tol = 1e-3
    median_diff = np.median(np.abs(optimal_value_grad - optimal_value_grad_via_fd))
    print(f"Median difference between value function gradient obtained by acados and via FD is {median_diff:.2e} should be < {test_tol:.2e}.")
    assert median_diff <= test_tol

if __name__ == "__main__":
    main()
