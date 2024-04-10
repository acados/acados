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
    ocp.solver_options.with_value_sens_wrt_params = True
    acados_ocp_solver = AcadosOcpSolver(ocp)

    optimal_value_grad = np.zeros(np_test)
    optimal_value = np.zeros(np_test)

    pi = np.zeros(np_test)
    for i, p in enumerate(p_test):

        for n in range(N_horizon+1):
            acados_ocp_solver.set(n, 'p', p)
        pi[i] = acados_ocp_solver.solve_for_x0(x0)[0]
        optimal_value[i] = acados_ocp_solver.get_cost()
        optimal_value_grad[i] = acados_ocp_solver.eval_and_get_optimal_value_gradient("params_global").item()

    # evaluate cost gradient
    optimal_value_grad_via_fd = np.gradient(optimal_value, delta_p)
    cost_reconstructed_np_grad = np.cumsum(optimal_value_grad_via_fd) * delta_p + optimal_value[0]

    plot_cost_gradient_results(p_test, optimal_value, optimal_value_grad, optimal_value_grad_via_fd, cost_reconstructed_np_grad)
    assert np.median(np.abs(optimal_value_grad - optimal_value_grad_via_fd)) <= 1e-1



if __name__ == "__main__":
    main()
