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
sys.path.insert(0, '../common')

from acados_template import AcadosOcp, AcadosOcpSolver, AcadosModel, latexify_plot
import numpy as np
import matplotlib.pyplot as plt
from utils import plot_pendulum

from policy_gradient_example import export_parametric_ocp, plot_cost_gradient_results



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

    sens_cost = np.zeros(np_test)
    cost_values = np.zeros(np_test)

    pi = np.zeros(np_test)
    for i, p in enumerate(p_test):

        for n in range(N_horizon+1):
            acados_ocp_solver.set(n, 'p', p)
        pi[i] = acados_ocp_solver.solve_for_x0(x0)[0]
        cost_values[i] = acados_ocp_solver.get_cost()
        sens_cost[i] = acados_ocp_solver.get_optimal_value_gradient("params_global").item()

    # evaluate cost gradient
    np_cost_grad = np.gradient(cost_values, delta_p)
    cost_reconstructed_np_grad = np.cumsum(np_cost_grad) * delta_p + cost_values[0]

    plot_cost_gradient_results(p_test, cost_values, sens_cost, np_cost_grad, cost_reconstructed_np_grad)
    assert np.allclose(sens_cost, np_cost_grad, atol=5e2, rtol=1e-2)


if __name__ == "__main__":
    main()
