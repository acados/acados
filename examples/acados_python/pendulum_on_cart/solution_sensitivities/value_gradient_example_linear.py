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

import numpy as np
from sensitivity_utils import plot_cost_gradient_results
from acados_template import AcadosModel, AcadosOcp, AcadosOcpSolver
from casadi.tools import entry, struct_symSX
import casadi as ca


PARAM_VALUE_DICT = {
    "A": np.array([[1.0, 0.25], [0.0, 1.0]]),
    "B": np.array([[0.03125], [0.25]]),
    "Q": np.identity(2),
    "R": np.identity(1),
    "b": np.array([[0.2], [0.2]]),
    "f": np.array([[0.0], [0.0], [0.0]]),
    "V_0": np.array([1e-3]),
}


def find_param_in_p_or_p_global(param_name: list[str], model: AcadosModel) -> list:
    if model.p == []:
        return {key: model.p_global[key] for key in param_name}
    elif model.p_global is None:
        return {key: model.p[key] for key in param_name}
    else:
        return {
            key: (model.p[key] if key in model.p.keys() else model.p_global[key])
            for key in param_name
        }


def disc_dyn_expr(model: AcadosModel):
    x = model.x
    u = model.u

    param = find_param_in_p_or_p_global(["A", "B", "b"], model)

    return param["A"] @ x + param["B"] @ u + param["b"]


def cost_expr_ext_cost(model: AcadosModel):
    x = model.x
    u = model.u
    param = find_param_in_p_or_p_global(["Q", "R", "f"], model)

    return 0.5 * (
        ca.transpose(x) @ param["Q"] @ x
        + ca.transpose(u) @ param["R"] @ u
        + ca.transpose(param["f"]) @ ca.vertcat(x, u)
    )


def cost_expr_ext_cost_0(model: AcadosModel):
    param = find_param_in_p_or_p_global(["V_0"], model)

    return param["V_0"] + cost_expr_ext_cost(model)


def cost_expr_ext_cost_e(model: AcadosModel, param: dict[str, np.ndarray]):

    x = model.x

    return 0.5 * ca.mtimes(
        [
            ca.transpose(x),
            param["Q"],
            x,
        ]
    )


def export_parametric_ocp(
    param: dict[str, np.ndarray],
    name: str = "lti",
    learnable_params: list[str] = [],
) -> AcadosOcp:
    ocp = AcadosOcp()

    ocp.model.name = name

    ocp.model.x = ca.SX.sym("x", 2)
    ocp.model.u = ca.SX.sym("u", 1)

    ocp.solver_options.N_horizon = 4
    ocp.solver_options.tf = 8
    ocp.solver_options.integrator_type = 'DISCRETE'

    # Add learnable parameters to p_global
    if len(learnable_params) != 0:
        ocp.model.p_global = struct_symSX(
            [entry(key, shape=param[key].shape) for key in learnable_params]
        )
        ocp.p_global_values = np.concatenate(
            [param[key].T.reshape(-1, 1) for key in learnable_params]
        ).flatten()

    # Add non_learnable parameters to p (stage-wise parameters)
    non_learnable_params = [key for key in param.keys() if key not in learnable_params]
    if len(non_learnable_params) != 0:
        ocp.model.p = struct_symSX(
            [entry(key, shape=param[key].shape) for key in non_learnable_params]
        )
        ocp.parameter_values = np.concatenate(
            [param[key].T.reshape(-1, 1) for key in non_learnable_params]
        ).flatten()

    print("learnable_params", learnable_params)
    print("non_learnable_params", non_learnable_params)

    ocp.model.disc_dyn_expr = disc_dyn_expr(ocp.model)

    ocp.cost.cost_type_0 = "EXTERNAL"
    ocp.model.cost_expr_ext_cost_0 = cost_expr_ext_cost_0(ocp.model)

    ocp.cost.cost_type = "EXTERNAL"
    ocp.model.cost_expr_ext_cost = cost_expr_ext_cost(ocp.model)

    ocp.cost.cost_type_e = "EXTERNAL"
    ocp.model.cost_expr_ext_cost_e = cost_expr_ext_cost_e(ocp.model, param)

    ocp.constraints.idxbx_0 = np.array([0, 1])
    ocp.constraints.lbx_0 = np.array([-1.0, -1.0])
    ocp.constraints.ubx_0 = np.array([1.0, 1.0])

    ocp.constraints.idxbx = np.array([0, 1])
    ocp.constraints.lbx = np.array([-0.0, -1.0])
    ocp.constraints.ubx = np.array([+1.0, +1.0])

    ocp.constraints.idxsbx = np.array([0])
    ocp.cost.zl = np.array([1e2])
    ocp.cost.zu = np.array([1e2])
    ocp.cost.Zl = np.diag([0])
    ocp.cost.Zu = np.diag([0])

    ocp.constraints.idxbu = np.array([0])
    ocp.constraints.lbu = np.array([-1.0])
    ocp.constraints.ubu = np.array([+1.0])

    if isinstance(ocp.model.p, struct_symSX):
        ocp.model.p = ocp.model.p.cat if ocp.model.p is not None else []
    if isinstance(ocp.model.p_global, struct_symSX):
        ocp.model.p_global = (
            ocp.model.p_global.cat if ocp.model.p_global is not None else None
        )

    return ocp


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
    ocp.solver_options.nlp_solver_type = "SQP"
    acados_ocp_solver = AcadosOcpSolver(ocp)

    optimal_value_grad = np.zeros((np_test,))
    optimal_value = np.zeros((np_test,))

    dt = ocp.solver_options.tf/ocp.solver_options.N_horizon
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
    print(f"Median difference between value function gradient obtained by acados and via FD is {median_diff:e} should be < {test_tol:e}.")
    assert median_diff <= test_tol

if __name__ == "__main__":
    main()
