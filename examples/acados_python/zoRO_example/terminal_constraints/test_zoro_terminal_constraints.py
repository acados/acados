# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# Licensed under the 2-Clause BSD License.

# Regression test for the zoRO custom update with terminal general constraints
# (linear `C_e x <= ug_e` and nonlinear `h_e(x) <= uh_e`).
# The terminal back-off has to be computed with the terminal constraint Jacobians,
# see https://github.com/acados/acados/issues/1968

import numpy as np
import casadi as ca
from acados_template import AcadosOcp, AcadosOcpSolver, AcadosModel, ZoroDescription

# scalar double integrator-like system: x_{k+1} = x_k + u_k, nx = nu = 1
NX = 1
NU = 1
N = 8
GAMMA = 2.0
W = 0.01
K_FB = 0.2
UG_E = 3.0
UH_E = 3.0


def compute_terminal_P():
    # zoRO uncertainty propagation with linear dynamics and fixed feedback gain u = -K x
    A_cl = 1.0 - K_FB  # A - B K
    P = 0.0
    for _ in range(N):
        P = A_cl * P * A_cl + W
    return P


def create_solver(terminal_type: str, build: bool = True) -> AcadosOcpSolver:
    ocp = AcadosOcp()

    model = AcadosModel()
    model.name = f"zoro_terminal_{terminal_type}"
    x = ca.SX.sym("x", NX)
    u = ca.SX.sym("u", NU)
    model.x = x
    model.u = u
    model.disc_dyn_expr = x + u
    if terminal_type == "h":
        # nonlinear terminal constraint which is linear in the region of interest
        model.con_h_expr_e = ca.vertcat(x, x**3)
    ocp.model = model

    ocp.solver_options.N_horizon = N
    ocp.solver_options.tf = float(N)
    ocp.solver_options.integrator_type = "DISCRETE"
    ocp.solver_options.nlp_solver_type = "SQP"
    ocp.solver_options.qp_solver = "PARTIAL_CONDENSING_HPIPM"
    ocp.solver_options.hessian_approx = "GAUSS_NEWTON"

    # cost pushes the terminal state upwards, against the terminal constraint
    ocp.cost.cost_type = "LINEAR_LS"
    ocp.cost.cost_type_e = "LINEAR_LS"
    ocp.cost.Vx = np.array([[1.0], [0.0]])
    ocp.cost.Vu = np.array([[0.0], [1.0]])
    ocp.cost.W = np.diag([1.0, 1e-2])
    ocp.cost.yref = np.array([10.0, 0.0])
    ocp.cost.Vx_e = np.array([[1.0]])
    ocp.cost.W_e = np.array([[1.0]])
    ocp.cost.yref_e = np.array([10.0])

    ocp.constraints.x0 = np.zeros(NX)

    # stage general constraint with a Jacobian that differs from the terminal one,
    # ng = 2 != ng_e = 1
    ocp.constraints.C = np.array([[0.0], [0.0]])
    ocp.constraints.D = np.array([[1.0], [-1.0]])
    ocp.constraints.lg = np.array([-100.0, -100.0])
    ocp.constraints.ug = np.array([100.0, 100.0])

    zoro_description = ZoroDescription()
    zoro_description.backoff_scaling_gamma = GAMMA
    zoro_description.P0_mat = np.zeros((NX, NX))
    zoro_description.fdbk_K_mat = np.array([[K_FB]])
    zoro_description.W_mat = np.array([[W]])

    if terminal_type == "g":
        ocp.constraints.C_e = np.array([[1.0]])
        ocp.constraints.lg_e = np.array([-100.0])
        ocp.constraints.ug_e = np.array([UG_E])
        zoro_description.idx_ug_e_t = [0]
    elif terminal_type == "h":
        # nh_e = 2 > nh = 0
        ocp.constraints.lh_e = np.array([-100.0, -1e6])
        ocp.constraints.uh_e = np.array([UH_E, 1e6])
        zoro_description.idx_uh_e_t = [0]
    else:
        raise ValueError(terminal_type)

    ocp.zoro_description = zoro_description

    ocp.solver_options.custom_update_filename = "custom_update_function.c"
    ocp.solver_options.custom_update_header_filename = "custom_update_function.h"
    ocp.solver_options.custom_update_copy = False
    ocp.solver_options.custom_templates = [
        ("custom_update_function_zoro_template.in.c", "custom_update_function.c"),
        ("custom_update_function_zoro_template.in.h", "custom_update_function.h"),
    ]
    ocp.code_gen_options.code_export_directory = f"c_generated_code_{model.name}"
    ocp.code_gen_options.json_file = f"{model.name}.json"

    return AcadosOcpSolver(ocp, build=build, generate=build)


def run_test(terminal_type: str):
    solver = create_solver(terminal_type)

    for _ in range(20):
        status = solver.solve()
        assert status == 0, f"solver returned status {status}"
        solver.custom_update(np.array([]))

    status = solver.solve()
    assert status == 0, f"solver returned status {status}"

    x_N = solver.get(N, "x")[0]
    P_N = compute_terminal_P()
    upper = UG_E if terminal_type == "g" else UH_E
    x_N_expected = upper - GAMMA * np.sqrt(P_N)

    print(f"terminal constraint {terminal_type}_e: x_N = {x_N:.6f}, expected {x_N_expected:.6f}")
    assert abs(x_N - x_N_expected) < 1e-6, \
        f"terminal constraint {terminal_type}_e not tightened correctly: got x_N = {x_N}, expected {x_N_expected}"


if __name__ == "__main__":
    run_test("g")
    run_test("h")
