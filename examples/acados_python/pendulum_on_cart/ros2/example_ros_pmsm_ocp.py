import argparse
import os
import sys

import numpy as np
import scipy.linalg

from acados_template import AcadosOcp, AcadosOcpSolver
from acados_template.ros2 import AcadosOcpRosOptions

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PMSM_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "..", "pmsm_example"))
if PMSM_DIR not in sys.path:
    sys.path.insert(0, PMSM_DIR)

import main as pmsm_main


def create_pmsm_ocp(export_dir: str) -> AcadosOcp:
    """Create a ROS2 OCP for the PMSM benchmark with phi, g, g_e, slack, and p constraints."""
    os.makedirs(export_dir, exist_ok=True)
    pmsm_main.FORMULATION = 1

    ocp = AcadosOcp()
    model = pmsm_main.export_pmsm_model()
    ocp.model = model

    nx = model.x.rows()
    nu = model.u.rows()
    ny = nx + nu
    ny_e = nx
    N = pmsm_main.N
    Ts = pmsm_main.Ts
    Tf = N * Ts

    ocp.solver_options.N_horizon = N
    ocp.solver_options.tf = Tf

    # LINEAR_LS cost matching pmsm_example (FORMULATION = 1)
    Q = np.diag([pmsm_main.wd * pmsm_main.Weight_TUNING, pmsm_main.wq * pmsm_main.Weight_TUNING])
    Q_e = np.diag([
        pmsm_main.wd * pmsm_main.Weight_E_TUNING * Tf / N,
        pmsm_main.wq * pmsm_main.Weight_E_TUNING * Tf / N,
    ])
    R = pmsm_main.INPUT_REG * np.eye(nu)

    ocp.cost.W = scipy.linalg.block_diag(Q, R)
    ocp.cost.W_e = Q_e

    Vu = np.zeros((ny, nu))
    Vu[2, 0] = 1.0
    Vu[3, 1] = 1.0
    ocp.cost.Vu = Vu

    Vx = np.zeros((ny, nx))
    Vx[:nx, :nx] = np.eye(nx)
    ocp.cost.Vx = Vx
    ocp.cost.Vx_e = np.eye(nx)

    ocp.cost.yref = np.zeros((ny,))
    ocp.cost.yref_e = np.zeros((ny_e,))

    # Stage & terminal polytopic constraints: lg <= D u + C x <= ug, lg_e <= C_e x <= ug_e
    con_dc = pmsm_main.get_general_constraints_DC()
    con_dc_e = pmsm_main.get_general_terminal_constraints_DC()
    ocp.constraints.D = con_dc["D"]
    ocp.constraints.C = con_dc["C"]
    ocp.constraints.lg = con_dc["lg"]
    ocp.constraints.ug = con_dc["ug"]
    ocp.constraints.C_e = con_dc_e["Ce"]
    ocp.constraints.lg_e = con_dc_e["lge"]
    ocp.constraints.ug_e = con_dc_e["uge"]

    # Convex-over-nonlinear (phi) torque & voltage constraints with slacks across 0, stage, e
    lphi = np.array([0.0, -1e9])
    uphi = np.array([0.0, pmsm_main.u_max**2 / 3.0])
    ocp.constraints.lphi_0 = lphi
    ocp.constraints.uphi_0 = uphi
    ocp.constraints.lphi = lphi
    ocp.constraints.uphi = uphi
    ocp.constraints.lphi_e = lphi
    ocp.constraints.uphi_e = uphi

    ocp.constraints.idxsphi_0 = np.array([0])
    ocp.constraints.idxsphi = np.array([0])
    ocp.constraints.idxsphi_e = np.array([0])

    # Slack penalties (ns_0 = ns = ns_e = 1)
    ocp.cost.zl_0 = pmsm_main.SLACK_E_TUNING * np.ones((1,))
    ocp.cost.zu_0 = pmsm_main.SLACK_E_TUNING * np.ones((1,))
    ocp.cost.Zl_0 = pmsm_main.SLACK_TUNINGHessian * np.ones((1,))
    ocp.cost.Zu_0 = pmsm_main.SLACK_TUNINGHessian * np.ones((1,))

    ocp.cost.zl = pmsm_main.SLACK_TUNING * np.ones((1,))
    ocp.cost.zu = pmsm_main.SLACK_TUNING * np.ones((1,))
    ocp.cost.Zl = pmsm_main.SLACK_TUNINGHessian * np.ones((1,))
    ocp.cost.Zu = pmsm_main.SLACK_TUNINGHessian * np.ones((1,))

    ocp.cost.zl_e = pmsm_main.SLACK_E_TUNING * np.ones((1,))
    ocp.cost.zu_e = pmsm_main.SLACK_E_TUNING * np.ones((1,))
    ocp.cost.Zl_e = pmsm_main.SLACK_TUNINGHessian * np.ones((1,))
    ocp.cost.Zu_e = pmsm_main.SLACK_TUNINGHessian * np.ones((1,))

    ocp.constraints.x0 = pmsm_main.x0Start.astype(float)
    ocp.parameter_values = np.array([pmsm_main.w_val, 0.0, 0.0, pmsm_main.tau_wal])

    # Solver options (SQP converges in a single call without offline RTI warmstart loop)
    ocp.solver_options.qp_solver = "FULL_CONDENSING_HPIPM"
    ocp.solver_options.hessian_approx = "GAUSS_NEWTON"
    ocp.solver_options.integrator_type = "IRK"
    ocp.solver_options.sim_method_num_stages = 1
    ocp.solver_options.qp_solver_tol_eq = 1e-4
    ocp.solver_options.qp_solver_tol_ineq = 1e-4
    ocp.solver_options.qp_solver_tol_comp = 1e-4
    ocp.solver_options.nlp_solver_type = "SQP"
    ocp.solver_options.nlp_solver_max_iter = 200

    # ROS2 options
    ocp.ros_opts = AcadosOcpRosOptions()
    ocp.ros_opts.package_name = "pmsm_ocp"
    ocp.ros_opts.generated_code_dir = export_dir
    ocp.ros_opts.publish_control_sequence = True
    ocp.ros_opts.threads = 2

    ocp.code_gen_options.code_export_directory = os.path.join(export_dir, "c_generated_code")
    ocp.code_gen_options.json_file = "acados_ocp.json"
    return ocp


def main():
    parser = argparse.ArgumentParser(
        description="Generate a ROS2 OCP package for the PMSM problem with phi, g, g_e, and slack constraints."
    )
    parser.add_argument(
        "--export-dir",
        default=os.path.join(SCRIPT_DIR, "ros_generated_pmsm_ocp"),
        help="Target directory for generated ROS2 packages and expected test data.",
    )
    args = parser.parse_args()

    export_dir = os.path.abspath(args.export_dir)
    ocp = create_pmsm_ocp(export_dir=export_dir)
    solver = AcadosOcpSolver(ocp)

    status = solver.solve()
    max_sl = max(float(np.max(solver.get(i, "sl"))) for i in range(ocp.solver_options.N_horizon + 1))
    max_su = max(float(np.max(solver.get(i, "su"))) for i in range(ocp.solver_options.N_horizon + 1))
    print(f"Solver status: {status}, SQP iters: {solver.get_stats('sqp_iter')}, residuals: {solver.get_residuals()}")
    print(f"Max slacks across horizon: sl={max_sl:.4e}, su={max_su:.4e}")

    iterate = solver.get_iterate()
    expected_u_file = os.path.join(export_dir, "expected_control_sequence.npy")
    np.save(expected_u_file, np.array(iterate.u))
    print(f"Saved expected control sequence to {expected_u_file}")


if __name__ == "__main__":
    main()
