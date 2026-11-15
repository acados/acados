import argparse
import os
import sys

import casadi as ca
import numpy as np
import scipy.linalg as la

from acados_template import AcadosOcp, AcadosOcpSolver
from acados_template.ros2 import AcadosOcpRosOptions

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
AIRCRAFT_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "..", "fixed_wing_aircraft"))
if AIRCRAFT_DIR not in sys.path:
    sys.path.insert(0, AIRCRAFT_DIR)

from model import export_aircraft_model
from utils import load_trajectory


def create_fixed_wing_ocp(
    export_dir: str,
    N_horizon: int = 20,
    dt: float = 0.1,
) -> AcadosOcp:
    """Create a complex fixed-wing aircraft OCP with slacks and ROS2 options."""
    os.makedirs(export_dir, exist_ok=True)
    tf = dt * N_horizon

    traj_file = os.path.join(AIRCRAFT_DIR, "data", "trajectory_tilt20deg_r50m.csv")
    aero_file = os.path.join(AIRCRAFT_DIR, "data", "aerodynamic_coefficients.json")
    x_ref, u_ref, _ = load_trajectory(traj_file, dt=dt, number_of_orbits=1)

    ocp = AcadosOcp()
    model = export_aircraft_model(aero_coeff_file=aero_file)
    ocp.model = model

    nx = model.x.rows()
    nu = model.u.rows()

    # Non-zero runtime wind parameters -> generates Parameters.msg & /ocp_params
    ocp.parameter_values = np.array([2.5, -1.5, 0.3])

    # NONLINEAR_LS costs across initial (_0), stage, and terminal (_e) nodes
    ocp.cost.cost_type_0 = "NONLINEAR_LS"
    ocp.cost.cost_type = "NONLINEAR_LS"
    ocp.cost.cost_type_e = "NONLINEAR_LS"

    ocp.model.cost_y_expr_0 = ca.vertcat(model.x, model.u)
    ocp.model.cost_y_expr = ca.vertcat(model.x, model.u)
    ocp.model.cost_y_expr_e = model.x

    q_diag = np.diag([0.1] * 3 + [0.1] * 3 + [1.0] * 9 + [1.0] * 3 + [0.1] * 2 + [1.0] * 2)
    r_diag = np.diag([0.1, 0.1, 50.0, 50.0])

    ocp.cost.W_0 = la.block_diag(q_diag, r_diag)
    ocp.cost.W = la.block_diag(q_diag, r_diag)
    ocp.cost.W_e = (tf / N_horizon) * q_diag

    ocp.cost.yref_0 = np.concatenate([x_ref[1], u_ref[0]])
    ocp.cost.yref = np.concatenate([x_ref[5], u_ref[5]])
    ocp.cost.yref_e = x_ref[10]

    # Control bounds with soft slack on dthrust[0]
    ocp.constraints.lbu = np.array([-3.0, -5.0, np.deg2rad(-19), np.deg2rad(-19)])
    ocp.constraints.ubu = np.array([3.0, 5.0, np.deg2rad(19), np.deg2rad(19)])
    ocp.constraints.idxbu = np.array([0, 1, 2, 3])
    ocp.constraints.idxsbu = np.array([0])
    ocp.constraints.lsbu = np.array([0.01])
    ocp.constraints.usbu = np.array([0.02])

    # Stage state bounds with tightened upper thrust[0] limit so stage slack activates
    ocp.constraints.lbx = np.array([-1.0, -1.0, -1.0, 0.0, 0.0, np.deg2rad(-32), np.deg2rad(-25)])
    ocp.constraints.ubx = np.array([1.0, 1.0, 1.0, 4.0, 20.0, np.deg2rad(32), np.deg2rad(25)])
    ocp.constraints.idxbx = np.array([15, 16, 17, 18, 19, 20, 21])
    ocp.constraints.idxsbx = np.array([0, 3])
    ocp.constraints.lsbx = np.array([0.01, 0.02])
    ocp.constraints.usbx = np.array([0.02, 0.01])

    # Terminal state bounds with soft slack on flap[0]
    ocp.constraints.lbx_e = np.array([np.deg2rad(-25), np.deg2rad(-20)])
    ocp.constraints.ubx_e = np.array([np.deg2rad(25), np.deg2rad(20)])
    ocp.constraints.idxbx_e = np.array([20, 21])
    ocp.constraints.idxsbx_e = np.array([0])
    ocp.constraints.lsbx_e = np.array([0.01])
    ocp.constraints.usbx_e = np.array([0.02])

    # Soft polytopic constraint: lg <= D u + C x <= ug
    d_mat = np.zeros((1, nu))
    d_mat[0, 0] = 1.0
    d_mat[0, 1] = -1.0
    c_mat = np.zeros((1, nx))
    c_mat[0, 18] = 0.05
    ocp.constraints.D = d_mat
    ocp.constraints.C = c_mat
    ocp.constraints.lg = np.array([-0.15])
    ocp.constraints.ug = np.array([0.15])
    ocp.constraints.idxsg = np.array([0])
    ocp.constraints.lsg = np.array([0.01])
    ocp.constraints.usg = np.array([0.02])

    # Initial nonlinear constraints (all 5 components softened)
    ocp.model.con_h_expr_0 = model.con_h_expr
    ocp.constraints.lh_0 = np.array([12.0, -5.0, -5.0, np.deg2rad(-5), np.deg2rad(-8)])
    ocp.constraints.uh_0 = np.array([18.0, 5.0, 5.0, np.deg2rad(10), np.deg2rad(8)])
    ocp.constraints.idxsh_0 = np.arange(5)
    ocp.constraints.lsh_0 = np.array([0.01, 0.0, 0.0, 0.0, 0.0])
    ocp.constraints.ush_0 = np.array([0.0, 0.02, 0.0, 0.0, 0.0])

    # Stage nonlinear constraints (tightened airspeed bound activates slack)
    ocp.constraints.lh = np.array([13.0, -5.0, -5.0, np.deg2rad(-5), np.deg2rad(-8)])
    ocp.constraints.uh = np.array([13.5, 5.0, 5.0, np.deg2rad(10), np.deg2rad(8)])
    ocp.constraints.idxsh = np.arange(5)
    ocp.constraints.lsh = np.array([0.02, 0.01, 0.0, 0.0, 0.0])
    ocp.constraints.ush = np.array([0.0, 0.0, 0.01, 0.0, 0.0])

    # Terminal nonlinear constraints (all 5 components softened)
    ocp.model.con_h_expr_e = model.con_h_expr
    ocp.constraints.lh_e = ocp.constraints.lh
    ocp.constraints.uh_e = ocp.constraints.uh
    ocp.constraints.idxsh_e = np.arange(5)
    ocp.constraints.lsh_e = np.array([0.01, 0.0, 0.0, 0.0, 0.0])
    ocp.constraints.ush_e = np.array([0.0, 0.01, 0.0, 0.0, 0.0])

    # Non-uniform L1 (z*) and L2 (Z*) slack penalties across all horizon regions:
    # ns_0 = 1 (sbu) + 1 (sg) + 5 (sh_0) = 7
    # ns   = 1 (sbu) + 2 (sbx) + 1 (sg) + 5 (sh) = 9
    # ns_e = 1 (sbx_e) + 5 (sh_e) = 6
    ocp.cost.zl_0 = np.linspace(30.0, 60.0, 7)
    ocp.cost.zu_0 = np.linspace(35.0, 65.0, 7)
    ocp.cost.Zl_0 = np.linspace(2.0, 8.0, 7)
    ocp.cost.Zu_0 = np.linspace(3.0, 9.0, 7)

    ocp.cost.zl = np.linspace(25.0, 65.0, 9)
    ocp.cost.zu = np.linspace(20.0, 60.0, 9)
    ocp.cost.Zl = np.linspace(1.5, 9.5, 9)
    ocp.cost.Zu = np.linspace(2.5, 10.5, 9)

    ocp.cost.zl_e = np.linspace(40.0, 70.0, 6)
    ocp.cost.zu_e = np.linspace(45.0, 75.0, 6)
    ocp.cost.Zl_e = np.linspace(4.0, 10.0, 6)
    ocp.cost.Zu_e = np.linspace(5.0, 11.0, 6)

    # Initial state and solver options
    ocp.constraints.x0 = x_ref[0]
    ocp.solver_options.N_horizon = N_horizon
    ocp.solver_options.tf = tf
    ocp.solver_options.nlp_solver_type = "SQP"
    ocp.solver_options.nlp_solver_max_iter = 50
    ocp.solver_options.integrator_type = "IRK"
    ocp.solver_options.collocation_type = "GAUSS_RADAU_IIA"
    ocp.solver_options.qp_solver = "PARTIAL_CONDENSING_HPIPM"
    ocp.solver_options.hessian_approx = "GAUSS_NEWTON"
    ocp.solver_options.levenberg_marquardt = 1e-4
    ocp.solver_options.sim_method_newton_iter = 10

    # ROS2 code generation options (multithreaded executor + control sequence topic)
    ocp.ros_opts = AcadosOcpRosOptions()
    ocp.ros_opts.package_name = "fixed_wing_aircraft_ocp"
    ocp.ros_opts.generated_code_dir = export_dir
    ocp.ros_opts.publish_control_sequence = True
    ocp.ros_opts.threads = 2

    ocp.code_gen_options.code_export_directory = os.path.join(export_dir, "c_generated_code")
    ocp.code_gen_options.json_file = "acados_ocp.json"
    return ocp


def main():
    parser = argparse.ArgumentParser(
        description="Generate a ROS2 OCP package for a slacked fixed-wing aircraft MPC problem."
    )
    parser.add_argument(
        "--export-dir",
        default=os.path.join(SCRIPT_DIR, "ros_generated_fixed_wing_ocp"),
        help="Target directory for generated ROS2 packages and expected test data.",
    )
    args = parser.parse_args()

    export_dir = os.path.abspath(args.export_dir)
    ocp = create_fixed_wing_ocp(export_dir=export_dir)
    solver = AcadosOcpSolver(ocp)

    status = solver.solve()
    max_sl = max(float(np.max(solver.get(i, "sl"))) for i in range(ocp.solver_options.N_horizon + 1))
    max_su = max(float(np.max(solver.get(i, "su"))) for i in range(ocp.solver_options.N_horizon + 1))
    print(f"Solver status: {status}, SQP iters: {solver.get_stats('sqp_iter')}, residuals: {solver.get_residuals()}")
    print(f"Max active slacks across horizon: sl={max_sl:.4e}, su={max_su:.4e}")

    iterate = solver.get_iterate()
    expected_u_file = os.path.join(export_dir, "expected_control_sequence.npy")
    np.save(expected_u_file, np.array(iterate.u))
    print(f"Saved expected control sequence to {expected_u_file}")


if __name__ == "__main__":
    main()
