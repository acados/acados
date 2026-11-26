#
# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# Licensed under the 2-Clause BSD License.

import argparse
import os
import sys

import numpy as np

from acados_template import AcadosOcp, AcadosOcpSolver
from acados_template.ros2 import AcadosOcpRosOptions

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
if SCRIPT_DIR not in sys.path:
    sys.path.insert(0, SCRIPT_DIR)

from setup_ocp_solver import setup_ocp
from utils import load_trajectory


def create_fixed_wing_ocp(
    export_dir: str,
    N_horizon: int = 20,
    dt: float = 0.1,
) -> AcadosOcp:
    """Create a ROS2 OCP package for the fixed-wing aircraft MPC problem."""
    os.makedirs(export_dir, exist_ok=True)
    tf = dt * N_horizon

    traj_file = os.path.join(SCRIPT_DIR, "data", "trajectory_tilt20deg_r50m.csv")
    aero_file = os.path.join(SCRIPT_DIR, "data", "aerodynamic_coefficients.json")
    x_ref, u_ref, _ = load_trajectory(traj_file, dt=dt, number_of_orbits=1)

    ocp = setup_ocp(
        x0=x_ref[0],
        N_horizon=N_horizon,
        Tf=tf,
        RTI=False,
        soft_constraints=True,
        aero_coeff_file=aero_file,
    )
    ocp.cost.yref = np.concatenate([x_ref[5], u_ref[5]])
    ocp.cost.yref_e = x_ref[10]
    ocp.solver_options.nlp_solver_max_iter = 50

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
        description="Generate a ROS2 OCP package for the fixed-wing aircraft MPC problem."
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
    max_sl = max(float(np.max(solver.get(i, "sl"), initial=0.0)) for i in range(ocp.solver_options.N_horizon + 1))
    max_su = max(float(np.max(solver.get(i, "su"), initial=0.0)) for i in range(ocp.solver_options.N_horizon + 1))
    print(f"Solver status: {status}, SQP iters: {solver.get_stats('sqp_iter')}, residuals: {solver.get_residuals()}")
    print(f"Max active slacks across horizon: sl={max_sl:.4e}, su={max_su:.4e}")

    iterate = solver.get_iterate()
    expected_u_file = os.path.join(export_dir, "expected_control_sequence.npy")
    np.save(expected_u_file, np.array(iterate.u))
    print(f"Saved expected control sequence to {expected_u_file}")


if __name__ == "__main__":
    main()
