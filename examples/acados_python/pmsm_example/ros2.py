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

from main import setup_pmsm_ocp


def create_pmsm_ocp(export_dir: str) -> AcadosOcp:
    """Create a ROS2 OCP for the PMSM benchmark with phi, g, g_e, slack, and p constraints."""
    os.makedirs(export_dir, exist_ok=True)

    # SQP converges in a single call without the offline RTI warmstart loop
    ocp = setup_pmsm_ocp(
        qp_solver="FULL_CONDENSING_HPIPM",
        nlp_solver_type="SQP",
    )
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
