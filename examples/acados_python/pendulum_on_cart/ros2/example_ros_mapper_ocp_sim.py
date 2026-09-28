#
# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# Licensed under the 2-Clause BSD License.

import sys
import os
script_dir = os.path.dirname(os.path.realpath(__file__))
common_path = os.path.join(script_dir, '..', 'common')
sys.path.insert(0, os.path.abspath(common_path))

from pprint import pprint

from acados_template import AcadosOcpSolver, AcadosSimSolver
from acados_template.ros2 import RosTopicMapper

from example_ros_minimal_ocp import create_minimal_ocp
from example_ros_minimal_sim import create_minimal_sim



def main():
    Fmax = 80
    Tf_ocp = 1.0
    Tf_sim = 0.05
    N = 20
    export_dir = os.path.join(script_dir, 'generated')
    c_generated_code_base = os.path.join(export_dir, "c_generated_code")

    ocp = create_minimal_ocp(export_dir, N, Tf_ocp, Fmax)
    ocp.code_export_directory = c_generated_code_base + "_ocp"
    ocp_solver = AcadosOcpSolver(ocp, json_file = str(os.path.join(export_dir, 'acados_ocp.json')))

    sim = create_minimal_sim(export_dir, Tf_sim)
    sim.code_export_directory = c_generated_code_base + "_sim"
    sim_solver = AcadosSimSolver(sim, json_file = str(os.path.join(export_dir, 'acados_sim.json')))

    ros_mapper = RosTopicMapper.from_instances(ocp_solver, sim_solver)
    ros_mapper.package_name = "sim_ocp_mapper"
    ros_mapper.generated_code_dir = export_dir

    ros_mapper.generate()

    
if __name__ == "__main__":
    main()