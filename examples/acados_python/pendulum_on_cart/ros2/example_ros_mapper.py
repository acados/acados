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

from acados_template import AcadosOcpSolver, AcadosOcp
from acados_template.ros2 import RosTopicMapper, RosTopicMsgOutput, build_default_control, build_default_state
from acados_template.ros2.default_msgs import GEOMETRY_MSGS_TWIST, GEOMETRY_MSGS_POSE

from example_ros_minimal_ocp import create_minimal_ocp


def main():
    Fmax = 80
    Tf_ocp = 1.0
    Tf_sim = 0.05
    N = 20
    export_dir = os.path.join(script_dir, 'generated')
    c_generated_code_base = os.path.join(export_dir, "c_generated_code")

    ocp = create_minimal_ocp(export_dir, N, Tf_ocp, Fmax)
    ocp.code_export_directory = c_generated_code_base + "_ocp"


    # --- SETUP MESSAGES ---
    # We create all messages here that the mapper should subscribe and publish.
    # e.g. /pose.position.x -> /ocp_state.x[0]
    #      /ocp_control.u[0] -> /cmd_vel.linear.x

    # setup control input messages
    in_ocp_ctrl  = build_default_control(ocp, direction_out=False)

    # setup state output messages
    out_ocp_state = build_default_state(ocp, direction_out=True)
    out_ocp_state.exec_topic = "/pose"
    out_ocp_state.mapping = [
        (f"pose.position.x", "x[0]"),
        (f"pose.position.y", "x[1]"),
        (f"pose.orientation.z", "x[2]")
    ]

    # setup twist msg
    twist = RosTopicMsgOutput.from_msg(GEOMETRY_MSGS_TWIST)
    twist.topic_name = "/cmd_vel"
    twist.exec_topic = ocp.ros_opts.control_topic
    twist.mapping = [
        (f"{in_ocp_ctrl.topic_name}.u[0]", "linear.x")
    ]

    # setup pose msg
    pose = GEOMETRY_MSGS_POSE
    pose.topic_name = "/pose"

    # --- GENERATE MAPPER ---
    ros_mapper = RosTopicMapper()
    ros_mapper.package_name = "ocp_mapper"
    ros_mapper.generated_code_dir = export_dir
    ros_mapper.in_msgs = [
        in_ocp_ctrl,
        pose
    ]
    ros_mapper.out_msgs = [
        out_ocp_state,
        twist
    ]

    json_file = str(os.path.join(export_dir, 'acados_ocp.json'))
    ocp.code_gen_options.json_file = json_file
    # NOTE: this is just to test: dump and load.
    ocp.dump_to_json()
    ocp = AcadosOcp.from_json(json_file)

    AcadosOcpSolver(ocp)
    ros_mapper.generate()


if __name__ == "__main__":
    main()