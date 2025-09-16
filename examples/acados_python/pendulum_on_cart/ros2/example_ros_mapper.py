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
import os
script_dir = os.path.dirname(os.path.realpath(__file__))
common_path = os.path.join(script_dir, '..', 'common')
sys.path.insert(0, os.path.abspath(common_path))

from acados_template import AcadosOcpSolver
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
    
    AcadosOcpSolver(ocp, json_file = str(os.path.join(export_dir, 'acados_ocp.json')))
    ros_mapper.generate()
    
    
if __name__ == "__main__":
    main()