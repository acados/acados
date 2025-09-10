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

from typing import Optional, Union
from .utils import ArchType, ControlLoopExec, AcadosRosBaseOptions


class RosTopicMsg:
    def __init__(self):
        self.__topic_name: str = ""
        self.__msg_type: str = ""
        self.__field_names: list[str] = list()
        self.__field_types: list[str] = list()

    @property
    def topic_name(self) -> str:
        return self.__topic_name

    @property
    def msg_type(self) -> str:
        return self.__msg_type

    @property
    def field_names(self) -> list[str]:
        return self.__field_names

    @property
    def field_types(self) -> list[str]:
        return self.__field_types
    
    @topic_name.setter
    def topic_name(self, topic_name: str):
        if not isinstance(topic_name, str):
            raise TypeError('Invalid topic_name value, expected str.\n')
        self.__topic_name = topic_name
        
    @msg_type.setter
    def msg_type(self, msg_type: str):
        if not isinstance(msg_type, str):
            raise TypeError('Invalid msg_type value, expected str.\n')
        if not ("/" in msg_type or "::" in msg_type):
            raise ValueError('Invalid msg_type format, expected "package/Message" or "package::Message".\n')
        self.__msg_type = msg_type
        
    @field_names.setter
    def field_names(self, field_names: list[str]):
        if not isinstance(field_names, list) or not all(isinstance(item, str) for item in field_names):
            raise TypeError('Invalid field_names value, expected list of str.\n')
        self.__field_names = field_names
        
    @field_types.setter
    def field_types(self, field_types: list[str]):
        if not isinstance(field_types, list) or not all(isinstance(item, str) for item in field_types):
            raise TypeError('Invalid field_types value, expected list of str.\n')
        self.__field_types = field_types
        
    def to_dict(self) -> dict:
        return {
            "topic_name": self.topic_name,
            "msg_type": self.msg_type,
            "field_names": self.field_names,
            "field_types": self.field_types,
        }


class RosTopicMapper(AcadosRosBaseOptions):
    def __init__(self):
        super().__init__()
        self._package_name: str = "ros_mapper"
        self._node_name: str = ""
        self._namespace: str = ""
        self._archtype: str = ArchType.NODE.value
        self._control_loop_executor: str = ControlLoopExec.TOPIC.value
        self.__header_includes: set[str] = set()
        self.__dependencies: set[str] = set()

        self.__input_topic_msg: RosTopicMsg = RosTopicMsg()
        self.__output_topic_msg: RosTopicMsg = RosTopicMsg()
        self.__value_mappings: dict[str, str] = dict()

        self.__ocp_json_file: str = ""
        self.__sim_json_file: str = ""
        self.__uses_auto: bool = False

    @property
    def input_topic_msg(self) -> Optional[RosTopicMsg]:
        return self.__input_topic_msg
    
    @property
    def output_topic_msg(self) -> Optional[RosTopicMsg]:
        return self.__output_topic_msg
    
    @property
    def value_mappings(self) -> Optional[dict[str, str]]:
        return self.__value_mappings

    @property
    def ocp_json_file(self) -> str:
        return self.__ocp_json_file

    @property
    def sim_json_file(self) -> str:
        return self.__sim_json_file

    def finalize(self):
        self.__uses_auto = (self.ocp_json_file and self.sim_json_file)
        
        non_values = [non_val for non_val in [self.input_topic_msg.topic_name, self.input_topic_msg.msg_type, 
                                     self.output_topic_msg.topic_name, self.output_topic_msg.msg_type, 
                                     self.value_mappings] if not non_val]

        if not self.__uses_auto and not non_values:
            input_msg_type = self.input_topic_msg.msg_type
            self.__add_types(input_msg_type)
            self.input_topic_msg.msg_type = self.__correct_msg_type(input_msg_type)

            output_msg_type = self.output_topic_msg.msg_type
            self.__add_types(output_msg_type)
            self.output_topic_msg.msg_type = self.__correct_msg_type(output_msg_type)
            
        else:
            raise ValueError("Either both ocp_json_file and sim_json_file must be set, or input_topic_msg, output_topic_msg, and value_mappings must be set. Currently missing: " + ", ".join(non_values))
            

    def to_dict(self):
        self.finalize()
        return super().to_dict() | {
            "input_topic_msg": self.input_topic_msg.to_dict(),
            "output_topic_msg": self.output_topic_msg.to_dict(),
            "ocp_json_file": self.ocp_json_file,
            "sim_json_file": self.sim_json_file,
            "header_includes": list(self.__header_includes),
            "dependencies": list(self.__dependencies),
        }
        
    def __add_types(self, msg_type: str):
        if "/" in msg_type:
            splitted = msg_type.split("/")
            self.__header_includes.add(f'{splitted[0]}/msg/{self.__camel_to_snake(splitted[-1])}.hpp')
            self.__dependencies.add(splitted[0])
        elif "::" in msg_type:
            splitted = msg_type.split("::")
            self.__header_includes.add(f'{splitted[0]}/msg/{self.__camel_to_snake(splitted[-1])}.hpp')
            self.__dependencies.add(splitted[0])

    @staticmethod
    def __correct_msg_type(msg_type: str) -> str:
        if "/" in msg_type:
            splitted = msg_type.split("/")
        elif "::" in msg_type:
            splitted = msg_type.split("::")
        else:
            raise ValueError(f"Invalid msg_type format: {msg_type}")
        return f'{splitted[0]}::msg::{splitted[-1]}'
    
    
if __name__ == "__main__":
    ros_mapper = RosTopicMapper()
    
    ros_mapper.node_name = "my_mapper_node"
    ros_mapper.control_loop_executor = ControlLoopExec.TOPIC
    ros_mapper.archtype = ArchType.NODE
    
    ros_mapper.input_topic_msg.topic_name = "ocp_state"
    ros_mapper.input_topic_msg.msg_type = "std_msgs/Float32MultiArray"
    
    ros_mapper.output_topic_msg.topic_name = "sim_state"
    ros_mapper.output_topic_msg.msg_type = "std_msgs/Float32MultiArray"
    
    # ros_mapper.ocp_json_file = "ocp_config.json"
    # ros_mapper.sim_json_file = "sim_config.json"
    
    print(ros_mapper.to_dict())