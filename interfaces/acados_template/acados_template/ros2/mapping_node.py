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
# POSSIBILITY OF SUCH DAMAGE.
#

import json
import os
import re

from typing import Optional, Union
from itertools import chain

from .utils import ArchType, ControlLoopExec, AcadosRosBaseOptions
from ..utils import make_object_json_dumpable, render_template
from ..acados_ocp import AcadosOcp
from ..acados_sim import AcadosSim



def _parse_msg_type(msg_type: str) -> tuple[str, str]:
    msg_type = msg_type.strip()
    if "/" in msg_type:
        # "pkg/Type"                     
        pkg, typ = msg_type.split("/", 1)
    elif "::" in msg_type:
        # "pkg::Type" oder "pkg::msg::Type"
        parts = [p for p in msg_type.split("::") if p]
        if len(parts) < 2:
            raise ValueError(f"Invalid msg_type format: {msg_type}")
        pkg, typ = parts[0], parts[-1]
    else:
        raise ValueError(f"Invalid msg_type format: {msg_type}. Valid types are e.g.: 'std_msgs/Header', 'std_msgs::Header' or 'std_msgs::msg::Header'")
    return pkg, typ


def _cpp_msg_type(msg_type: str) -> str:
    pkg, typ = _parse_msg_type(msg_type)
    return f'{pkg}::msg::{typ}'



class RosField:
    def __init__(
        self, 
        name: str, 
        ftype: str, 
        is_array: bool = False, 
        array_size: Optional[int] = None, 
        children: list['RosField'] | None = None
    ):
        if not isinstance(name, str):
            raise TypeError("RosField.name must be str")
        if not isinstance(ftype, str):
            raise TypeError("RosField.ftype must be str")
        if not isinstance(is_array, bool):
            raise TypeError("RosField.is_array must be bool")
        if array_size is not None and not isinstance(array_size, int):
            raise TypeError("RosField.array_size must be int")
        if children is None:
            children = []
        if not isinstance(children, list) or not all(isinstance(ch, RosField) for ch in children):
            raise TypeError("RosField.children must be list[RosField]")
        self.name = name
        self.ftype = ftype
        self.is_array = is_array
        self.array_size = array_size
        self.children = children
        self.cpp_type = ""
        self.needs_storage = False

    def to_dict(self) -> dict:
        self.cpp_type = self.__get_cpp_type(self.ftype)
        return {
            "name": self.name,
            "cpp_type": self.cpp_type,
            "is_array": self.is_array,
            "array_size": self.array_size,
            "needs_storage": self.needs_storage,
            "children": [c.to_dict() for c in self.children],
        }

    def flatten(self, field: Optional['RosField'] = None) -> list['RosField']:
        new_name = self.name if field is None else f"{field.name}.{self.name}"
        field_copy = RosField(new_name, self.ftype, self.is_array, self.array_size, self.children)
        
        if not self.children:
            return [field_copy]
        
        out: list['RosField'] = []
        for c in self.children:
            out.extend(c.flatten(field_copy))
        return out
    
    @staticmethod
    def __get_cpp_type(ftype: str):
        match ftype:
            case "float64":
                return "double"
            case "float32":
                return "float"
            case "int8":
                return "int8_t"
            case "int16":
                return "int16_t"
            case "int32":
                return "int32_t"
        
        if "/" in ftype or "::" in ftype:
            return _cpp_msg_type(ftype)
        
        return ftype
    

class RosTopicMsg:
    def __init__(self):
        self.__topic_name: str = ""
        self.__msg_type: str = ""
        self.__field_tree: list[RosField] = list()
        self._flat_field_tree: list[RosField] = list()

    @property
    def topic_name(self) -> str:
        return self.__topic_name

    @property
    def msg_type(self) -> str:
        return self.__msg_type
    
    @property
    def field_tree(self) -> list[RosField]:
        return self.__field_tree

    @topic_name.setter
    def topic_name(self, value: str):
        if not isinstance(value, str):
            raise TypeError('Invalid topic_name value, expected str.\n')
        self.__topic_name = value
        
    @msg_type.setter
    def msg_type(self, value: str):
        if not isinstance(value, str):
            raise TypeError('Invalid msg_type value, expected str.\n')
        if not ("/" in value or "::" in value):
            raise ValueError('Invalid msg_type format, expected "package/Message" or "package::Message".\n')
        self.__msg_type = value

    @field_tree.setter
    def field_tree(self, value: list[RosField]):
        if not isinstance(value, list) or not all(isinstance(item, RosField) for item in value):
            raise TypeError('Invalid field_tree value, expected list of RosField.\n')
        self.__field_tree = value
        
    def flatten_field_tree(self):
        self._flat_field_tree = list(chain.from_iterable(f.flatten() for f in self.field_tree))

    def to_dict(self) -> dict:
        return {
            "topic_name": self.topic_name,
            "msg_type": self.msg_type,
            # "field_tree": [field.to_dict() for field in self.field_tree],
            "flat_field_tree": [field.to_dict() for field in self._flat_field_tree]
        }


class RosTopicMsgOutput(RosTopicMsg):
    _NOT_IMPLEMENTED_ARCHTYPES: set[ArchType] = {
        ArchType.LIFECYCLE_NODE,
        ArchType.ROS2_CONTROLLER,
        ArchType.NAV2_CONTROLLER}
    _NOT_IMPLEMENTED_EXECUTORS: set[ControlLoopExec] = {
        ControlLoopExec.TIMER,
        ControlLoopExec.SRV,
        ControlLoopExec.ACTION}
    
    def __init__(self):
        super().__init__()
        self.__mapping: list[dict] = list()
        self.__exec_topic: str = ""
        self._needs_publish_lock: bool = False
        
    @property
    def mapping(self) -> list[dict]:
        return self.__mapping

    @property
    def exec_topic(self) -> str:
        return self.__exec_topic
            
    @mapping.setter
    def mapping(self, value: list[tuple[str, str]]):
        if not isinstance(value, list) or not all(isinstance(item, tuple) and len(item) == 2 and all(isinstance(i, str) for i in item) for item in value):
            raise TypeError('Invalid mapping value, expected list of tuples (str, str).\n')
        self.__mapping = [self.__parse_mapping_pair(src, dest) for src, dest in value]

    @exec_topic.setter
    def exec_topic(self, value: str):
        if not isinstance(value, str):
            raise TypeError('Invalid exec_topic value, expected str.\n')
        self.__exec_topic = value

    def to_dict(self) -> dict:
        return super().to_dict() | {
            "mapping": self.mapping,
            "exec_topic": self.exec_topic,
            "needs_publish_lock": self._needs_publish_lock}
        
    @staticmethod
    def __parse_mapping_string(map_str: str) -> dict:
        """Parses a string like 'field.name[index]' into a structured dict."""
        match = re.match(r"^(.*?)(?:\[(\d+)\])?$", map_str)
        if not match:
            raise ValueError(f"Invalid mapping format: '{map_str}'")
        
        base_name, index_str = match.groups()
        index = int(index_str) if index_str is not None else None
        
        return {"full_name": map_str, "base_name": base_name, "index": index}

    def __parse_mapping_pair(self, source_str: str, dest_str: str) -> dict:
        """Parses a source-destination pair into a structured dict."""
        source_parts = self.__parse_mapping_string(source_str)
        dest_parts = self.__parse_mapping_string(dest_str)

        # Split source into topic and field
        base_split = source_parts['base_name'].split('.', 1)
        if len(base_split) != 2:
            raise ValueError(f"Source mapping '{source_str}' must be in 'topic_name.field_name' format.")
        
        return {
            "source": {
                "full": source_str,
                "topic": base_split[0],
                "field": base_split[1],
                "index": source_parts["index"]
            },
            "dest": {
                "full": dest_str,
                "field": dest_parts["base_name"],
                "index": dest_parts["index"]
            }
        }



class RosTopicMapper(AcadosRosBaseOptions):
    def __init__(self):
        super().__init__()
        self.package_name: str = "ros_mapper"
        self.node_name: str = ""
        self.namespace: str = ""
        self.archtype: str = ArchType.NODE.value
        self.control_loop_executor: str = ControlLoopExec.TOPIC.value
        self.__header_includes: set[str] = set()
        self.__dependencies: set[str] = set()

        self.__in_msgs: list[RosTopicMsg] = []
        self.__out_msgs: list[RosTopicMsgOutput] = []

        self.__ocp_instance: Optional[AcadosOcp] = None
        self.__sim_instance:  Optional[AcadosSim] = None
        self.__mapper_json_file: str = "ros_mapper.json"

    @property
    def in_msgs(self) -> list[RosTopicMsg]:
        return self.__in_msgs
    
    @property
    def out_msgs(self) -> list[RosTopicMsgOutput]:
        return self.__out_msgs
    
    @property
    def ocp_instance(self) -> Optional[AcadosOcp]:
        return self.__ocp_instance

    @property
    def sim_instance(self) -> Optional[AcadosSim]:
        return self.__sim_instance

    @property
    def mapper_json_file(self) -> str:
        return self.__mapper_json_file

    @in_msgs.setter
    def in_msgs(self, value: list[RosTopicMsg]):
        if not isinstance(value, list) or not all(isinstance(item, RosTopicMsg) for item in value):
            raise TypeError('Invalid in_msg value, expected list of RosTopicMsg.\n')
        self.__in_msgs = value

    @out_msgs.setter
    def out_msgs(self, value: list[RosTopicMsgOutput]):
        if not isinstance(value, list ) or not all(isinstance(item, RosTopicMsgOutput) for item in value):
            raise TypeError('Invalid out_msg value, expected list of RosTopicMsgOutput.\n')
        self.__out_msgs = value

    @ocp_instance.setter
    def ocp_instance(self, value: Optional[AcadosOcp]):
        if value is not None and not isinstance(value, AcadosOcp):
            raise TypeError('Invalid ocp_instance value, expected AcadosOcp.\n')
        self.__ocp_instance = value

    @sim_instance.setter
    def sim_instance(self, value: Optional[AcadosSim]):
        if value is not None and not isinstance(value, AcadosSim):
            raise TypeError('Invalid sim_instance value, expected AcadosSim.\n')
        self.__sim_instance = value

    @mapper_json_file.setter
    def mapper_json_file(self, value: str):
        if not isinstance(value, str):
            raise TypeError('Invalid mapper_json_file value, expected str.\n')
        self.__mapper_json_file = value


    def check_consistency(self):
        in_topic_index = {m.topic_name: m for m in self.in_msgs}
        
        for out_msg in self.out_msgs:
            out_field_index = {f.name: f for f in out_msg._flat_field_tree}
            
            for mapping in out_msg.mapping:
                source = mapping['source']
                dest = mapping['dest']

                # Check source topic
                in_msg = in_topic_index.get(source['topic'])
                if in_msg is None:
                    raise ValueError(f"Source topic '{source['topic']}' in mapping '{source['full']}' not found in input messages.")
                
                in_field_index = {f.name: f for f in in_msg._flat_field_tree}

                # Check source field
                source_field_obj = in_field_index.get(source['field'])
                if source_field_obj is None:
                    raise ValueError(f"Source field '{source['field']}' not found in topic '{source['topic']}'. Available: {list(in_field_index.keys())}")

                # Check source indexing
                if source['index'] is not None:
                    if not source_field_obj.is_array:
                        raise ValueError(f"Source field '{source['field']}' is not an array, but is being indexed in '{source['full']}'.")
                    if source_field_obj.array_size > 0 and source['index'] >= source_field_obj.array_size:
                        raise ValueError(f"Source index in '{source['full']}' is out of bounds. Size is {source_field_obj.array_size}.")

                # Check destination field
                dest_field_obj = out_field_index.get(dest['field'])
                if dest_field_obj is None:
                    raise ValueError(f"Destination field '{dest['field']}' not found in output message '{out_msg.topic_name}'. Available: {list(out_field_index.keys())}")

                # Check destination indexing
                if dest['index'] is not None:
                    if not dest_field_obj.is_array:
                        raise ValueError(f"Destination field '{dest['field']}' is not an array, but is being indexed in '{dest['full']}'.")
                    if dest_field_obj.array_size > 0 and dest['index'] >= dest_field_obj.array_size:
                        raise ValueError(f"Destination index in '{dest['full']}' is out of bounds. Size is {dest_field_obj.array_size}.")
                    

    def finalize(self):
        self.mapper_json_file = os.path.join(self.generated_code_dir, self.mapper_json_file)
        
        # map ocp to sim
        if self.ocp_instance and self.sim_instance:
            self.map_ocp_to_sim()
            
        # topic normalize + flatten
        for msg in self.in_msgs + self.out_msgs:
            if msg.topic_name.startswith("/"):
                msg.topic_name = msg.topic_name[1:]
            msg.flatten_field_tree()
            self.__add_types(msg.msg_type)
            msg.msg_type = _cpp_msg_type(msg.msg_type)
            
        # execution topic normalize
        for om in self.out_msgs:
            if om.exec_topic.startswith("/"):
                om.exec_topic = om.exec_topic[1:]
        
        self.check_none_values()
        self.check_consistency()
        
        # compute if field needs_storage 
        for in_msg in self.in_msgs:
            for in_field in in_msg._flat_field_tree:
                full_field_name = f"{in_msg.topic_name}.{in_field.name}"
                for out_msg in self.out_msgs:
                    is_used = any(
                        f"{mapping['source']['topic']}.{mapping['source']['field']}" == full_field_name
                        for mapping in out_msg.mapping
                    )
                    if is_used and out_msg.exec_topic != in_msg.topic_name:
                        in_field.needs_storage = True
                        break
        
        # compute if out msg needs lock 
        for out_msg in self.out_msgs:
            out_msg._needs_publish_lock = any(
                m["source"]["topic"] != out_msg.exec_topic
                for m in out_msg.mapping
            )


    def to_dict(self):
        return super().to_dict() | {
            "in_msgs": [msg.to_dict() for msg in self.in_msgs],
            "out_msgs": [msg.to_dict() for msg in self.out_msgs],
            "ocp_json_file": self.ocp_instance,
            "sim_json_file": self.sim_instance,
            "mapper_json_file": self.mapper_json_file,
            "header_includes": list(self.__header_includes),
            "dependencies": list(self.__dependencies),
        }


    def dump_to_json(self) -> None:
        self.finalize()
        os.makedirs(os.path.dirname(self.mapper_json_file), exist_ok=True)
        with open(self.mapper_json_file, 'w') as f:
            json.dump(self.to_dict(), f, default=make_object_json_dumpable, indent=4, sort_keys=True)
            
    
    def render_templates(self):
        if not os.path.exists(self.mapper_json_file):
            raise FileNotFoundError(f"{self.mapper_json_file} not found!")

        template_list = self._get_ros_template_list()

        # Render templates
        for tup in template_list:
            output_dir = self.generated_code_dir if len(tup) <= 2 else tup[2]
            render_template(tup[0], tup[1], output_dir, self.mapper_json_file)
            
    
    def check_none_values(self):
        non_values: list[str] = []
        
        if not self.in_msgs:
            non_values.append(f"input")
        if not self.out_msgs:
            non_values.append(f"output")
        
        for i, msg in enumerate(self.in_msgs):
            if not msg.topic_name:
                non_values.append(f"input[{i}].topic_name")
            if not msg.msg_type:
                non_values.append(f"input[{i}].msg_type")
            if not msg._flat_field_tree:
                non_values.append(f"input[{i}].flat_field_tree")

        for i, msg in enumerate(self.out_msgs):
            if not msg.topic_name:
                non_values.append(f"output[{i}].topic_name")
            if not msg.msg_type:
                non_values.append(f"output[{i}].msg_type")
            if not msg._flat_field_tree:
                non_values.append(f"output[{i}].flat_field_tree")
            if msg.mapping is None or len(msg.mapping) == 0:
                non_values.append(f"output[{i}].mapping")
            if not msg.exec_topic:
                non_values.append(f"output[{i}].exec_topic")

        if non_values:
            raise ValueError("The fields 'in_msgs' and 'out_msgs' must be non-empty. Currently missing: " + ", ".join(non_values))  
    
    
    def map_ocp_to_sim(self):
        """
        Automatically map the states and inputs from the OCP to the simulator,
        via the given json's.
        """
        pass
    
    
    def _get_ros_template_list(self) -> list:
        template_list = []

        # --- Simulator Package --- 
        ros_pkg_dir = os.path.join('ros2', 'ros_mapper_templates')
        package_dir = os.path.join(self.generated_code_dir, self.package_name)
        template_file = os.path.join(ros_pkg_dir, 'README.in.md')
        template_list.append((template_file, 'README.md', package_dir))
        template_file = os.path.join(ros_pkg_dir, 'CMakeLists.in.txt')
        template_list.append((template_file, 'CMakeLists.txt', package_dir))
        template_file = os.path.join(ros_pkg_dir, 'package.in.xml')
        template_list.append((template_file, 'package.xml', package_dir))

        # # Header
        include_dir = os.path.join(package_dir, 'include', self.package_name)
        template_file = os.path.join(ros_pkg_dir, 'utils.in.hpp')
        template_list.append((template_file, 'utils.hpp', include_dir))
        template_file = os.path.join(ros_pkg_dir, 'node.in.h')
        template_list.append((template_file, 'node.h', include_dir))

        # Source
        src_dir = os.path.join(package_dir, 'src')
        template_file = os.path.join(ros_pkg_dir, 'node.in.cpp')
        template_list.append((template_file, 'node.cpp', src_dir))
        
        # # Test
        # test_dir = os.path.join(package_dir, 'test')
        # template_file = os.path.join(ros_pkg_dir, 'test.launch.in.py')
        # template_list.append((template_file, f'test_{self.package_name}.launch.py', test_dir))
        return template_list
        
        
    def __add_types(self, msg_type: str):
        pkg, typ = _parse_msg_type(msg_type)
        self.__header_includes.add(f"{pkg}/msg/{self.camel_to_snake(typ)}.hpp")
        self.__dependencies.add(pkg)
        
    

    
if __name__ == "__main__":
    from pathlib import Path
    ros_mapper = RosTopicMapper()
    
    # ros_mapper.package_name = "my_mapper_node"
    ros_mapper.control_loop_executor = ControlLoopExec.TOPIC
    ros_mapper.archtype = ArchType.NODE
    
    in_msg = RosTopicMsg()
    in_msg.topic_name = "next_state"
    in_msg.msg_type = "sim_package_interface/State"
    in_msg.field_tree += [
        RosField(name="header", ftype="std_msgs/Header"),
        RosField(name="x", ftype="float64", is_array=True, array_size=5),
        RosField(name="status", ftype="int8")]
    ros_mapper.in_msgs += [in_msg]
    
    
    out_msg = RosTopicMsgOutput()
    out_msg.topic_name = "state"
    out_msg.msg_type = "ocp_package_interface/State"
    out_msg.field_tree += [
        RosField(name="header", ftype="std_msgs/Header"),
        RosField(name="x", ftype="float64", is_array=True, array_size=4)]
    out_msg.mapping = [
        ("next_state.x[2]", "x[0]"),
        ("next_state.x[3]", "x[1]"),
        ("next_state.x[4]", "x[2]"),
        ("control_input.u[0]","x[3]")]
    out_msg.exec_topic = "next_state"
    ros_mapper.out_msgs += [out_msg]
    
    
    in_msg = RosTopicMsg()
    in_msg.topic_name = "control_input"
    in_msg.msg_type = "ocp_package_interface/ControlInput"
    in_msg.field_tree += [
        RosField(name="header", ftype="std_msgs/Header"),
        RosField(name="u", ftype="float64", is_array=True, array_size=2),
        RosField(name="status", ftype="int8"),
        RosField(name="test_val", ftype="Pose", children=[
            RosField(name="x", ftype="float64"),
            RosField(name="y", ftype="float64"),
            RosField(name="z", ftype="float64")])]
    ros_mapper.in_msgs += [in_msg]
    
    out_msg = RosTopicMsgOutput()
    out_msg.topic_name = "sim_control"
    out_msg.msg_type = "sim_package_interface/State"
    out_msg.field_tree += [
        RosField(name="header", ftype="std_msgs/Header"),
        RosField(name="u", ftype="float64", is_array=True, array_size=2)]
    out_msg.mapping = [
        ("control_input.u", "u")]
    out_msg.exec_topic = "control_input"
    ros_mapper.out_msgs += [out_msg]
    
    ros_mapper.generated_code_dir = str(Path(__file__).resolve().parent / "generated")
    
    from pprint import pprint
    pprint(ros_mapper.to_dict())
    
    ros_mapper.dump_to_json()
    ros_mapper.render_templates()