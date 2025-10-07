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

from typing import Optional, Union, TYPE_CHECKING, Literal, List, Tuple
from itertools import chain
from casadi import SX, MX

from .utils import ArchType, AcadosRosBaseOptions
from ..utils import make_object_json_dumpable, render_template

# Avoid circular imports at runtime; only import these for type checking
if TYPE_CHECKING:
    from ..acados_ocp import AcadosOcp
    from ..acados_sim import AcadosSim
    from ..acados_ocp_solver import AcadosOcpSolver
    from ..acados_sim_solver import AcadosSimSolver



def _parse_msg_type(msg_type: str) -> Tuple[str, str]:
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
        children: Optional[List['RosField']] = None
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
            raise TypeError("RosField.children must be List[RosField]")
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

    def flatten(self, field: Optional['RosField'] = None) -> List['RosField']:
        new_name = self.name if field is None else f"{field.name}.{self.name}"
        field_copy = RosField(new_name, self.ftype, self.is_array, self.array_size, self.children)

        if not self.children:
            return [field_copy]

        out: List['RosField'] = []
        for c in self.children:
            out.extend(c.flatten(field_copy))
        return out

    @staticmethod
    def __get_cpp_type(ftype: str):
        if ftype == "float64":
            return "double"
        elif ftype == "float32":
            return "float"
        elif ftype == "int8":
            return "int8_t"
        elif ftype == "int16":
            return "int16_t"
        elif ftype == "int32":
            return "int32_t"

        if "/" in ftype or "::" in ftype:
            return _cpp_msg_type(ftype)

        return ftype


class RosTopicMsg:
    def __init__(self):
        self.__topic_name: str = ""
        self.__msg_type: str = ""
        self.__field_tree: List[RosField] = list()
        self._flat_field_tree: List[RosField] = list()

    @property
    def topic_name(self) -> str:
        return self.__topic_name

    @property
    def msg_type(self) -> str:
        return self.__msg_type

    @property
    def field_tree(self) -> List[RosField]:
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
    def field_tree(self, value: List[RosField]):
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
    def __init__(self):
        super().__init__()
        self.__mapping: List[dict] = list()
        self.__exec_topic: str = ""
        self._needs_publish_lock: bool = False

    @property
    def mapping(self) -> List[dict]:
        return self.__mapping

    @property
    def exec_topic(self) -> str:
        return self.__exec_topic

    @mapping.setter
    def mapping(self, value: List[Tuple[str, str]]):
        if not isinstance(value, list) or not all(isinstance(item, tuple) and len(item) == 2 and all(isinstance(i, str) for i in item) for item in value):
            raise TypeError('Invalid mapping value, expected list of tuples (str, str).\n')
        self.__mapping = [self.__parse_mapping_pair(src, dest) for src, dest in value]

    @exec_topic.setter
    def exec_topic(self, value: str):
        if not isinstance(value, str):
            raise TypeError('Invalid exec_topic value, expected str.\n')
        self.__exec_topic = value

    def to_dict(self) -> dict:
        d = super().to_dict()
        d.update({
            "mapping": self.mapping,
            "exec_topic": self.exec_topic,
            "needs_publish_lock": self._needs_publish_lock
        })
        return d

    @classmethod
    def from_msg(
        cls,
        base_msg: RosTopicMsg
    ) -> 'RosTopicMsgOutput':
        """
        Create a RosTopicMsgOutput instance from a base RosTopicMsg template.
        """
        output_msg = cls()
        output_msg.msg_type = base_msg.msg_type
        output_msg.field_tree = base_msg.field_tree
        output_msg.topic_name = base_msg.topic_name
        return output_msg

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
        self.__header_includes: set[str] = set()
        self.__dependencies: set[str] = set()

        self.__in_msgs: List[RosTopicMsg] = []
        self.__out_msgs: List[RosTopicMsgOutput] = []

        self.__ocp_json_file: str = ""
        self.__sim_json_file: str = ""
        self.__mapper_json_file = "ros_mapper.json"

    @property
    def in_msgs(self) -> List[RosTopicMsg]:
        return self.__in_msgs

    @property
    def out_msgs(self) -> List[RosTopicMsgOutput]:
        return self.__out_msgs

    @property
    def ocp_json_file(self) -> str:
        return self.__ocp_json_file

    @property
    def sim_json_file(self) -> str:
        return self.__sim_json_file

    @property
    def mapper_json_file(self) -> str:
        return self.__mapper_json_file

    @in_msgs.setter
    def in_msgs(self, value: List[RosTopicMsg]):
        if not isinstance(value, list) or not all(isinstance(item, RosTopicMsg) for item in value):
            raise TypeError('Invalid in_msg value, expected list of RosTopicMsg.\n')
        self.__in_msgs = value

    @out_msgs.setter
    def out_msgs(self, value: List[RosTopicMsgOutput]):
        if not isinstance(value, list ) or not all(isinstance(item, RosTopicMsgOutput) for item in value):
            raise TypeError('Invalid out_msg value, expected list of RosTopicMsgOutput.\n')
        self.__out_msgs = value

    @ocp_json_file.setter
    def ocp_json_file(self, value: str):
        if not isinstance(value, str):
            raise TypeError('Invalid ocp_json_file value, expected str.\n')
        self.__ocp_json_file = value

    @sim_json_file.setter
    def sim_json_file(self, value: str):
        if not isinstance(value, str):
            raise TypeError('Invalid sim_json_file value, expected str.\n')
        self.__sim_json_file = value

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
            "ocp_json_file": self.ocp_json_file,
            "sim_json_file": self.sim_json_file,
            "mapper_json_file": self.mapper_json_file,
            "header_includes": list(self.__header_includes),
            "dependencies": list(self.__dependencies),
        }


    def dump_to_json(self) -> None:
        dir_name = os.path.dirname(self.mapper_json_file)
        if not dir_name:
            self.mapper_json_file = os.path.join(self.generated_code_dir, self.mapper_json_file)

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
            template_glob = None if len(tup) <= 3 else tup[3]
            render_template(tup[0], tup[1], output_dir, self.mapper_json_file, template_glob=template_glob)


    def check_none_values(self):
        non_values: List[str] = []

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


    def generate(self):
        self.finalize()
        self.dump_to_json()
        self.render_templates()


    @classmethod
    def from_instances(cls, ocp_solver: 'AcadosOcpSolver', sim_solver: 'AcadosSimSolver'):
        """
        Build a :py:class:`RosTopicMapper` by wiring OCP -> SIM controls (``u``) and
        SIM -> OCP states (``x``) using the default interface message types.

        The mapper creates two input messages (SIM ``State``, OCP ``ControlInput``)
        and two output messages (OCP ``State``, SIM ``ControlInput``). It then derives
        field mappings based on symbolic element names of the respective CasADi
        vectors and assigns execution topics so that each output is published when
        its corresponding input arrives.

        :param ocp_solver: OCP solver instance whose OCP and model define the OCP side.
        :type ocp_solver: acados_template.acados_ocp_solver.AcadosOcpSolver
        :param sim_solver: Simulator solver instance whose simulator and model define the SIM side.
        :type sim_solver: acados_template.acados_sim_solver.AcadosSimSolver
        :returns: A configured topic mapper with populated inputs, outputs, mappings,
                  and references to the OCP/SIM JSON files.
        :rtype: RosTopicMapper

        :raises ValueError: If required ROS options (e.g. ``package_name``, topics) are
            missing in either solver, or if message field trees are incomplete at finalize-time.
        """
        ocp = ocp_solver.acados_ocp
        sim = sim_solver.acados_sim
        obj = cls()

        # State and Control names
        ocp_x_names = _elem_names(ocp.model.x)
        sim_x_names = _elem_names(sim.model.x)
        ocp_u_names = _elem_names(ocp.model.u)
        sim_u_names = _elem_names(sim.model.u)

        # --- Build input messages ---
        in_sim_state = build_default_state(sim, direction_out=False)
        in_ocp_ctrl  = build_default_control(ocp, direction_out=False)
        obj.in_msgs = [in_sim_state, in_ocp_ctrl]

        # --- Build output messages ---
        out_ocp_state = build_default_state(ocp, direction_out=True)
        out_sim_ctrl  = build_default_control(sim, direction_out=True)

        # Map SIM x -> OCP x
        x_map_pairs = _compute_mapping(
            src_topic=in_sim_state.topic_name,
            src_name="x",
            src_labels=sim_x_names,
            dst_name="x",
            dst_labels=ocp_x_names,
        )
        out_ocp_state.mapping = x_map_pairs
        out_ocp_state.exec_topic = in_sim_state.topic_name

        # Map OCP u -> SIM u
        u_map_pairs = _compute_mapping(
            src_topic=in_ocp_ctrl.topic_name,
            src_name="u",
            src_labels=ocp_u_names,
            dst_name="u",
            dst_labels=sim_u_names,
        )
        out_sim_ctrl.mapping = u_map_pairs
        out_sim_ctrl.exec_topic = in_ocp_ctrl.topic_name

        obj.out_msgs = [out_ocp_state, out_sim_ctrl]
        obj.ocp_json_file = ocp.json_file
        obj.sim_json_file = sim.json_file
        return obj


    def _get_ros_template_list(self) -> list:
        template_list = []

        acados_template_path = os.path.dirname(os.path.dirname(__file__))
        ros_template_glob = os.path.join(acados_template_path, 'ros2_templates', '**', '*')

        # --- Simulator Package ---
        ros_pkg_dir = os.path.join('ros_mapper_templates')
        package_dir = os.path.join(self.generated_code_dir, self.package_name)
        template_file = os.path.join(ros_pkg_dir, 'README.in.md')
        template_list.append((template_file, 'README.md', package_dir, ros_template_glob))
        template_file = os.path.join(ros_pkg_dir, 'CMakeLists.in.txt')
        template_list.append((template_file, 'CMakeLists.txt', package_dir, ros_template_glob))
        template_file = os.path.join(ros_pkg_dir, 'package.in.xml')
        template_list.append((template_file, 'package.xml', package_dir, ros_template_glob))

        # # Header
        include_dir = os.path.join(package_dir, 'include', self.package_name)
        template_file = os.path.join(ros_pkg_dir, 'utils.in.hpp')
        template_list.append((template_file, 'utils.hpp', include_dir, ros_template_glob))
        template_file = os.path.join(ros_pkg_dir, 'node.in.h')
        template_list.append((template_file, 'node.h', include_dir, ros_template_glob))

        # Source
        src_dir = os.path.join(package_dir, 'src')
        template_file = os.path.join(ros_pkg_dir, 'node.in.cpp')
        template_list.append((template_file, 'node.cpp', src_dir, ros_template_glob))

        # Test
        test_dir = os.path.join(package_dir, 'test')
        template_file = os.path.join(ros_pkg_dir, 'test.launch.in.py')
        template_list.append((template_file, f'test_{self.package_name}.launch.py', test_dir, ros_template_glob))
        return template_list


    def __add_types(self, msg_type: str):
        pkg, typ = _parse_msg_type(msg_type)
        self.__header_includes.add(f"{pkg}/msg/{self.camel_to_snake(typ)}.hpp")
        self.__dependencies.add(pkg)



# --- From Instances Helpers ---
def _size(sym: Union[SX, MX]) -> int:
    return int(sym.numel())


def _elem_names(sym: Union[SX, MX]) -> List[str]:
    n = _size(sym)
    return [str(sym[i].name()) for i in range(n)]


def _compute_mapping(
        src_topic: str,
        src_name: str,
        src_labels: List[str],
        dst_name: str,
        dst_labels: List[str]
) -> List[Tuple[str, str]]:
    """Return mapping pairs [(f"{src_topic}.{src_name}[i]","{dst_name}[j]")]."""
    if src_labels == dst_labels:
        return [(f"{src_topic}.{src_name}", f"{dst_name}")]

    pairs: List[Tuple[str, str]] = []

    # Label-based match
    for j, lbl in enumerate(dst_labels):
        i = src_labels.get(lbl.lower())
        if i is not None:
            pairs.append((f"{src_topic}.{src_name}[{i}]", f"{dst_name}[{j}]"))

    return pairs


def build_default_state(
        solver_instance: Union['AcadosOcp', 'AcadosSim'],
        direction_out: bool = False
) -> Union[RosTopicMsg, RosTopicMsgOutput]:
    """
    Creates a standard ocp- or sim-interface state message of an
    in- or outgoing message, dependent on the setting.
    """
    if not (hasattr(solver_instance, "ros_opts") and getattr(solver_instance, "ros_opts") is not None):
        raise ValueError(f"Field 'ros_opts' is not set in the solver {solver_instance.__class__.__name__}")

    m = RosTopicMsgOutput() if direction_out else RosTopicMsg()
    m.topic_name = solver_instance.ros_opts.state_topic
    m.msg_type = f"{solver_instance.ros_opts.package_name}_interface/State"
    m.field_tree = [
        RosField(name="header", ftype="std_msgs/Header"),
        RosField(name="x", ftype="float64", is_array=True, array_size=solver_instance.model.x.rows()),
    ]
    if direction_out:
        m.field_tree.append(RosField(name="status", ftype="int8"))
    m.flatten_field_tree()
    return m


def build_default_control(
        solver_instance: Union['AcadosOcp', 'AcadosSim'],
        direction_out: bool = False
) -> Union[RosTopicMsg, RosTopicMsgOutput]:
    """
    Creates a standard ocp- or sim-interface control message of an
    in- or outgoing message, dependent on the setting.
    """
    if not (hasattr(solver_instance, "ros_opts") and getattr(solver_instance, "ros_opts") is not None):
        raise ValueError(f"Field 'ros_opts' is not set in the solver {solver_instance.__class__.__name__}")

    m = RosTopicMsgOutput() if direction_out else RosTopicMsg()
    m.topic_name = solver_instance.ros_opts.control_topic
    m.msg_type = f"{solver_instance.ros_opts.package_name}_interface/ControlInput"
    m.field_tree = [
        RosField(name="header", ftype="std_msgs/Header"),
        RosField(name="u", ftype="float64", is_array=True, array_size=solver_instance.model.u.rows()),
    ]
    if not direction_out:
        m.field_tree.append(RosField(name="status", ftype="int8"))
    m.flatten_field_tree()
    return m
