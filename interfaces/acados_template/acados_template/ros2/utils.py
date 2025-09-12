import re

from enum import Enum
from typing import Union


class ControlLoopExec(str, Enum):
    TOPIC = "topic"
    TIMER = "timer"
    SRV = "srv"
    ACTION = "action"

class ArchType(str, Enum):
    NODE = "node"
    LIFECYCLE_NODE = "lifecycle_node"
    ROS2_CONTROLLER = "ros2_controller"
    NAV2_CONTROLLER = "nav2_controller"


class AcadosRosBaseOptions:
    def __init__(self):
        self.__package_name: str = "acados_base"
        self.__node_name: str = "acados_base_node"
        self.__namespace: str = ""
        self.__archtype: str = ArchType.NODE.value
        self.__control_loop_executor: str = ControlLoopExec.TIMER.value

    @property
    def package_name(self) -> str:
        return self.__package_name

    @property
    def node_name(self) -> str:
        return self.__node_name

    @property
    def namespace(self) -> str:
        return self.__namespace

    @property
    def archtype(self) -> str:
        return self.__archtype

    @property
    def control_loop_executor(self) -> str:
        return self.__control_loop_executor

    @package_name.setter
    def package_name(self, package_name: str):
        if not isinstance(package_name, str):
            raise TypeError('Invalid package_name value, expected str.\n')
        self.__package_name = package_name

    @node_name.setter
    def node_name(self, node_name: str):
        if not isinstance(node_name, str):
            raise TypeError('Invalid node_name value, expected str.\n')
        self.__node_name = node_name

    @namespace.setter
    def namespace(self, namespace: str):
        if not isinstance(namespace, str):
            raise TypeError('Invalid namespace value, expected str.\n')
        self.__namespace = namespace

    @archtype.setter
    def archtype(self, node_archtype: Union[ArchType, str]):
        try:
            if isinstance(node_archtype, ArchType):
                archtype_enum = node_archtype
            elif isinstance(node_archtype, str):
                archtype_enum = ArchType(node_archtype)
            else:
                raise TypeError()
        except (ValueError, TypeError):
            valid_types = [e.value for e in ArchType]
            raise TypeError(f"Invalid node_archtype. Expected one of {valid_types} or a ArchType enum member.")

        not_implemented = [
            ArchType.ROS2_CONTROLLER,
            ArchType.NAV2_CONTROLLER,
            ArchType.LIFECYCLE_NODE
        ]
        if archtype_enum in not_implemented:
            raise NotImplementedError(f"Archtype '{archtype_enum.value}' is not implemented yet.")
        self.__archtype = archtype_enum.value

    @control_loop_executor.setter
    def control_loop_executor(self, control_loop_executor: Union[ControlLoopExec, str]) -> None:
        try:
            if isinstance(control_loop_executor, ControlLoopExec):
                control_loop_executor_enum = control_loop_executor
            elif isinstance(control_loop_executor, str):
                control_loop_executor_enum = ControlLoopExec(control_loop_executor)
            else:
                raise TypeError()
        except (ValueError, TypeError):
            valid_types = [e.value for e in ControlLoopExec]
            raise TypeError(f"Invalid control_loop_executor. Expected one of {valid_types} or a ControlLoopExec enum member.")
        
        not_implemented = [
            ControlLoopExec.SRV,
            ControlLoopExec.ACTION
        ]
        if control_loop_executor_enum in not_implemented:
            raise NotImplementedError(f"Control loop executor '{control_loop_executor_enum.value}' is not implemented yet.")
        self.__control_loop_executor = control_loop_executor_enum.value

    def to_dict(self) -> dict:
        if self.node_name == "":
            self.node_name = self.package_name + "_node"
        package_name_snake  = self.camel_to_snake(self.package_name)
        node_name_snake     = self.camel_to_snake(self.node_name)
        namespace_snake     = self.camel_to_snake(self.namespace)
        return {
            "package_name":             package_name_snake,
            "node_name":                node_name_snake,
            "namespace":                namespace_snake,
            "archtype":                 self.archtype,
            "control_loop_executor":    self.control_loop_executor
        }

    @staticmethod
    def camel_to_snake(name: str) -> str:
        s1 = re.sub(r'(.)([A-Z][a-z]+)', r'\1_\2', name)
        return re.sub(r'([a-z0-9])([A-Z])', r'\1_\2', s1).lower()
