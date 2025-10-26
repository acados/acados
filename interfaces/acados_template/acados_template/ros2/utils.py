import re

from enum import Enum
from typing import Union, Set


class ArchType(str, Enum):
    NODE = "node"
    LIFECYCLE_NODE = "lifecycle_node"
    ROS2_CONTROLLER = "ros2_controller"
    NAV2_CONTROLLER = "nav2_controller"


class AcadosRosBaseOptions:
    _NOT_IMPLEMENTED_ARCHTYPES: Set[ArchType] = {
        ArchType.LIFECYCLE_NODE,
        ArchType.ROS2_CONTROLLER,
        ArchType.NAV2_CONTROLLER}

    def __init__(self):
        self.__package_name: str = "acados_base"
        self.__node_name: str = ""
        self.__namespace: str = ""
        self.__generated_code_dir: str = "ros_generated_code"
        self.__archtype: str = ArchType.NODE.value

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
    def generated_code_dir(self) -> str:
        return self.__generated_code_dir

    @property
    def archtype(self) -> str:
        return self.__archtype

    @package_name.setter
    def package_name(self, value: str):
        if not isinstance(value, str):
            raise TypeError('Invalid package_name value, expected str.\n')
        self.__package_name = value

    @node_name.setter
    def node_name(self, value: str):
        if not isinstance(value, str):
            raise TypeError('Invalid node_name value, expected str.\n')
        self.__node_name = value

    @generated_code_dir.setter
    def generated_code_dir(self, value: str):
        if not isinstance(value, str):
            raise TypeError('Invalid generated_code_dir value, expected str.\n')
        self.__generated_code_dir = value

    @namespace.setter
    def namespace(self, value: str):
        if not isinstance(value, str):
            raise TypeError('Invalid namespace value, expected str.\n')
        self.__namespace = value

    @archtype.setter
    def archtype(self, value: Union[ArchType, str]):
        try:
            if isinstance(value, ArchType):
                archtype_enum = value
            elif isinstance(value, str):
                archtype_enum = ArchType(value)
            else:
                raise TypeError()
        except (ValueError, TypeError):
            valid_types = [e.value for e in ArchType]
            raise TypeError(f"Invalid node_archtype. Expected one of {valid_types} or a ArchType enum member.")

        if archtype_enum in type(self)._NOT_IMPLEMENTED_ARCHTYPES:
            raise NotImplementedError(f"Archtype '{archtype_enum.value}' is not implemented for {type(self).__name__} yet.")
        self.__archtype = archtype_enum.value

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
            "generated_code_dir":       self.generated_code_dir,
            "archtype":                 self.archtype,
        }

    @staticmethod
    def camel_to_snake(name: str) -> str:
        s1 = re.sub(r'(.)([A-Z][a-z]+)', r'\1_\2', name)
        return re.sub(r'([a-z0-9])([A-Z])', r'\1_\2', s1).lower()
