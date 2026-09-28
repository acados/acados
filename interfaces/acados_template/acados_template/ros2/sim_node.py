# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# Licensed under the 2-Clause BSD License.

from .utils import ArchType, AcadosRosBaseOptions


# --- Ros Options ---
class AcadosSimRosOptions(AcadosRosBaseOptions):
    def __init__(self):
        super().__init__()
        self.package_name: str = "acados_sim"
        self.node_name: str = ""
        self.namespace: str = ""
        self.archtype: str = ArchType.NODE.value

        self.__control_topic = "sim_control"
        self.__state_topic = "sim_state"

    @property
    def control_topic(self) -> str:
        return self.__control_topic

    @property
    def state_topic(self) -> str:
        return self.__state_topic

    @control_topic.setter
    def control_topic(self, value: str):
        if not isinstance(value, str):
            raise TypeError('Invalid control_topic value, expected str.\n')
        self.__control_topic = value

    @state_topic.setter
    def state_topic(self, value: str):
        if not isinstance(value, str):
            raise TypeError('Invalid state_topic value, expected str.\n')
        self.__state_topic = value

    def to_dict(self) -> dict:
        return super().to_dict() | {
            "control_topic": self.control_topic,
            "state_topic": self.state_topic,
        }