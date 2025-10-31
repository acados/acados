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

from .utils import ArchType, AcadosRosBaseOptions

# --- Ros Options ---
class AcadosOcpRosOptions(AcadosRosBaseOptions):
    def __init__(self):
        super().__init__()
        self.package_name: str = "acados_ocp"
        self.node_name: str = ""
        self.namespace: str = ""
        self.archtype: str = ArchType.NODE.value

        # Topics
        self.__control_topic: str = "ocp_control"
        self.__state_topic: str = "ocp_state"
        self.__parameters_topic: str = "ocp_params"
        self.__reference_topic: str = "ocp_reference"

        # Other options
        self.__threads: int = 1
        self.__publish_control_sequence: bool = False

    @property
    def control_topic(self) -> str:
        return self.__control_topic

    @property
    def state_topic(self) -> str:
        return self.__state_topic

    @property
    def parameters_topic(self) -> str:
        return self.__parameters_topic

    @property
    def reference_topic(self) -> str:
        return self.__reference_topic

    @property
    def threads(self) -> str:
        return self.__threads

    @property
    def publish_control_sequence(self) -> bool:
        return self.__publish_control_sequence

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

    @parameters_topic.setter
    def parameters_topic(self, value: str):
        if not isinstance(value, str):
            raise TypeError('Invalid parameters_topic value, expected str.\n')
        self.__parameters_topic = value

    @reference_topic.setter
    def reference_topic(self, value: str):
        if not isinstance(value, str):
            raise TypeError('Invalid reference_topic value, expected str.\n')
        self.__reference_topic = value

    @threads.setter
    def threads(self, value: int):
        if not isinstance(value, int):
            raise TypeError('Invalid threads value, expected int.\n')
        self.__threads = value

    @publish_control_sequence.setter
    def publish_control_sequence(self, value: bool):
        if not isinstance(value, bool):
            raise TypeError('Invalid publish_control_sequence value, expected bool.\n')
        self.__publish_control_sequence = value

    def to_dict(self) -> dict:
        return super().to_dict() | {
            "control_topic": self.control_topic,
            "state_topic": self.state_topic,
            "parameters_topic": self.parameters_topic,
            "reference_topic": self.reference_topic,
            "threads": self.threads,
            "publish_control_sequence": self.publish_control_sequence,
        }