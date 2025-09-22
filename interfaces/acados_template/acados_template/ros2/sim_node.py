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