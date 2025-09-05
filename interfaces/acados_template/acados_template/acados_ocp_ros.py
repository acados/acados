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

import os
import re

from enum import Enum
    
class ControlLoopExec(str, Enum):
    TOPIC = "topic"
    TIMER = "timer"
    
# --- Ros Options ---
class AcadosOcpRos:
    def __init__(self):
        self.__package_name: str = "acados_ros"
        self.__node_name: str = "acados_solver_node"
        self.__namespace: str = ""
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

    @control_loop_executor.setter
    def control_loop_executor(self, control_loop_executor: ControlLoopExec | str):
        if isinstance(control_loop_executor, ControlLoopExec):
            self.__control_loop_executor = control_loop_executor.value
        elif isinstance(control_loop_executor, str) and control_loop_executor in [e.value for e in ControlLoopExec]:
            self.__control_loop_executor = control_loop_executor
        else:
            raise TypeError('Invalid control_loop_executor value, expected ControlLoopExec or str.\n')
    
    def to_dict(self) -> dict:
        self.package_name = self.__camel_to_snake(self.package_name)
        self.node_name = self.__camel_to_snake(self.node_name)
        self.namespace = self.__camel_to_snake(self.namespace)
        return {
            "package_name": self.package_name,
            "node_name": self.node_name,
            "namespace": self.namespace,
            "control_loop_executor": self.control_loop_executor
        }
        
    @staticmethod
    def __camel_to_snake(name: str) -> str:
        s1 = re.sub(r'(.)([A-Z][a-z]+)', r'\1_\2', name)
        return re.sub(r'([a-z0-9])([A-Z])', r'\1_\2', s1).lower()


if __name__ == "__main__":
    ros_opt = AcadosOcpRos()
    
    ros_opt.node_name = "my_node"
    ros_opt.control_loop_executor = ControlLoopExec.TIMER
    
    print(ros_opt.to_dict())