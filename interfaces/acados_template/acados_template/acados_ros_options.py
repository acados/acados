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

from enum import Enum
from typing import Any

    
class ParamType(str, Enum):
    DOUBLE = "double"
    FLOAT = "float"
    INT = "int"
    BOOL = "bool"
    STRING = "str"
    DOUBLE_VECTOR = "std::vector<double>"
    INT_VECTOR = "std::vector<int>"
    
class ControlLoopExec(str, Enum):
    TOPIC = "topic"
    TIMER = "timer"

    

# --- Ros Topic ---
class AcadosRosTopic:
    def __init__(self):
        self.__name: str            = ""
        self.__msg_type: str        = ""
        self.__qos_profile: str     = "default"
        self.__description: str     = ""
        
    @property
    def name(self):
        """Name of the topic.
        Default: ''.
        """
        return self.__name
    
    @property
    def msg_type(self):
        """Message type of the topic.
        Default: ''.
        """
        return self.__msg_type
    
    @property
    def qos_profile(self):
        """QoS profile of the topic.
        Default: 'default'.
        """
        return self.__qos_profile
    
    @property
    def description(self):
        """Description of the topic.
        Default: 'default'.
        """
        return self.__description
    
    @name.setter
    def name(self, name: str):
        if not isinstance(name, str):
            raise TypeError('Invalid name value, expected str.\n')
        self.__name = name
        
    @msg_type.setter
    def msg_type(self, msg_type: str):
        if not isinstance(msg_type, str):
            raise TypeError('Invalid msg_type value, expected str.\n')
        self.__msg_type = msg_type

    @qos_profile.setter
    def qos_profile(self, qos_profile: str):
        if not isinstance(qos_profile, str):
            raise TypeError('Invalid qos_profile value, expected str.\n')
        self.__qos_profile = qos_profile
        
    @description.setter
    def description(self, description: str):
        if not isinstance(description, str):
            raise TypeError('Invalid description value, expected str.\n')
        self.__description = description

    def __repr__(self) -> str:
        return (f"{self.__class__.__name__}(name={self.name!r}, msg_type={self.msg_type!r}, "
                f"qos_profile={self.qos_profile!r}, description={self.description!r})")

    def __str__(self) -> str:
        return f"Topic '{self.name}': type={self.msg_type}"
        
    def to_dict(self) -> dict:
        if not self.name or not self.msg_type or not self.qos_profile:
            raise ValueError("Topic name, msg_type, and qos_profile must be set.")

        return {
            "name": self.name,
            "msg_type": self.msg_type,
            "qos_profile": self.qos_profile,
            "description": self.description
        }
        

# --- Ros Publisher ---
class AcadosRosPublisher(AcadosRosTopic):
    def __init__(self):
        super().__init__()
        self.__queue_size: str = 10
        
    @property
    def queue_size(self):
        """Queue size of the publisher.
        Default: 10.
        """
        return self.__queue_size
    
    @queue_size.setter
    def queue_size(self, queue_size: int):
        if not isinstance(queue_size, int):
            raise TypeError('Invalid queue_size value, expected int.\n')
        self.__queue_size = queue_size
        
    def __repr__(self) -> str:
        return (super().__repr__()[:-1] +
                f", queue_size={self.queue_size!r})")
    
    def __str__(self) -> str:
        return f"Publisher '{self.name}': type={self.msg_type}"
        
    def to_dict(self) -> dict:
        return super().to_dict() | {
            "queue_size": self.queue_size
        }
    
    
# --- Ros Subscriber ---
class AcadosRosSubscriber(AcadosRosTopic):
    def __init__(self):
        super().__init__()

    def __str__(self) -> str:
        return f"Subscriber '{self.name}': type={self.msg_type}"
        
    def to_dict(self) -> dict:
        return super().to_dict()
    
        
# --- Ros Services ---
class AcadosRosService:
    def __init__(self):
        self.__name: str = "srv_name"
        self.__srv_type: str = "srv_type"
        self.__description: str = ""
        self.__template_hint: str = "template_hint"

    @property
    def name(self):
        """Name of the service.
        Default: 'srv_name'.
        """
        return self.__name

    @property
    def srv_type(self):
        """Service type.
        Default: 'srv_type'.
        """
        return self.__srv_type

    @property
    def purpose(self):
        """Purpose/description of the service.
        Default: ''.
        """
        return self.__description

    @property
    def template_hint(self):
        """Template hint for the service generation.
        Default: 'template_hint'.
        """
        return self.__template_hint

    @name.setter
    def name(self, name: str):
        if not isinstance(name, str):
            raise TypeError('Invalid name value, expected str.\n')
        self.__name = name

    @srv_type.setter
    def srv_type(self, srv_type: str):
        if not isinstance(srv_type, str):
            raise TypeError('Invalid srv_type value, expected str.\n')
        self.__srv_type = srv_type

    @purpose.setter
    def purpose(self, purpose: str):
        if not isinstance(purpose, str):
            raise TypeError('Invalid purpose value, expected str.\n')
        self.__description = purpose

    @template_hint.setter
    def template_hint(self, template_hint: str):
        if not isinstance(template_hint, str):
            raise TypeError('Invalid template_hint value, expected str.\n')
        self.__template_hint = template_hint

    def __repr__(self) -> str:
        return (f"{self.__class__.__name__}(name={self.name!r}, srv_type={self.srv_type!r}, "
                f"purpose={self.purpose!r}, template_hint={self.template_hint!r})")

    def __str__(self) -> str:
        return (f"Service '{self.name}': type={self.srv_type}, purpose={self.purpose}")
    
    def to_dict(self) -> dict:
        return {
            "name": self.name,
            "srv_type": self.srv_type,
            "purpose": self.purpose,
            "template_hint": self.template_hint
        }


# --- Ros Parameters ---
class AcadosRosParameter:
    def __init__(self):
        self.__name: str = ""
        self.__p_type: str = ParamType.STRING.value
        self.__default: str = ""
        self.__description: str = ""

    @property
    def name(self):
        """Name of the parameter.
        Default: ''.
        """
        return self.__name

    @property
    def p_type(self):
        """Parameter type (as str or ParamType enum).
        Default: ''.
        """
        return self.__p_type

    @property
    def default(self):
        """Default value for the parameter.
        """
        return self.__default

    @property
    def description(self):
        """Description of the parameter.
        Default: ''.
        """
        return self.__description

    @name.setter
    def name(self, name):
        if not isinstance(name, str):
            raise TypeError('Invalid name value, expected str.\n')
        self.__name = name

    @p_type.setter
    def p_type(self, p_type: ParamType | str):
        if isinstance(p_type, ParamType):
            self.__p_type = p_type.value
        elif isinstance(p_type, str) and p_type in [e.value for e in ParamType]:
            self.__p_type = p_type
        else:
            raise TypeError('Invalid p_type value, expected ParamType enum or str.\n')

    @default.setter
    def default(self, default: Any):
        # allow any type for default (could be numeric, str, list, etc.)
        self.__default = default

    @description.setter
    def description(self, description: str):
        if not isinstance(description, str):
            raise TypeError('Invalid description value, expected str.\n')
        self.__description = description

    def __repr__(self) -> str:
        return (f"{self.__class__.__name__}(name={self.name!r}, p_type={self.p_type!r}, "
                f"default={self.default!r}, description={self.description!r})")

    def __str__(self) -> str:
        return (f"Parameter '{self.name}': type={self.p_type}, default={self.default}")
    
    def to_dict(self) -> dict:
        if not self.name or not self.p_type:
            raise ValueError("Parameter name and p_type must be set.")
        return {
            "name": self.name,
            "type": self.p_type,
            "default": self.default,
            "description": self.description
        }
        
# --- Ros Package ---
class AcadosRosPackage:
    def __init__(self):
        self.__package_name: str        = "acados_solver"
        self.__version: str             = "0.1.0"
        self.__description: str         = "ACADOS ROS 2 Interface"
        self.__author_name: str         = "Your Name"
        self.__author_email: str        = "your.name@email.com"
        self.__license: str             = "MY LICENSE"
        
    @property
    def package_name(self):
        """Name of the package.
        Default: 'acados_solver'.
        """
        return self.__package_name

    @property
    def version(self):
        """Version of the package.
        Default: '0.1.0'.
        """
        return self.__version

    @property
    def description(self):
        """Short description of the package.
        Default: 'ACADOS ROS 2 Interface'.
        """
        return self.__description

    @property
    def author_name(self):
        """Author name of the package.
        Default: 'Your Name'.
        """
        return self.__author_name

    @property
    def author_email(self):
        """Author email of the package.
        Default: 'your.name@email.com'.
        """
        return self.__author_email

    @property
    def license(self):
        """License of the package.
        Default: 'MY LICENSE'.
        """
        return self.__license
    
    @package_name.setter
    def package_name(self, package_name: str):
        if not isinstance(package_name, str):
            raise TypeError('Invalid package_name value, expected str.\n')
        self.__package_name = package_name

    @version.setter
    def version(self, version: str):
        if not isinstance(version, str):
            raise TypeError('Invalid version value, expected str.\n')
        self.__version = version

    @description.setter
    def description(self, description: str):
        if not isinstance(description, str):
            raise TypeError('Invalid description value, expected str.\n')
        self.__description = description

    @author_name.setter
    def author_name(self, author_name: str):
        if not isinstance(author_name, str):
            raise TypeError('Invalid author_name value, expected str.\n')
        self.__author_name = author_name

    @author_email.setter
    def author_email(self, author_email: str):
        if not isinstance(author_email, str):
            raise TypeError(f'Invalid author_email value, expected str.\n')
        self.__author_email = author_email

    @license.setter
    def license(self, license: str):
        if not isinstance(license, str):
            raise TypeError('Invalid license value, expected str.\n')
        self.__license = license

    def __repr__(self) -> str:
        return (f"{self.__class__.__name__}(name={self.package_name!r}, version={self.version!r}, "
                f"author={self.author_name!r}, email={self.author_email!r}, "
                f"license={self.license!r})")

    def __str__(self) -> str:
        return (f"Package {self.package_name} v{self.version} by {self.author_name} <{self.author_email}>")
    
    def to_dict(self) -> dict:
        return {
            "package_name": self.package_name,
            "version": self.version,
            "description": self.description,
            "author_name": self.author_name,
            "author_email": self.author_email,
            "license": self.license
        }

# --- Ros Options ---
class AcadosRosOptions:
    def __init__(self):
        self.__package: AcadosRosPackage = AcadosRosPackage()
        self.__node_name: str = "acados_solver_node"
        self.__namespace: str = ""
        self.__publishers: list[AcadosRosPublisher] = []
        self.__subscribers: list[AcadosRosSubscriber] = []
        self.__services: list[AcadosRosService] = []
        self.__parameters: list[AcadosRosParameter] = []
        self.__control_loop_executor: str = ControlLoopExec.TIMER.value
        self.__extra_include_dirs: set[str] = set()
        self.__extra_link_libs: set[str] = set()

    @property
    def package(self) -> AcadosRosPackage:
        return self.__package

    @property
    def node_name(self) -> str:
        return self.__node_name

    @property
    def namespace(self) -> str:
        return self.__namespace

    @property
    def publishers(self) -> list[AcadosRosPublisher]:
        return self.__publishers
    
    @property
    def subscribers(self) -> list[AcadosRosSubscriber]:
        return self.__subscribers

    @property
    def services(self) -> list[AcadosRosService]:
        return self.__services

    @property
    def parameters(self) -> list[AcadosRosParameter]:
        return self.__parameters

    @property
    def control_loop_executor(self) -> str:
        return self.__control_loop_executor

    @property
    def extra_include_dirs(self) -> set[str]:
        return self.__extra_include_dirs

    @property
    def extra_link_libs(self) -> set[str]:
        return self.__extra_link_libs
    
    @package.setter
    def package(self, package: str):
        if not isinstance(package, str):
            raise TypeError('Invalid package value, expected str.\n')
        self.__package = package
        
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

    @publishers.setter
    def publishers(self, publishers: list[AcadosRosPublisher]):
        if not all(isinstance(p, AcadosRosPublisher) for p in publishers):
            raise TypeError('Invalid publishers value, expected iterable of AcadosRosPublisher.\n')
        self.__publishers = publishers
        
    @subscribers.setter
    def subscribers(self, subscribers: list[AcadosRosSubscriber]):
        if not all(isinstance(s, AcadosRosSubscriber) for s in subscribers):
            raise TypeError('Invalid subscribers value, expected iterable of AcadosRosSubscriber.\n')
        self.__subscribers = subscribers

    @services.setter
    def services(self, services: list[AcadosRosService]):
        if not all(isinstance(s, AcadosRosService) for s in services):
            raise TypeError('Invalid services value, expected iterable of AcadosRosService.\n')
        self.__services = services

    @parameters.setter
    def parameters(self, parameters: list[AcadosRosParameter]):
        if not all(isinstance(p, AcadosRosParameter) for p in parameters):
            raise TypeError('Invalid parameters value, expected iterable of AcadosRosParameter.\n')
        self.__parameters = parameters

    @control_loop_executor.setter
    def control_loop_executor(self, control_loop_executor: ControlLoopExec | str):
        if isinstance(control_loop_executor, ControlLoopExec):
            self.__control_loop_executor = control_loop_executor.value
        elif isinstance(control_loop_executor, str) and control_loop_executor in [e.value for e in ControlLoopExec]:
            self.__control_loop_executor = control_loop_executor
        else:
            raise TypeError('Invalid control_loop_executor value, expected ControlLoopExec or str.\n')

    @extra_include_dirs.setter
    def extra_include_dirs(self, extra_include_dirs: set[str]):
        if not all(isinstance(d, str) for d in extra_include_dirs):
            raise TypeError('Invalid extra_include_dirs value, expected iterable of str.\n')
        self.__extra_include_dirs = set(extra_include_dirs)

    @extra_link_libs.setter
    def extra_link_libs(self, extra_link_libs: set[str]):
        if not all(isinstance(l, str) for l in extra_link_libs):
            raise TypeError('Invalid extra_link_libs value, expected iterable of str.\n')
        self.__extra_link_libs = set(extra_link_libs)
        
    def __repr__(self):
        pkg_name = self.package.package_name if hasattr(self.package, 'package_name') else repr(self.package)
        includes = sorted(list(self.extra_include_dirs)) if self.extra_include_dirs else []
        links = sorted(list(self.extra_link_libs)) if self.extra_link_libs else []
        return (
            f"{self.__class__.__name__}(package={pkg_name!r}, node_name={self.node_name!r}, "
            f"namespace={self.namespace!r}, publishers={len(self.publishers)}, "
            f"subscribers={len(self.subscribers)}, services={len(self.services)}, "
            f"parameters={len(self.parameters)}, control_loop_executor={self.control_loop_executor!r}, "
            f"include_dirs={includes!r}, link_libs={links!r})"
        )

    def __str__(self) -> str:
        lines = [f"{self.__class__.__name__}:",
                 f"  package: {self.package}",
                 f"  node: {self.node_name}",
                 f"  namespace: {self.namespace}",
                 f"  control_loop_executor: {self.control_loop_executor}",
                 f"  include_dirs: {sorted(self.extra_include_dirs)}",
                 f"  link_libs: {sorted(self.extra_link_libs)}",
                 f"  publishers:"]
        for p in self.publishers:
            lines.append(f"    - {p}")
        lines.append("  subscribers:")
        for s in self.subscribers:
            lines.append(f"    - {s}")
        lines.append("  services:")
        for s in self.services:
            lines.append(f"    - {s}")
        lines.append("  parameters:")
        for p in self.parameters:
            lines.append(f"    - {p}")
        return "\n".join(lines)
    
    def to_dict(self) -> dict:
        return {
            "package": self.package.to_dict() if hasattr(self.package, 'to_dict') else str(self.package),
            "node_name": self.node_name,
            "namespace": self.namespace,
            "publishers": [p.to_dict() for p in self.publishers],
            "subscribers": [s.to_dict() for s in self.subscribers],
            "services": [s.to_dict() for s in self.services],
            "parameters": [p.to_dict() for p in self.parameters],
            "control_loop_executor": self.control_loop_executor,
            "extra_include_dirs": sorted(list(self.extra_include_dirs)),
            "extra_link_libs": sorted(list(self.extra_link_libs))
        }


if __name__ == "__main__":
    ros_opt = AcadosRosOptions()
    
    ros_opt.node_name = "my_node"
    ros_opt.publishers = [AcadosRosPublisher()]
    ros_opt.subscribers = [AcadosRosSubscriber()]
    ros_opt.control_loop_executor = ControlLoopExec.TIMER
    ros_opt.extra_include_dirs = ["include", "2"]
    ros_opt.extra_link_libs = ["libacados.so"]
    
    ros_opt.publishers[0].name = "pub1"
    ros_opt.publishers[0].msg_type = "std_msgs/msg/Float64"
    ros_opt.publishers[0].queue_size = 5
    
    ros_opt.subscribers[0].name = "sub1"
    ros_opt.subscribers[0].msg_type = "std_msgs/msg/Float64"

    print(ros_opt.to_dict())