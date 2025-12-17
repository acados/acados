#
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
import json
import numpy as np
from typing import Dict

from .utils import get_shared_lib_ext, get_acados_path, get_os_str
from sysconfig import get_paths

class AcadosCodeGenOpts:
    def __init__(self) -> None:
        acados_path = get_acados_path()

        self.__acados_include_path = os.path.join(acados_path, 'include').replace(os.sep, '/')
        self.__acados_link_libs: Dict[str, str] = {}

        self.__shared_lib_ext = get_shared_lib_ext()
        self.__os = get_os_str()
        self.__cython_include_dirs = [np.get_include(), get_paths()['include']]

        self.__acados_lib_path = os.path.join(acados_path, 'lib')
        self.__json_file: str = ''
        self.__code_export_directory = 'c_generated_code'

    # read-only properties
    @property
    def acados_include_path(self) -> str:
        return self.__acados_include_path

    @property
    def acados_link_libs(self) -> Dict[str, str]:
        return self.__acados_link_libs

    @property
    def shared_lib_ext(self) -> str:
        return self.__shared_lib_ext

    @property
    def os(self) -> str:
        return self.__os

    @property
    def cython_include_dirs(self) -> list:
        return self.__cython_include_dirs

    # public properties with setters
    @property
    def acados_lib_path(self) -> str:
        return self.__acados_lib_path

    @acados_lib_path.setter
    def acados_lib_path(self, path: str) -> None:
        if not isinstance(path, str):
            raise TypeError("acados_lib_path must be a string")
        self.__acados_lib_path = os.path.abspath(path)

    @property
    def json_file(self) -> str:
        return self.__json_file

    @json_file.setter
    def json_file(self, filename: str) -> None:
        if not isinstance(filename, str):
            raise TypeError("json_file must be a string")
        self.__json_file = filename

    @property
    def code_export_directory(self) -> str:
        return self.__code_export_directory

    @code_export_directory.setter
    def code_export_directory(self, directory: str) -> None:
        if not isinstance(directory, str):
            raise TypeError("code_export_directory must be a string")
        # store as absolute path for consistency
        self.__code_export_directory = os.path.abspath(directory)

    def make_consistent(self) -> None:
        """
        Load link_libs.json from acados_lib_path and store ordered dict
        into acados_link_libs.
        """
        json_path = os.path.join(self.acados_lib_path, 'link_libs.json')
        with open(json_path) as f:
            self.__acados_link_libs = json.load(f)
