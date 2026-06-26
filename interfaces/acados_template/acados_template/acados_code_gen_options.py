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
import inspect
import warnings

import casadi as ca
import numpy as np
from typing import Dict

from .utils import get_shared_lib_ext, get_acados_path, get_os_str
from sysconfig import get_paths

class AcadosCodeGenOptions:
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
        self.__acados_version = None
        self.__casadi_code_gen_options = {"mex": False, "casadi_int": "int", "casadi_real": "double", "force_canonical": False}

        env = os.environ
        self.__ext_fun_compile_flags = '-O2' if 'ACADOS_EXT_FUN_COMPILE_FLAGS' not in env else env['ACADOS_EXT_FUN_COMPILE_FLAGS']
        self.__ext_fun_expand_constr = False
        self.__ext_fun_expand_cost = False
        self.__ext_fun_expand_precompute = False
        self.__ext_fun_expand_dyn = False
        self.__model_external_shared_lib_dir = None
        self.__model_external_shared_lib_name = None

        self.__with_solution_sens_wrt_params = False
        self.__with_value_sens_wrt_params = False
        self.__generate_hess = False

        self.__sens_forw_p = False

    # read-only properties
    @property
    def acados_include_path(self) -> str:
        return self.__acados_include_path

    @acados_include_path.setter
    def acados_include_path(self, path: str) -> None:
        if not isinstance(path, str):
            raise TypeError("acados_include_path must be a string")
        self.__acados_include_path = os.path.abspath(path).replace(os.sep, '/')

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

    @property
    def acados_version(self) -> str:
        """The acados version is detected automatically. It is required for verifying if code reuse is possible."""
        return self.__acados_version

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

    @property
    def casadi_code_gen_options(self):
        """
        Options to be passed to CasADi code generation.
        Default: {"mex": False, "casadi_int": "int", "casadi_real": "double"}.
        """
        return self.__casadi_code_gen_options

    @casadi_code_gen_options.setter
    def casadi_code_gen_options(self, opts):
        if not isinstance(opts, dict):
            raise TypeError("casadi_code_gen_options must be a dictionary")
        self.__casadi_code_gen_options = opts

    @property
    def ext_fun_compile_flags(self):
        """
        String with compiler flags for external function compilation.
        Default: '-O2' if environment variable ACADOS_EXT_FUN_COMPILE_FLAGS is not set, else ACADOS_EXT_FUN_COMPILE_FLAGS is used as default.
        """
        return self.__ext_fun_compile_flags

    @ext_fun_compile_flags.setter
    def ext_fun_compile_flags(self, ext_fun_compile_flags):
        if isinstance(ext_fun_compile_flags, str):
            self.__ext_fun_compile_flags = ext_fun_compile_flags
        else:
            raise TypeError('Invalid ext_fun_compile_flags value, expected a string.\n')

    @property
    def ext_fun_expand_constr(self):
        """
        Flag indicating whether CasADi.MX should be expanded to CasADi.SX before code generation for constraint functions.
        Default: False
        """
        return self.__ext_fun_expand_constr

    @ext_fun_expand_constr.setter
    def ext_fun_expand_constr(self, ext_fun_expand_constr):
        if not isinstance(ext_fun_expand_constr, bool):
            raise TypeError('Invalid ext_fun_expand_constr value, expected bool.\n')
        self.__ext_fun_expand_constr = ext_fun_expand_constr

    @property
    def ext_fun_expand_cost(self):
        """
        Flag indicating whether CasADi.MX should be expanded to CasADi.SX before code generation for cost functions.
        Default: False
        """
        return self.__ext_fun_expand_cost

    @ext_fun_expand_cost.setter
    def ext_fun_expand_cost(self, ext_fun_expand_cost):
        if not isinstance(ext_fun_expand_cost, bool):
            raise TypeError('Invalid ext_fun_expand_cost value, expected bool.\n')
        self.__ext_fun_expand_cost = ext_fun_expand_cost

    @property
    def ext_fun_expand_dyn(self):
        """
        Flag indicating whether CasADi.MX should be expanded to CasADi.SX before code generation for dynamics functions.
        Default: False
        """
        return self.__ext_fun_expand_dyn

    @ext_fun_expand_dyn.setter
    def ext_fun_expand_dyn(self, ext_fun_expand_dyn):
        if not isinstance(ext_fun_expand_dyn, bool):
            raise TypeError('Invalid ext_fun_expand_dyn value, expected bool.\n')
        self.__ext_fun_expand_dyn = ext_fun_expand_dyn

    @property
    def ext_fun_expand_precompute(self):
        """
        Flag indicating whether CasADi.MX should be expanded to CasADi.SX before code generation for the precompute function.
        Default: False
        """
        return self.__ext_fun_expand_precompute

    @ext_fun_expand_precompute.setter
    def ext_fun_expand_precompute(self, ext_fun_expand_precompute):
        if not isinstance(ext_fun_expand_precompute, bool):
            raise TypeError('Invalid ext_fun_expand_precompute value, expected bool.\n')
        self.__ext_fun_expand_precompute = ext_fun_expand_precompute

    @property
    def model_external_shared_lib_dir(self):
        """Path to the .so lib"""
        return self.__model_external_shared_lib_dir

    @model_external_shared_lib_dir.setter
    def model_external_shared_lib_dir(self, model_external_shared_lib_dir):
        if isinstance(model_external_shared_lib_dir, str) :
            self.__model_external_shared_lib_dir = model_external_shared_lib_dir
        else:
            raise TypeError('Invalid model_external_shared_lib_dir value. Str expected.' \
            + '.\n\nYou have: ' + type(model_external_shared_lib_dir) + '.\n\n')

    @property
    def model_external_shared_lib_name(self):
        """Name of the .so lib"""
        return self.__model_external_shared_lib_name

    @model_external_shared_lib_name.setter
    def model_external_shared_lib_name(self, model_external_shared_lib_name):
        if isinstance(model_external_shared_lib_name, str) :
            if model_external_shared_lib_name[-3:] == '.so' :
                raise ValueError('Invalid model_external_shared_lib_name value. Remove the .so extension.' \
            + '.\n\nYou have: ' + type(model_external_shared_lib_name) + '.\n\n')
            else :
                self.__model_external_shared_lib_name = model_external_shared_lib_name
        else:
            raise TypeError('Invalid model_external_shared_lib_name value. Str expected.'
            + '.\n\nYou have: ' + type(model_external_shared_lib_name) + '.\n\n')

    @property
    def with_solution_sens_wrt_params(self):
        """
        Flag indicating whether solution sensitivities wrt. parameters can be computed.
        """
        return self.__with_solution_sens_wrt_params

    @with_solution_sens_wrt_params.setter
    def with_solution_sens_wrt_params(self, with_solution_sens_wrt_params):
        if isinstance(with_solution_sens_wrt_params, bool):
            self.__with_solution_sens_wrt_params = with_solution_sens_wrt_params
        else:
            raise TypeError('Invalid with_solution_sens_wrt_params value. Expected bool.')

    @property
    def with_value_sens_wrt_params(self):
        """
        Flag indicating whether value function sensitivities wrt. parameters can be computed.
        """
        return self.__with_value_sens_wrt_params

    @with_value_sens_wrt_params.setter
    def with_value_sens_wrt_params(self, with_value_sens_wrt_params):
        if isinstance(with_value_sens_wrt_params, bool):
            self.__with_value_sens_wrt_params = with_value_sens_wrt_params
        else:
            raise TypeError('Invalid with_value_sens_wrt_params value. Expected bool.')

    @property
    def generate_hess(self):
        """
        Flag indicating whether Hessian code should be generated.
        """
        return self.__generate_hess

    @generate_hess.setter
    def generate_hess(self, generate_hess):
        if isinstance(generate_hess, bool):
            self.__generate_hess = generate_hess
        else:
            raise TypeError('Invalid generate_hess value. Expected bool.')

    @property
    def sens_forw_p(self):
        """Boolean determining if integrator should support forward sensitivities with respect to parameters. Default: False"""
        return self.__sens_forw_p

    @sens_forw_p.setter
    def sens_forw_p(self, sens_forw_p):
        if isinstance(sens_forw_p, bool):
            self.__sens_forw_p = sens_forw_p
        else:
            raise TypeError('Invalid sens_forw_p value. Expected bool.')


    def make_consistent(self,) -> None:
        """
        Load link_libs.json from acados_lib_path and store ordered dict
        into acados_link_libs.
        """

        json_path = os.path.join(self.acados_lib_path, 'link_libs.json')
        with open(json_path) as f:
            self.__acados_link_libs = json.load(f)

            git_hash_file = os.path.join(self.acados_lib_path, 'git_commit_hash')
            try:
                with open(git_hash_file, 'r') as f:
                    self.__acados_version = f.read().strip()
            except Exception:
                warnings.warn("Could not read acados version from git_commit_hash file.")
                self.__acados_version = None

        self.code_export_directory = os.path.abspath(self.code_export_directory)

        # CasADi codegen options
        if self.casadi_code_gen_options.get("mex") is not False:
            warnings.warn("casadi_code_gen_options['mex'] is set to True, this is not supported by acados. Setting it to False.")
            self.casadi_code_gen_options["mex"] = False

        if self.casadi_code_gen_options.get("casadi_int") != 'int':
            warnings.warn("casadi_code_gen_options['casadi_int'] is set to a value other than 'int', this is not supported by acados. Setting it to 'int'.")
            self.casadi_code_gen_options["casadi_int"] = 'int'

        if self.casadi_code_gen_options.get("casadi_real") != 'double':
            warnings.warn("casadi_code_gen_options['casadi_real'] is set to a value other than 'double', this is not supported by acados. Setting it to 'double'.")
            self.casadi_code_gen_options["casadi_real"] = 'double'

        for k, v in list(self.casadi_code_gen_options.items()):
            try:
                ca.CodeGenerator("foo", {k: v})
            except:
                warnings.warn(f"CasADi codegen option {k} is not supported by this version of CasADi, removing it from casadi_code_gen_options.")
                del self.casadi_code_gen_options[k]

    @classmethod
    def from_dict(cls, dict):
        """
        Load all properties from a given dictionary (obtained from loading a generated json).
        Values that correspond to the empty list are ignored.
        """

        options = cls()

        # loop over all properties
        for attr, _ in inspect.getmembers(type(options), lambda v: isinstance(v, property)):

            value = dict.get(attr)

            if value is None:
                warnings.warn(f"Attribute {attr} not in dictionary.")
            else:
                try:
                    # check whether value is not the empty list
                    if not (isinstance(value, list) and not value):
                        setattr(options, attr, value)
                except Exception as e:
                    ValueError("Failed to load attribute {attr} from dictionary:\n" + repr(e))

        return options
