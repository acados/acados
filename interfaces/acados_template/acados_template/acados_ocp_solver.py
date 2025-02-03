# -*- coding: future_fstrings -*-
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

import importlib
import json
import os
import shutil
import sys
import time

from ctypes import (POINTER, byref, c_char_p, c_double, c_int, c_bool,
                    c_void_p, cast)
if os.name == 'nt':
    from ctypes import wintypes
    from ctypes import WinDLL as DllLoader
else:
    from ctypes import CDLL as DllLoader
from datetime import datetime
from typing import Union, Optional, List, Tuple, Sequence, Dict

import numpy as np
import scipy.linalg
from .builders import CMakeBuilder
from .acados_ocp import AcadosOcp
from .acados_multiphase_ocp import AcadosMultiphaseOcp
from .gnsf.detect_gnsf_structure import detect_gnsf_structure
from .utils import (get_shared_lib_ext, get_shared_lib_prefix, get_shared_lib_dir, get_shared_lib,
                    make_object_json_dumpable, set_up_imported_gnsf_model, verbose_system_call,
                    acados_lib_is_compiled_with_openmp, is_empty, set_directory)
from .acados_ocp_iterate import AcadosOcpIterate, AcadosOcpIterates, AcadosOcpFlattenedIterate


class AcadosOcpSolver:
    """
    Class to interact with the acados ocp solver C object.

        :param acados_ocp: type :py:class:`~acados_template.acados_ocp.AcadosOcp` or :py:class:`~acados_template.acados_multiphase_ocp.AcadosMultiphaseOcp` - description of the OCP for acados
        :param json_file: name for the json file used to render the templated code - default: acados_ocp_nlp.json
    """
    if os.name == 'nt':
        dlclose = DllLoader('kernel32', use_last_error=True).FreeLibrary
        dlclose.argtypes = [wintypes.HMODULE]
        winmode = 8 # why 8? what does that mean?
    else:
        dlclose = DllLoader(None).dlclose
        dlclose.argtypes = [c_void_p]
        winmode = None

    @property
    def acados_lib_uses_omp(self,):
        """`acados_lib_uses_omp` - flag indicating whether the acados library has been compiled with openMP."""
        return self.__acados_lib_uses_omp

    @property
    def acados_lib(self,):
        """`acados_lib` - acados shared library"""
        return self.__acados_lib

    @property
    def shared_lib(self,):
        """`shared_lib` - solver shared library"""
        return self.__shared_lib

    # TODO move this to AcadosOcp
    @classmethod
    def generate(cls, acados_ocp: Union[AcadosOcp, AcadosMultiphaseOcp], json_file: str, simulink_opts=None, cmake_builder: CMakeBuilder = None):
        """
        Generates the code for an acados OCP solver, given the description in acados_ocp.
            :param acados_ocp: type Union[AcadosOcp, AcadosMultiphaseOcp] - description of the OCP for acados
            :param json_file: name for the json file used to render the templated code - default: `acados_ocp_nlp.json`
            :param simulink_opts: Options to configure Simulink S-function blocks, mainly to activate possible inputs and
                   outputs; default: `None`
            :param cmake_builder: type :py:class:`~acados_template.builders.CMakeBuilder` generate a `CMakeLists.txt` and use
                   the `CMake` pipeline instead of a `Makefile` (`CMake` seems to be the better option in conjunction with
                   `MS Visual Studio`); default: `None`
        """
        acados_ocp.code_export_directory = os.path.abspath(acados_ocp.code_export_directory)
        acados_ocp.simulink_opts = simulink_opts

        # add kwargs to acados_ocp
        acados_ocp.json_file = json_file

        # make consistent
        acados_ocp.make_consistent()

        # module dependent post processing
        if acados_ocp.solver_options.integrator_type == 'GNSF':
            if 'gnsf_model' in acados_ocp.__dict__:
                set_up_imported_gnsf_model(acados_ocp)
            else:
                detect_gnsf_structure(acados_ocp)

        if acados_ocp.solver_options.qp_solver == 'PARTIAL_CONDENSING_QPDUNES':
            acados_ocp.remove_x0_elimination()

        if acados_ocp.solver_options.qp_solver in ['FULL_CONDENSING_QPOASES', 'PARTIAL_CONDENSING_QPDUNES', 'PARTIAL_CONDENSING_OSQP']:
            print(f"NOTE: The selected QP solver {acados_ocp.solver_options.qp_solver} does not support one-sided constraints yet.")

        # generate code (external functions and templated code)
        acados_ocp.generate_external_functions()
        acados_ocp.dump_to_json()
        acados_ocp.render_templates(cmake_builder=cmake_builder)

        # copy custom update function
        if acados_ocp.solver_options.custom_update_filename != "" and acados_ocp.solver_options.custom_update_copy:
            target_location = os.path.join(acados_ocp.code_export_directory, acados_ocp.solver_options.custom_update_filename)
            shutil.copyfile(acados_ocp.solver_options.custom_update_filename, target_location)
        if acados_ocp.solver_options.custom_update_header_filename != "" and acados_ocp.solver_options.custom_update_copy:
            target_location = os.path.join(acados_ocp.code_export_directory, acados_ocp.solver_options.custom_update_header_filename)
            shutil.copyfile(acados_ocp.solver_options.custom_update_header_filename, target_location)


    @classmethod
    def build(cls, code_export_dir, with_cython=False, cmake_builder: CMakeBuilder = None, verbose: bool = True):
        """
        Builds the code for an acados OCP solver, that has been generated in code_export_dir
            :param code_export_dir: directory in which acados OCP solver has been generated, see generate()
            :param with_cython: option indicating if the cython interface is build, default: False.
            :param cmake_builder: type :py:class:`~acados_template.builders.CMakeBuilder` generate a `CMakeLists.txt` and use
                   the `CMake` pipeline instead of a `Makefile` (`CMake` seems to be the better option in conjunction with
                   `MS Visual Studio`); default: `None`
            :param verbose: indicating if build command is printed
        """
        code_export_dir = os.path.abspath(code_export_dir)

        with set_directory(code_export_dir):
            if os.name == 'nt':
                make_cmd = 'mingw32-make'
            else:
                make_cmd = 'make'

            if with_cython:
                verbose_system_call([make_cmd, 'clean_all'], verbose)
                verbose_system_call([make_cmd, 'ocp_cython'], verbose)
            else:
                if cmake_builder is not None:
                    cmake_builder.exec(code_export_dir, verbose)
                else:
                    verbose_system_call([make_cmd, 'clean_ocp_shared_lib'], verbose)
                    verbose_system_call([make_cmd, 'ocp_shared_lib'], verbose)


    @classmethod
    def create_cython_solver(cls, json_file):
        """
        Returns an `AcadosOcpSolverCython` object.

        This is an alternative Cython based Python wrapper to the acados OCP solver in C.
        This offers faster interaction with the solver, because getter and setter calls, which include shape checking are done in compiled C code.

        The default wrapper `AcadosOcpSolver` is based on ctypes.
        """
        with open(json_file, 'r') as f:
            acados_ocp_json = json.load(f)
        code_export_directory = acados_ocp_json['code_export_directory']

        importlib.invalidate_caches()
        sys.path.append(os.path.dirname(code_export_directory))
        acados_ocp_solver_pyx = importlib.import_module(f'{os.path.split(code_export_directory)[1]}.acados_ocp_solver_pyx')

        AcadosOcpSolverCython = getattr(acados_ocp_solver_pyx, 'AcadosOcpSolverCython')
        return AcadosOcpSolverCython(acados_ocp_json['model']['name'],
                    acados_ocp_json['solver_options']['nlp_solver_type'],
                    acados_ocp_json['dims']['N'])

    @property
    def save_p_global(self) -> bool:
        return self.__save_p_global

    def __init__(self, acados_ocp: Union[AcadosOcp, AcadosMultiphaseOcp], json_file=None, simulink_opts=None, build=True, generate=True, cmake_builder: CMakeBuilder = None, verbose=True, save_p_global=False):

        self.solver_created = False
        self.__save_p_global = save_p_global
        if save_p_global:
            self.__p_global_values = acados_ocp.p_global_values
        else:
            self.__p_global_values = np.array([])

        if not (isinstance(acados_ocp, AcadosOcp) or isinstance(acados_ocp, AcadosMultiphaseOcp)):
            raise Exception('acados_ocp should be of type AcadosOcp or AcadosMultiphaseOcp.')

        if json_file is not None:
            acados_ocp.json_file = json_file

        if generate:
            self.generate(acados_ocp, json_file=acados_ocp.json_file, simulink_opts=simulink_opts, cmake_builder=cmake_builder)
        else:
            acados_ocp.make_consistent()

        # load json, store options in object
        with open(acados_ocp.json_file, 'r') as f:
            acados_ocp_json = json.load(f)
        if isinstance(acados_ocp, AcadosOcp):
            self.N = acados_ocp_json['dims']['N']
        elif isinstance(acados_ocp, AcadosMultiphaseOcp):
            self.N = acados_ocp_json['N_horizon']
        self.__solver_options = acados_ocp_json['solver_options']
        self.name = acados_ocp_json['name']

        acados_lib_path = acados_ocp_json['acados_lib_path']
        code_export_directory = acados_ocp_json['code_export_directory']

        if build:
            self.build(code_export_directory, with_cython=False, cmake_builder=cmake_builder, verbose=verbose)

        # prepare library loading
        lib_ext = get_shared_lib_ext()
        lib_prefix = get_shared_lib_prefix()
        lib_dir = get_shared_lib_dir()

        # Load acados library to avoid unloading the library.
        # This is necessary if acados was compiled with OpenMP, since the OpenMP threads can't be destroyed.
        # Unloading a library which uses OpenMP results in a segfault (on any platform?).
        # see [https://stackoverflow.com/questions/34439956/vc-crash-when-freeing-a-dll-built-with-openmp]
        # or [https://python.hotexamples.com/examples/_ctypes/-/dlclose/python-dlclose-function-examples.html]
        libacados_name = f'{lib_prefix}acados{lib_ext}'
        libacados_filepath = os.path.join(acados_lib_path, '..', lib_dir, libacados_name)
        self.__acados_lib = get_shared_lib(libacados_filepath, self.winmode)

        # find out if acados was compiled with OpenMP
        self.__acados_lib_uses_omp = acados_lib_is_compiled_with_openmp(self.__acados_lib, verbose)

        libacados_ocp_solver_name = f'{lib_prefix}acados_ocp_solver_{self.name}{lib_ext}'
        self.shared_lib_name = os.path.join(code_export_directory, libacados_ocp_solver_name)

        # get shared_lib
        self.__shared_lib = get_shared_lib(self.shared_lib_name, self.winmode)

        # create capsule
        getattr(self.__shared_lib, f"{self.name}_acados_create_capsule").restype = c_void_p
        self.capsule = getattr(self.__shared_lib, f"{self.name}_acados_create_capsule")()

        # create solver
        getattr(self.__shared_lib, f"{self.name}_acados_create").argtypes = [c_void_p]
        getattr(self.__shared_lib, f"{self.name}_acados_create").restype = c_int
        assert getattr(self.__shared_lib, f"{self.name}_acados_create")(self.capsule)==0
        self.solver_created = True

        self.acados_ocp = acados_ocp

        # get pointers solver
        self.__get_pointers_solver()

        self.status = 0
        self.time_solution_sens_solve = 0.0
        self.time_solution_sens_lin = 0.0

        # gettable fields
        self.__qp_dynamics_fields = ['A', 'B', 'b']
        self.__qp_cost_fields = ['Q', 'R', 'S', 'q', 'r', 'zl', 'zu', 'Zl', 'Zu']
        self.__qp_constraint_fields = ['C', 'D', 'lg', 'ug', 'lbx', 'ubx', 'lbu', 'ubu']
        self.__qp_constraint_int_fields = ['idxs', 'idxb']
        self.__qp_pc_hpipm_fields = ['P', 'K', 'Lr', 'p']
        self.__qp_pc_fields = ['pcond_Q', 'pcond_R', 'pcond_S']

        # set arg and res types
        self.__acados_lib.ocp_nlp_dims_get_from_attr.argtypes = [c_void_p, c_void_p, c_void_p, c_int, c_char_p]
        self.__acados_lib.ocp_nlp_dims_get_from_attr.restype = c_int
        self.__acados_lib.ocp_nlp_eval_params_jac.argtypes = [c_void_p, c_void_p, c_void_p]
        self.__acados_lib.ocp_nlp_eval_lagrange_grad_p.argtypes = [c_void_p, c_void_p, c_char_p, POINTER(c_double)]
        self.__acados_lib.ocp_nlp_out_get.argtypes = [c_void_p, c_void_p, c_void_p, c_int, c_char_p, c_void_p]
        self.__acados_lib.ocp_nlp_in_get.argtypes = [c_void_p, c_void_p, c_void_p, c_int, c_char_p, c_void_p]

        self.__acados_lib.ocp_nlp_eval_param_sens.argtypes = [c_void_p, c_char_p, c_int, c_int, c_void_p]
        self.__acados_lib.ocp_nlp_eval_param_sens.restype = None

        self.__acados_lib.ocp_nlp_eval_solution_sens_adj_p.argtypes = [c_void_p, c_void_p, c_void_p, c_char_p, c_int, c_void_p]
        self.__acados_lib.ocp_nlp_eval_solution_sens_adj_p.restype = None

        self.__acados_lib.ocp_nlp_solver_opts_set.argtypes = [c_void_p, c_void_p, c_char_p, c_void_p]
        self.__acados_lib.ocp_nlp_get.argtypes = [c_void_p, c_char_p, c_void_p]

        self.__acados_lib.ocp_nlp_eval_cost.argtypes = [c_void_p, c_void_p, c_void_p]
        self.__acados_lib.ocp_nlp_eval_residuals.argtypes = [c_void_p, c_void_p, c_void_p]

        self.__acados_lib.ocp_nlp_constraints_model_set.argtypes = [c_void_p, c_void_p, c_void_p, c_int, c_char_p, c_void_p]
        self.__acados_lib.ocp_nlp_constraints_model_get.argtypes = [c_void_p, c_void_p, c_void_p, c_int, c_char_p, c_void_p]
        self.__acados_lib.ocp_nlp_cost_model_set.argtypes =  [c_void_p, c_void_p, c_void_p, c_int, c_char_p, c_void_p]
        self.__acados_lib.ocp_nlp_cost_model_get.argtypes =  [c_void_p, c_void_p, c_void_p, c_int, c_char_p, c_void_p]

        self.__acados_lib.ocp_nlp_out_set.argtypes = [c_void_p, c_void_p, c_void_p, c_int, c_char_p, c_void_p]
        self.__acados_lib.ocp_nlp_set.argtypes = [c_void_p, c_int, c_char_p, c_void_p]

        self.__acados_lib.ocp_nlp_cost_dims_get_from_attr.argtypes = [c_void_p, c_void_p, c_void_p, c_int, c_char_p, POINTER(c_int)]
        self.__acados_lib.ocp_nlp_cost_dims_get_from_attr.restype = c_int

        self.__acados_lib.ocp_nlp_constraint_dims_get_from_attr.argtypes = [c_void_p, c_void_p, c_void_p, c_int, c_char_p, POINTER(c_int)]
        self.__acados_lib.ocp_nlp_constraint_dims_get_from_attr.restype = c_int

        self.__acados_lib.ocp_nlp_qp_dims_get_from_attr.argtypes = [c_void_p, c_void_p, c_void_p, c_int, c_char_p, POINTER(c_int)]
        self.__acados_lib.ocp_nlp_qp_dims_get_from_attr.restype = c_int

        self.__acados_lib.ocp_nlp_get_at_stage.argtypes = [c_void_p, c_int, c_char_p, c_void_p]

        self.__acados_lib.ocp_nlp_get_from_iterate.argtypes = [c_void_p, c_int, c_int, c_char_p, c_void_p]
        self.__acados_lib.ocp_nlp_get_from_iterate.restypes = c_void_p

        self.__acados_lib.ocp_nlp_dims_get_total_from_attr.argtypes = [c_void_p, c_void_p, c_void_p, c_char_p]
        self.__acados_lib.ocp_nlp_dims_get_total_from_attr.restype = c_int

        self.__acados_lib.ocp_nlp_get_all.argtypes = [c_void_p, c_void_p, c_void_p, c_char_p, c_void_p]
        self.__acados_lib.ocp_nlp_get_all.restype = None

        self.__acados_lib.ocp_nlp_set_all.argtypes = [c_void_p, c_void_p, c_void_p, c_char_p, c_void_p]
        self.__acados_lib.ocp_nlp_set_all.restype = None

        self.__acados_lib.ocp_nlp_out_set_values_to_zero.argtypes = [c_void_p, c_void_p, c_void_p]

        getattr(self.shared_lib, f"{self.name}_acados_solve").argtypes = [c_void_p]
        getattr(self.shared_lib, f"{self.name}_acados_solve").restype = c_int

        getattr(self.shared_lib, f"{self.name}_acados_reset").argtypes = [c_void_p, c_int]
        getattr(self.shared_lib, f"{self.name}_acados_reset").restype = c_int

        getattr(self.shared_lib, f"{self.name}_acados_custom_update").argtypes = [c_void_p, POINTER(c_double), c_int]
        getattr(self.shared_lib, f"{self.name}_acados_custom_update").restype = c_int

        getattr(self.shared_lib, f"{self.name}_acados_create_with_discretization").argtypes = [c_void_p, c_int, c_void_p]
        getattr(self.shared_lib, f"{self.name}_acados_create_with_discretization").restype = c_int

        getattr(self.shared_lib, f"{self.name}_acados_free").argtypes = [c_void_p]
        getattr(self.shared_lib, f"{self.name}_acados_free").restype = c_int

        getattr(self.shared_lib, f"{self.name}_acados_free_capsule").argtypes = [c_void_p]
        getattr(self.shared_lib, f"{self.name}_acados_free_capsule").restype = c_int

        getattr(self.shared_lib, f"{self.name}_acados_update_params_sparse").argtypes = [c_void_p, c_int, POINTER(c_int), POINTER(c_double), c_int]
        getattr(self.shared_lib, f"{self.name}_acados_update_params_sparse").restype = c_int

        getattr(self.shared_lib, f"{self.name}_acados_update_params").argtypes = [c_void_p, c_int, POINTER(c_double), c_int]
        getattr(self.shared_lib, f"{self.name}_acados_update_params").restype = c_int

        getattr(self.shared_lib, f"{self.name}_acados_set_p_global_and_precompute_dependencies").argtypes = [c_void_p, POINTER(c_double), c_int]
        getattr(self.shared_lib, f"{self.name}_acados_set_p_global_and_precompute_dependencies").restype = c_int

        # these do not work for multi phase OCPs
        if isinstance(self.acados_ocp, AcadosOcp):
            getattr(self.shared_lib, f'{self.name}_acados_update_qp_solver_cond_N').argtypes = [c_void_p, c_int]
            getattr(self.shared_lib, f'{self.name}_acados_update_qp_solver_cond_N').restype = c_int
            getattr(self.shared_lib, f"{self.name}_acados_update_time_steps").argtypes = [c_void_p, c_int, c_void_p]
            getattr(self.shared_lib, f"{self.name}_acados_update_time_steps").restype = c_int

        return

    def __get_pointers_solver(self):
        """
        Private function to get the pointers for solver
        """
        # get pointers solver
        getattr(self.shared_lib, f"{self.name}_acados_get_nlp_opts").argtypes = [c_void_p]
        getattr(self.shared_lib, f"{self.name}_acados_get_nlp_opts").restype = c_void_p
        self.nlp_opts = getattr(self.shared_lib, f"{self.name}_acados_get_nlp_opts")(self.capsule)

        getattr(self.shared_lib, f"{self.name}_acados_get_nlp_dims").argtypes = [c_void_p]
        getattr(self.shared_lib, f"{self.name}_acados_get_nlp_dims").restype = c_void_p
        self.nlp_dims = getattr(self.shared_lib, f"{self.name}_acados_get_nlp_dims")(self.capsule)

        getattr(self.shared_lib, f"{self.name}_acados_get_nlp_config").argtypes = [c_void_p]
        getattr(self.shared_lib, f"{self.name}_acados_get_nlp_config").restype = c_void_p
        self.nlp_config = getattr(self.shared_lib, f"{self.name}_acados_get_nlp_config")(self.capsule)

        getattr(self.shared_lib, f"{self.name}_acados_get_nlp_out").argtypes = [c_void_p]
        getattr(self.shared_lib, f"{self.name}_acados_get_nlp_out").restype = c_void_p
        self.nlp_out = getattr(self.shared_lib, f"{self.name}_acados_get_nlp_out")(self.capsule)

        getattr(self.shared_lib, f"{self.name}_acados_get_sens_out").argtypes = [c_void_p]
        getattr(self.shared_lib, f"{self.name}_acados_get_sens_out").restype = c_void_p
        self.sens_out = getattr(self.shared_lib, f"{self.name}_acados_get_sens_out")(self.capsule)

        getattr(self.shared_lib, f"{self.name}_acados_get_nlp_in").argtypes = [c_void_p]
        getattr(self.shared_lib, f"{self.name}_acados_get_nlp_in").restype = c_void_p
        self.nlp_in = getattr(self.shared_lib, f"{self.name}_acados_get_nlp_in")(self.capsule)

        getattr(self.shared_lib, f"{self.name}_acados_get_nlp_solver").argtypes = [c_void_p]
        getattr(self.shared_lib, f"{self.name}_acados_get_nlp_solver").restype = c_void_p
        self.nlp_solver = getattr(self.shared_lib, f"{self.name}_acados_get_nlp_solver")(self.capsule)


    def solve_for_x0(self, x0_bar, fail_on_nonzero_status=True, print_stats_on_failure=True):
        """
        Wrapper around `solve()` which sets initial state constraint, solves the OCP, and returns u0.
        """
        self.set(0, "lbx", x0_bar)
        self.set(0, "ubx", x0_bar)

        status = self.solve()

        if status != 0:
            if print_stats_on_failure:
                self.print_statistics()
            if fail_on_nonzero_status:
                raise Exception(f'acados acados_ocp_solver returned status {status}')
            elif print_stats_on_failure:
                print(f'Warning: acados acados_ocp_solver returned status {status}')

        u0 = self.get(0, "u")
        return u0


    def solve(self) -> int:
        """
        Solve the ocp with current input.

        :return: status of the solver
        """
        self.status = getattr(self.shared_lib, f"{self.name}_acados_solve")(self.capsule)

        return self.status


    def get_dim_flat(self, field: str):
        """
        Get dimension of flattened iterate.
        """
        if field not in ['x', 'u', 'z', 'pi', 'lam', 'sl', 'su', 'p']:
            raise Exception(f'AcadosOcpSolver.get_dim_flat(field={field}): \'{field}\' is an invalid argument.')

        return self.__acados_lib.ocp_nlp_dims_get_total_from_attr(self.nlp_config, self.nlp_dims, self.nlp_out, field.encode('utf-8'))


    def custom_update(self, data_: np.ndarray):
        """
        A custom function that can be implemented by a user to be called between solver calls.
        By default this does nothing.
        The idea is to have a convenient wrapper to do complex updates of parameters and numerical data efficiently in C,
        in a function that is compiled into the solver library and can be conveniently used in the Python environment.
        """
        data = np.ascontiguousarray(data_, dtype=np.float64)
        c_data = cast(data.ctypes.data, POINTER(c_double))
        data_len = len(data)

        status = getattr(self.shared_lib, f"{self.name}_acados_custom_update")(self.capsule, c_data, data_len)

        return status


    def reset(self, reset_qp_solver_mem=1):
        """
        Sets current iterate to all zeros.
        """
        getattr(self.shared_lib, f"{self.name}_acados_reset")(self.capsule, reset_qp_solver_mem)


    def set_new_time_steps(self, new_time_steps):
        """
        Set new time steps.
        Recreates the solver if N changes.

            :param new_time_steps: 1 dimensional np array of new time steps for the solver

            .. note:: This allows for different use-cases: either set a new size of time_steps or a new distribution of
                      the shooting nodes without changing the number, e.g., to reach a different final time. Both cases
                      do not require a new code export and compilation.
        """
        if isinstance(self.acados_ocp, AcadosMultiphaseOcp):
            raise Exception('This function can only be used for single phase OCPs!')

        # unlikely but still possible
        if not self.solver_created:
            raise Exception('Solver was not yet created!')

        # check if time steps really changed in value
        if np.array_equal(self.__solver_options['time_steps'], new_time_steps):
            return

        N = new_time_steps.size
        new_time_steps_data = cast(new_time_steps.ctypes.data, POINTER(c_double))

        # check if recreation of acados is necessary (no need to recreate acados if sizes are identical)
        if len(self.__solver_options['time_steps']) == N:
            assert getattr(self.shared_lib, f"{self.name}_acados_update_time_steps")(self.capsule, N, new_time_steps_data) == 0
        else:  # recreate the solver with the new time steps
            self.solver_created = False

            # delete old memory (analog to __del__)
            getattr(self.shared_lib, f"{self.name}_acados_free")(self.capsule)

            # create solver with new time steps
            assert getattr(self.shared_lib, f"{self.name}_acados_create_with_discretization")(self.capsule, N, new_time_steps_data) == 0

            self.solver_created = True

            # get pointers solver
            self.__get_pointers_solver()

        # store time_steps, N
        self.__solver_options['time_steps'] = new_time_steps
        self.N = N
        self.__solver_options['Tsim'] = self.__solver_options['time_steps'][0]


    def update_qp_solver_cond_N(self, qp_solver_cond_N: int):
        """
        Recreate solver with new value `qp_solver_cond_N` with a partial condensing QP solver.
        This function is relevant for code reuse, i.e., if either `set_new_time_steps(...)` is used or
        the influence of a different `qp_solver_cond_N` is studied without code export and compilation.

            :param qp_solver_cond_N: new number of condensing stages for the solver

            .. note:: This function can only be used in combination with a partial condensing QP solver.

            .. note:: After `set_new_time_steps(...)` is used and depending on the new number of time steps it might be
                      necessary to change `qp_solver_cond_N` as well (using this function), i.e., typically
                      `qp_solver_cond_N < N`.
        """
        # unlikely but still possible
        if not self.solver_created:
            raise Exception('Solver was not yet created!')
        if self.N < qp_solver_cond_N:
            raise Exception('Setting qp_solver_cond_N to be larger than N does not work!')
        if self.__solver_options['qp_solver_cond_N'] != qp_solver_cond_N:
            self.solver_created = False

            # recreate the solver
            assert getattr(self.shared_lib, f'{self.name}_acados_update_qp_solver_cond_N')(self.capsule, qp_solver_cond_N) == 0

            # store the new value
            self.__solver_options['qp_solver_cond_N'] = qp_solver_cond_N
            self.solver_created = True

            # get pointers solver
            self.__get_pointers_solver()


    def eval_and_get_optimal_value_gradient(self, with_respect_to: str = "initial_state") -> np.ndarray:
        """
        Returns the gradient of the optimal value function w.r.t. what is specified in `with_respect_to`.

        Disclaimer: This function only returns reasonable values if the solver has converged for the current problem instance.

        Notes:
        - for field `initial_state`, the gradient is the Lagrange multiplier of the initial state constraint.
        The gradient computation consists of adding the Lagrange multipliers corresponding to the upper and lower bound of the initial state.

        - for field `initial_control`, the gradient is the Lagrange multiplier of the initial control constraint.
        The gradient computation consists of adding the Lagrange multipliers corresponding to the upper and lower bound of the initial control.
        This requires the OCP to have control bounds with lbu = ubu at the first stage, i.e. the gradient of the state-action value function or Q-function is computed.

        - for field `p_global`, the gradient of the Lagrange function w.r.t. the global parameters is computed.

        :param with_respect_to: string in ["initial_state", "initial_control", "p_global"]
        """

        if with_respect_to == "params_global":
            print("Deprecation warning: 'params_global' is deprecated and has been renamed to 'p_global'.")
            with_respect_to = "p_global"

        if with_respect_to == "initial_state":
            if not self.acados_ocp.constraints.has_x0:
                raise Exception("OCP does not have an initial state constraint.")

            nx = self.__acados_lib.ocp_nlp_dims_get_from_attr(self.nlp_config, self.nlp_dims, self.nlp_out, 0, "x".encode('utf-8'))
            nbu = self.__acados_lib.ocp_nlp_dims_get_from_attr(self.nlp_config, self.nlp_dims, self.nlp_out, 0, "lbu".encode('utf-8'))
            ns = self.__acados_lib.ocp_nlp_dims_get_from_attr(self.nlp_config, self.nlp_dims, self.nlp_out, 0, "s".encode('utf-8'))

            lam = self.get(0, 'lam')
            nlam_non_slack = lam.shape[0]//2 - ns
            grad = lam[nbu:nbu+nx] - lam[nlam_non_slack+nbu : nlam_non_slack+nbu+nx]

        elif with_respect_to == "initial_control":
            nu = self.__acados_lib.ocp_nlp_dims_get_from_attr(self.nlp_config, self.nlp_dims, self.nlp_out, 0, "u".encode('utf-8'))
            nbu = self.__acados_lib.ocp_nlp_dims_get_from_attr(self.nlp_config, self.nlp_dims, self.nlp_out, 0, "lbu".encode('utf-8'))
            ns = self.__acados_lib.ocp_nlp_dims_get_from_attr(self.nlp_config, self.nlp_dims, self.nlp_out, 0, "s".encode('utf-8'))
            lbu = self.get_from_qp_in(0, 'lbu')
            ubu = self.get_from_qp_in(0, 'ubu')

            if not (nbu == nu and np.all(lbu == ubu) and self.acados_ocp.dims.nsbu == 0):
                raise Exception("OCP does not have an equality constraint on the initial control.")

            lam = self.get(0, 'lam')
            nlam_non_slack = lam.shape[0]//2 - ns
            grad = lam[:nbu] - lam[nlam_non_slack : nlam_non_slack+nbu]

        elif with_respect_to == "p_global":
            np_global = self.__acados_lib.ocp_nlp_dims_get_from_attr(self.nlp_config, self.nlp_dims, self.nlp_out, 0, "p_global".encode('utf-8'))

            field = "p_global".encode('utf-8')
            t0 = time.time()
            grad = np.zeros((np_global,), dtype=np.float64, order='C')
            c_grad_p = cast(grad.ctypes.data, POINTER(c_double))
            self.__acados_lib.ocp_nlp_eval_lagrange_grad_p(self.nlp_solver, self.nlp_in, field, c_grad_p)
            self.time_value_grad = time.time() - t0

        else:
            raise Exception(f"AcadosOcpSolver.eval_and_get_optimal_value_gradient(): Unknown field: with_respect_to = {with_respect_to}")
        return grad


    def get_optimal_value_gradient(self, with_respect_to: str = "initial_state") -> np.ndarray:
        print("Deprecation warning: get_optimal_value_gradient() is deprecated and has been renamed to eval_and_get_optimal_value_gradient().")
        return self.eval_and_get_optimal_value_gradient(with_respect_to)


    def _sanity_check_solution_sensitivities(self, parametric=True) -> None:
        if not (self.acados_ocp.solver_options.qp_solver == 'FULL_CONDENSING_HPIPM' or
                self.acados_ocp.solver_options.qp_solver == 'PARTIAL_CONDENSING_HPIPM'):
            raise Exception("Parametric sensitivities are only available with HPIPM as QP solver.")

        if not (
            self.acados_ocp.solver_options.hessian_approx == 'EXACT' and
            self.acados_ocp.solver_options.regularize_method == 'NO_REGULARIZE' and
            self.acados_ocp.solver_options.levenberg_marquardt == 0 and
            self.acados_ocp.solver_options.exact_hess_constr == 1 and
            self.acados_ocp.solver_options.exact_hess_cost == 1 and
            self.acados_ocp.solver_options.exact_hess_dyn == 1 and
            self.acados_ocp.solver_options.fixed_hess == 0 and
            is_empty(self.acados_ocp.model.cost_expr_ext_cost_custom_hess_0) and
            is_empty(self.acados_ocp.model.cost_expr_ext_cost_custom_hess) and
            is_empty(self.acados_ocp.model.cost_expr_ext_cost_custom_hess_e)
        ):
            raise Exception("Parametric sensitivities are only correct if an exact Hessian is used!")

        if parametric and not self.acados_ocp.solver_options.with_solution_sens_wrt_params:
            raise Exception("Parametric sensitivities are only available if with_solution_sens_wrt_params is set to True.")


    def eval_solution_sensitivity(self,
                                  stages: Union[int, List[int]],
                                  with_respect_to: str,
                                  return_sens_x: bool = True,
                                  return_sens_u: bool = True,
                                  return_sens_pi: bool = False,
                                  return_sens_lam: bool = False,
                                  return_sens_su: bool = False,
                                  return_sens_sl: bool = False,
                                  ) \
                -> Dict:
        """
        Evaluate the sensitivity of the current solution x_i, u_i with respect to the initial state or the parameters for all stages i in `stages`.

            :param stages: stages for which the sensitivities are returned, int or list of int
            :param with_respect_to: string in ["initial_state", "p_global"]
            :param return_sens_x: Flag indicating whether sensitivities of x should be returned. Default: True.
            :param return_sens_u: Flag indicating whether sensitivities of u should be returned. Default: True.
            :param return_sens_pi: Flag indicating whether sensitivities of pi should be returned. Default: False.
            :param return_sens_lam: Flag indicating whether sensitivities of lam should be returned. Default: False.
            :param return_sens_su: Flag indicating whether sensitivities of su should be returned. Default: False.
            :param return_sens_sl: Flag indicating whether sensitivities of sl should be returned. Default: False.
            :returns: A dictionary with the solution sensitivities with fields sens_x, sens_u, sens_pi, sens_lam, sens_su, sens_sl if corresponding flags were set.
                    If stages is a list, sens_x, sens_lam, sens_su, sens_sl is a list of the same length.
                    For sens_u, sens_pi, the list has length len(stages) or len(stages)-1 depending on whether N is included or not.
                    If stages is a scalar, the returned sensitivities are np.ndarrays of shape (nfield[stages], ngrad).

        .. note::  Correct computation of sensitivities requires \n

        (1) HPIPM as QP solver, \n

        (2) the usage of an exact Hessian, \n

        (3) positive definiteness of the full-space Hessian if the square-root version of the Riccati recursion is used
            OR positive definiteness of the reduced Hessian if the classic Riccati recursion is used (compare: `solver_options.qp_solver_ric_alg`), \n

        (4) the solution of at least one QP in advance to evaluation of the sensitivities as the factorization is reused.

        .. note:: Timing of the sensitivities computation consists of time_solution_sens_lin, time_solution_sens_solve.
        .. note:: Solution sensitivities with respect to parameters are currently implemented assuming the parameter vector p is global within the OCP, i.e. p=p_i with i=0, ..., N.
        .. note:: Solution sensitivities with respect to parameters are currently implemented only for parametric discrete dynamics and parametric external costs (in particular, parametric constraints are not covered).
        """

        if with_respect_to == "params_global":
            print("Deprecation warning: 'params_global' is deprecated and has been renamed to 'p_global'.")
            with_respect_to = "p_global"

        stages_is_list = isinstance(stages, list)
        stages_ = stages if stages_is_list else [stages]

        sens_x = []
        sens_u = []
        sens_pi = []
        sens_lam = []
        sens_sl = []
        sens_su = []

        N = self.acados_ocp.solver_options.N_horizon

        for s in stages_:
            if not isinstance(s, int) or s < 0 or s > N:
                raise Exception(f"AcadosOcpSolver.eval_solution_sensitivity(): stages need to be int or list[int] and in [0, N], got stages = {stages_}.")

        if with_respect_to == "initial_state":
            nx = self.__acados_lib.ocp_nlp_dims_get_from_attr(self.nlp_config, self.nlp_dims, self.nlp_out, 0, "x".encode('utf-8'))
            ngrad = nx
            field = "ex"
            self._sanity_check_solution_sensitivities(parametric=False)

        elif with_respect_to == "p_global":
            np_global = self.__acados_lib.ocp_nlp_dims_get_from_attr(self.nlp_config, self.nlp_dims, self.nlp_out, 0, "p_global".encode('utf-8'))
            ngrad = np_global
            field = "p_global"
            self._sanity_check_solution_sensitivities()

            # compute jacobians wrt params in all modules
            t0 = time.time()
            self.__acados_lib.ocp_nlp_eval_params_jac(self.nlp_solver, self.nlp_in, self.nlp_out)
            self.time_solution_sens_lin = time.time() - t0

        else:
            raise Exception(f"AcadosOcpSolver.eval_solution_sensitivity(): Unknown field: with_respect_to = {with_respect_to}")

        # initialize jacobians with zeros
        for s in stages_:

            if return_sens_x:
                nx = self.__acados_lib.ocp_nlp_dims_get_from_attr(self.nlp_config, self.nlp_dims, self.nlp_out, s, "x".encode('utf-8'))
                sens_x.append(np.zeros((nx, ngrad)))

            if return_sens_lam:
                nlam = self.__acados_lib.ocp_nlp_dims_get_from_attr(self.nlp_config, self.nlp_dims, self.nlp_out, s, "lam".encode('utf-8'))
                sens_lam.append(np.zeros((nlam, ngrad)))

            if return_sens_sl:
                ns = self.__acados_lib.ocp_nlp_dims_get_from_attr(self.nlp_config, self.nlp_dims, self.nlp_out, s, "s".encode('utf-8'))
                sens_sl.append(np.zeros((ns, ngrad)))

            if return_sens_su:
                ns = self.__acados_lib.ocp_nlp_dims_get_from_attr(self.nlp_config, self.nlp_dims, self.nlp_out, s, "s".encode('utf-8'))
                sens_su.append(np.zeros((ns, ngrad)))

            if s < N:
                if return_sens_u:
                    nu = self.__acados_lib.ocp_nlp_dims_get_from_attr(self.nlp_config, self.nlp_dims, self.nlp_out, s, "u".encode('utf-8'))
                    sens_u.append(np.zeros((nu, ngrad)))

                if return_sens_pi:
                    npi = self.__acados_lib.ocp_nlp_dims_get_from_attr(self.nlp_config, self.nlp_dims, self.nlp_out, s, "pi".encode('utf-8'))
                    sens_pi.append(np.zeros((npi, ngrad)))

        self.time_solution_sens_solve = 0.0
        for k in range(ngrad):
            # evaluate sensitivity
            self.__acados_lib.ocp_nlp_eval_param_sens(self.nlp_solver, field.encode('utf-8'), 0, k, self.sens_out)

            # get timing
            self.time_solution_sens_solve += self.get_stats("time_solution_sensitivities")

            # extract sensitivities
            for n, s in enumerate(stages_):

                if return_sens_x:
                    sens_x[n][:, k] = self.get(s, "sens_x")
                if return_sens_lam:
                    sens_lam[n][:, k] = self.get(s, "sens_lam")
                if return_sens_sl:
                    sens_sl[n][:, k] = self.get(s, "sens_sl")
                if return_sens_su:
                    sens_su[n][:, k] = self.get(s, "sens_su")

                if s < N:
                    if return_sens_u:
                        sens_u[n][:, k] = self.get(s, "sens_u")
                    if return_sens_pi:
                        sens_pi[n][:, k] = self.get(s, "sens_pi")

        out = {}

        if return_sens_x:
            out["sens_x"] = sens_x if stages_is_list else sens_x[0]

        if return_sens_u:
            out["sens_u"] = sens_u if stages_is_list else sens_u[0]

        if return_sens_pi:
            out["sens_pi"] = sens_pi if stages_is_list else sens_pi[0]

        if return_sens_lam:
            out["sens_lam"] = sens_lam if stages_is_list else sens_lam[0]

        if return_sens_sl:
            out["sens_sl"] = sens_sl if stages_is_list else sens_sl[0]

        if return_sens_su:
            out["sens_su"] = sens_su if stages_is_list else sens_su[0]

        return out


    def eval_adjoint_solution_sensitivity(self,
                                          seed_x: Optional[Sequence[Tuple[int, np.ndarray]]],
                                          seed_u: Optional[Sequence[Tuple[int, np.ndarray]]],
                                          with_respect_to: str = "p_global",
                                          sanity_checks: bool = True,
                                          ) -> np.ndarray:
        """
        Evaluate the adjoint sensitivity of the solution with respect to the parameters.
            :param seed_x : Sequence of tuples of the form (stage: int, seed_vec: np.ndarray).
                    The stage is the stage at which the seed_vec is applied, and seed_vec is the seed for the states at that stage with shape (nx, n_seeds)
            :param seed_u : Sequence of tuples of the form (stage: int, seed_vec: np.ndarray).
                    The stage is the stage at which the seed_vec is applied, and seed_vec is the seed for the controls at that stage with shape (nu, n_seeds).
            :param with_respect_to : string in ["p_global"]
            :param sanity_checks : bool - whether to perform sanity checks, turn off for minimal overhead, default: True
        """

        # get n_seeds
        if seed_x is None:
            seed_x = []
        elif not isinstance(seed_x, Sequence):
            raise Exception(f"seed_x should be a Sequence, got {type(seed_x)}")

        if seed_u is None:
            seed_u = []
        elif not isinstance(seed_u, Sequence):
            raise Exception(f"seed_u should be a Sequence, got {type(seed_u)}")

        if len(seed_x) == 0 and len(seed_u) == 0:
            raise Exception("seed_x and seed_u cannot both be empty.")
        if len(seed_x) > 0:
            if not isinstance(seed_x[0], tuple) or len(seed_x[0]) != 2:
                raise Exception(f"seed_x[0] should be tuple of length 2, got seed_x[0] {seed_x[0]}")
            s = seed_x[0][1]
            if not isinstance(s, np.ndarray):
                raise Exception(f"seed_x[0][1] should be np.ndarray, got {type(s)}")
            n_seeds = seed_x[0][1].shape[1]
        if len(seed_u) > 0:
            if not isinstance(seed_u[0], tuple) or len(seed_u[0]) != 2:
                raise Exception(f"seed_u[0] should be tuple of length 2, got seed_u[0] {seed_u[0]}")
            s = seed_u[0][1]
            if not isinstance(s, np.ndarray):
                raise Exception(f"seed_u[0][1] should be np.ndarray, got {type(s)}")
            n_seeds = seed_u[0][1].shape[1]

        if sanity_checks:
            N_horizon = self.acados_ocp.solver_options.N_horizon
            self._sanity_check_solution_sensitivities()
            nx = self.__acados_lib.ocp_nlp_dims_get_from_attr(self.nlp_config, self.nlp_dims, self.nlp_out, 0, "x".encode('utf-8'))
            nu = self.__acados_lib.ocp_nlp_dims_get_from_attr(self.nlp_config, self.nlp_dims, self.nlp_out, 0, "u".encode('utf-8'))

            # check seeds
            for seed, name, dim in [(seed_x, "seed_x", nx), (seed_u, "seed_u", nu)]:
                for stage, seed_stage in seed:
                    if not isinstance(stage, int) or stage < 0 or stage > N_horizon:
                        raise Exception(f"AcadosOcpSolver.eval_adjoint_solution_sensitivity(): stage {stage} for {name} is not valid.")
                    if not isinstance(seed_stage, np.ndarray):
                        raise Exception(f"{name} for stage {stage} should be np.ndarray, got {type(seed_stage)}")
                    if seed_stage.shape != (dim, n_seeds):
                        raise Exception(f"{name} for stage {stage} should have shape (dim, n_seeds) = ({dim}, {n_seeds}), got {seed_stage.shape}.")

        if with_respect_to == "p_global":
            field = "p_global".encode('utf-8')

            nparam = self.__acados_lib.ocp_nlp_dims_get_from_attr(self.nlp_config, self.nlp_dims, self.nlp_out, 0, field)

            grad_p = np.zeros((n_seeds, nparam), order='C', dtype=np.float64)

            # compute jacobian wrt params
            t0 = time.time()
            self.__acados_lib.ocp_nlp_eval_params_jac(self.nlp_solver, self.nlp_in, self.nlp_out)
            self.time_solution_sens_lin = time.time() - t0

            self.time_solution_sens_solve = 0.0
            for i_seed in range(n_seeds):
                # set seed:
                self.reset_sens_out()
                for (stage, sx) in seed_x:
                    self.set(stage, 'sens_x', sx[:, i_seed])
                for (stage, su) in seed_u:
                    self.set(stage, 'sens_u', su[:, i_seed])

                c_grad_p = cast(grad_p[i_seed, :].ctypes.data, POINTER(c_double))

                # solve adjoint sensitivities
                self.__acados_lib.ocp_nlp_eval_solution_sens_adj_p(self.nlp_solver, self.nlp_in, self.sens_out, field, 0, c_grad_p)
                self.time_solution_sens_solve += self.get_stats("time_solution_sensitivities")

            return grad_p
        else:
            raise NotImplementedError(f"with_respect_to {with_respect_to} not implemented.")



    def eval_param_sens(self, index: int, stage: int=0, field="ex"):
        """
        Calculate the sensitivity of the current solution with respect to the initial state component of index.

        NOTE: Correct computation of sensitivities requires

        (1) HPIPM as QP solver,

        (2) the usage of an exact Hessian,

        (3) positive definiteness of the full-space Hessian if the square-root version of the Riccati recursion is used
            OR positive definiteness of the reduced Hessian if the classic Riccati recursion is used (compare: `solver_options.qp_solver_ric_alg`),
        (4) the solution of at least one QP in advance to evaluation of the sensitivities as the factorization is reused.

            :param index: integer corresponding to initial state index in range(nx)
        """

        print("WARNING: eval_param_sens() is deprecated. Please use eval_solution_sensitivity() instead!")

        self._sanity_check_solution_sensitivities(False)

        field = field.encode('utf-8')

        if not isinstance(index, int):
            raise Exception('AcadosOcpSolver.eval_param_sens(): index must be Integer.')

        if field == "ex":
            if not stage == 0:
                raise Exception('AcadosOcpSolver.eval_param_sens(): only stage == 0 is supported.')
            nx = self.__acados_lib.ocp_nlp_dims_get_from_attr(self.nlp_config, self.nlp_dims, self.nlp_out, stage, "x".encode('utf-8'))

            if index < 0 or index > nx:
                raise Exception(f'AcadosOcpSolver.eval_param_sens(): index must be in [0, nx-1], got: {index}.')

        elif field == "p_global":
            nparam = self.__acados_lib.ocp_nlp_dims_get_from_attr(self.nlp_config, self.nlp_dims, self.nlp_out, 0, "p".encode('utf-8'))

            if index < 0 or index > nparam:
                raise Exception(f'AcadosOcpSolver.eval_param_sens(): index must be in [0, nparam-1], got: {index}.')

        # actual eval_param
        self.__acados_lib.ocp_nlp_eval_param_sens(self.nlp_solver, field, stage, index, self.sens_out)

        return


    def get(self, stage_: int, field_: str):
        """
        Get the last solution of the solver:

            :param stage: integer corresponding to shooting node
            :param field: string in ['x', 'u', 'z', 'pi', 'lam', 'sl', 'su', 'p', 'sens_u', 'sens_pi', 'sens_x', 'sens_lam', 'sens_sl', 'sens_su']

            .. note:: regarding lam: \n
                    the inequalities are internally organized in the following order: \n
                    [ lbu lbx lg lh lphi ubu ubx ug uh uphi; \n
                      lsbu lsbx lsg lsh lsphi usbu usbx usg ush usphi]

            .. note:: pi: multipliers for dynamics equality constraints \n
                      lam: multipliers for inequalities \n
                      t: slack variables corresponding to evaluation of all inequalities (at the solution) \n
                      sl: slack variables of soft lower inequality constraints \n
                      su: slack variables of soft upper inequality constraints \n
        """

        out_fields = ['x', 'u', 'z', 'pi', 'lam', 'sl', 'su']
        in_fields = ['p']
        sens_fields = ['sens_u', 'sens_x', 'sens_pi', 'sens_lam', 'sens_sl', 'sens_su']
        all_fields = out_fields + in_fields + sens_fields

        if (field_ not in all_fields):
            raise Exception(f'AcadosOcpSolver.get(stage={stage_}, field={field_}): \'{field_}\' is an invalid argument.\
                    \n Possible values are {all_fields}.')

        if not isinstance(stage_, int):
            raise Exception(f'AcadosOcpSolver.get(stage={stage_}, field={field_}): stage index must be an integer, got type {type(stage_)}.')

        if stage_ < 0 or stage_ > self.N:
            raise Exception(f'AcadosOcpSolver.get(stage={stage_}, field={field_}): stage index must be in [0, {self.N}], got: {stage_}.')

        if stage_ == self.N and field_ == 'pi':
            raise Exception(f'AcadosOcpSolver.get(stage={stage_}, field={field_}): field \'{field_}\' does not exist at final stage {stage_}.')

        field = field_.replace('sens_', '') if field_ in sens_fields else field_
        field = field.encode('utf-8')

        dims = self.__acados_lib.ocp_nlp_dims_get_from_attr(self.nlp_config, self.nlp_dims, self.nlp_out, stage_, field)

        out = np.zeros((dims,), dtype=np.float64, order="C")
        out_data = cast(out.ctypes.data, POINTER(c_double))

        if field_ in in_fields:
            self.__acados_lib.ocp_nlp_in_get(self.nlp_config, self.nlp_dims, self.nlp_in, stage_, field, out_data)
        else:
            out_pointer = self.nlp_out if field_ in out_fields else self.sens_out
            self.__acados_lib.ocp_nlp_out_get(self.nlp_config, self.nlp_dims, out_pointer, stage_, field, out_data)

        return out


    def get_flat(self, field_: str) -> np.ndarray:
        """
        Get concatenation of all stages of last solution of the solver.

            :param field: string in ['x', 'u', 'z', 'pi', 'lam', 'sl', 'su', 'p', 'p_global']

            .. note:: The parameter 'p_global' has no stage-wise structure and is processed in a memory saving manner by default. \n
                  In order to read the 'p_global' parameter, the option 'save_p_global' must be set to 'True' upon instantiation. \n
        """
        if field_ not in ['x', 'u', 'z', 'pi', 'lam', 'sl', 'su', 'p', 'p_global']:
            raise Exception(f'AcadosOcpSolver.get_flat(field={field_}): \'{field_}\' is an invalid argument.')

        if field_ == 'p_global':
            if not self.__save_p_global:
                raise Exception(f'The field \'{field_}\' is not stored within the solver by default. Please set the option \'save_p_global=True\' when instantiating the solver.')
            return self.__p_global_values

        field = field_.encode('utf-8')

        dims = self.__acados_lib.ocp_nlp_dims_get_total_from_attr(self.nlp_config, self.nlp_dims, self.nlp_out, field)

        out = np.zeros((dims,), dtype=np.float64, order="C")
        out_data = cast(out.ctypes.data, POINTER(c_double))

        self.__acados_lib.ocp_nlp_get_all(self.nlp_solver, self.nlp_in, self.nlp_out, field, out_data)

        return out


    def set_flat(self, field_: str, value_: np.ndarray) -> None:
        """
        Set concatenation solver initialization .

            :param field: string in ['x', 'u', 'z', 'pi', 'lam', 'sl', 'su', 'p']
        """
        field = field_.encode('utf-8')
        if field_ not in ['x', 'u', 'z', 'pi', 'lam', 'sl', 'su', 'p']:
            raise Exception(f'AcadosOcpSolver.get_flat(field={field_}): \'{field_}\' is an invalid argument.')
        dims = self.__acados_lib.ocp_nlp_dims_get_total_from_attr(self.nlp_config, self.nlp_dims, self.nlp_out, field)

        if len(value_) != dims:
            raise Exception(f'AcadosOcpSolver.set_flat(field={field_}, value): value has wrong length, expected {dims}, got {len(value_)}.')

        value_ = value_.astype(float)
        value_data = cast(value_.ctypes.data, POINTER(c_double))
        value_data_p = cast((value_data), c_void_p)

        self.__acados_lib.ocp_nlp_set_all(self.nlp_solver, self.nlp_in, self.nlp_out, field, value_data_p)
        return


    def print_statistics(self):
        """
        prints statistics of previous solver run as a table:
            - iter: iteration number
            - res_stat: stationarity residual
            - res_eq: residual wrt equality constraints (dynamics)
            - res_ineq: residual wrt inequality constraints (constraints)
            - res_comp: residual wrt complementarity conditions
            - qp_stat: status of QP solver
            - qp_iter: number of QP iterations
            - alpha: SQP step size
            - qp_res_stat: stationarity residual of the last QP solution
            - qp_res_eq: residual wrt equality constraints (dynamics) of the last QP solution
            - qp_res_ineq: residual wrt inequality constraints (constraints)  of the last QP solution
            - qp_res_comp: residual wrt complementarity conditions of the last QP solution
        """
        stat = self.get_stats("statistics")

        if self.__solver_options['nlp_solver_type'] == 'SQP':
            print('\niter\tres_stat\tres_eq\t\tres_ineq\tres_comp\tqp_stat\tqp_iter\talpha')
            if stat.shape[0]>8:
                print('\tqp_res_stat\tqp_res_eq\tqp_res_ineq\tqp_res_comp')
            for jj in range(stat.shape[1]):
                print(f'{int(stat[0][jj]):d}\t{stat[1][jj]:e}\t{stat[2][jj]:e}\t{stat[3][jj]:e}\t' +
                      f'{stat[4][jj]:e}\t{int(stat[5][jj]):d}\t{int(stat[6][jj]):d}\t{stat[7][jj]:e}\t')
                if stat.shape[0]>8:
                    print('\t{:e}\t{:e}\t{:e}\t{:e}'.format( \
                        stat[8][jj], stat[9][jj], stat[10][jj], stat[11][jj]))
            print('\n')
        elif self.__solver_options['nlp_solver_type'] == 'SQP_RTI':
            header = '\niter\tqp_stat\tqp_iter'
            if self.__solver_options['nlp_solver_ext_qp_res'] == 1:
                header += '\tqp_res_stat\tqp_res_eq\tqp_res_ineq\tqp_res_comp'
            if self.__solver_options['rti_log_residuals'] == 1:
                header += '\tres_stat\tres_eq\t\tres_ineq\tres_comp'
            print(header)
            for jj in range(stat.shape[1]):
                line = '{:d}\t{:d}\t{:d}'.format( int(stat[0][jj]), int(stat[1][jj]), int(stat[2][jj]))
                offset = 2
                if self.__solver_options['nlp_solver_ext_qp_res'] == 1:
                    line += '\t{:e}\t{:e}\t{:e}\t{:e}'.format( \
                         stat[offset+1][jj], stat[offset+2][jj], stat[offset+3][jj], stat[offset+4][jj])
                    offset += 4
                if self.__solver_options['rti_log_residuals'] == 1:
                    line += '\t{:e}\t{:e}\t{:e}\t{:e}'.format( \
                         stat[offset+1][jj], stat[offset+2][jj], stat[offset+3][jj], stat[offset+4][jj])
                print(line)
            print('\n')
        elif self.__solver_options['nlp_solver_type'] == 'DDP':
            for jj in range(stat.shape[1]):
                if jj % 10 == 0:
                    # print('\niter\tres_stat\tres_eq\t\tqp_stat\tqp_iter\talpha')
                    print(("{iter:>6} | {obj:^10} | {inf:^10} | {stat:^10} | "
                   "{alpha:^10} | {gamma:^10} | {qp_status:^10} | {qp_iter:^10}").format(
                        obj='objective',
                        iter='iter.',
                        inf='res_eq',
                        stat='res_stat',
                        alpha='alpha',
                        gamma='LM_reg.',
                        qp_status='qp_status',
                        qp_iter='qp_iter.'))
                # print(f'{int(stat[0][jj]):d}\t{stat[1][jj]:e}\t{stat[2][jj]:e}\t{int(stat[5][jj]):d}\t{int(stat[6][jj]):d}\t{stat[7][jj]:e}\t')
                print(("{iter:>6} | {obj:^10.4e} | {inf:^10.4e} | {stat:^10.4e} | "
                   "{alpha:^10.4e} | {gamma:^10.4e} | {qp_status:^10} | {qp_iter:^10}").format(
                     iter=int(stat[0][jj]),
                     stat=stat[1][jj],
                     inf=stat[2][jj],
                     obj=stat[3][jj],
                     gamma=stat[4][jj],
                     qp_status=int(stat[5][jj]),
                     qp_iter=int(stat[6][jj]),
                     alpha=stat[7][jj]))
            print('\n')

        return


    def store_iterate(self, filename: str = '', overwrite: bool = False, verbose: bool = True):
        """
        Stores the current iterate of the OCP solver in a json file.
        Note: This does not contain the iterate of the integrators, and the parameters.

            :param filename: if not set, use f'{self.name}_iterate.json'
            :param overwrite: if false and filename exists add timestamp to filename
        """
        if filename == '':
            filename = f'{self.name}_iterate.json'

        if not overwrite:
            # append timestamp
            if os.path.isfile(filename):
                filename = filename[:-5]
                filename += datetime.now().strftime('%Y-%m-%d-%H:%M:%S.%f') + '.json'

        # get iterate:
        solution = dict()

        lN = len(str(self.N+1))
        for i in range(self.N+1):
            i_string = f'{i:0{lN}d}'
            solution['x_'+i_string] = self.get(i,'x')
            solution['u_'+i_string] = self.get(i,'u')
            solution['z_'+i_string] = self.get(i,'z')
            solution['lam_'+i_string] = self.get(i,'lam')
            solution['sl_'+i_string] = self.get(i, 'sl')
            solution['su_'+i_string] = self.get(i, 'su')
            if i < self.N:
                solution['pi_'+i_string] = self.get(i,'pi')

        for k in list(solution.keys()):
            if len(solution[k]) == 0:
                del solution[k]

        # save
        with open(filename, 'w') as f:
            json.dump(solution, f, default=make_object_json_dumpable, indent=4, sort_keys=True)

        if verbose:
            print("stored current iterate in ", os.path.join(os.getcwd(), filename))


    def qp_diagnostics(self, hessian_type: str = 'FULL_HESSIAN'):
            """
            Compute some diagnostic values for the last QP.
            result = ocp_solver.qp_diagnostics(hessian_type). Possible values are
            'FULL_HESSIAN' or 'PROJECTED_HESSIAN'

            returns a dictionary with the following fields:
            - min_eigv_stage: dict with minimum eigenvalue for each Hessian block.
            - max_eigv_stage: dict with maximum eigenvalue for each Hessian block.
            - condition_number_stage: dict with condition number for each Hessian block.
            - condition_number_global: condition number for the full Hessian.
            - min_eigv_global: minimum eigenvalue for the full Hessian.
            - min_abs_eigv_global: minimum absolute eigenvalue for the full Hessian.
            - max_eigv_global: maximum eigenvalue for the full Hessian.

            for the 'PROJECTED_HESSIAN' it also includes
            - min_eigv_P_global: minimum eigenvalue of P matrices
            - min_abs_eigv_P_global: minimum absolute eigenvalue of P matrices
            """
            if hessian_type not in ['FULL_HESSIAN', 'PROJECTED_HESSIAN']:
                raise TypeError("Input should be string with value FULL_HESSIAN, PROJECTED_HESSIAN")

            qp_diagnostic = {}
            N_horizon = self.N

            min_eigv_global = np.inf
            max_eigv_global = -np.inf
            min_abs_eigv = np.inf
            max_abs_eigv = -np.inf

            min_eig_P_global = np.inf
            min_abs_eig_P_global = np.inf

            max_eigv_stage = []
            min_eigv_stage = []
            condition_number_stage = []

            for i in range(N_horizon+1):
                if hessian_type == "FULL_HESSIAN":
                    hess_block = self.get_hessian_block(i)

                elif hessian_type == "PROJECTED_HESSIAN":
                    P_mat = self.get_from_qp_in(i, 'P')
                    B_mat = self.get_from_qp_in(i-1, 'B')
                    # Lr: lower triangular decomposition of R within Riccati != R in qp_in!
                    Lr = self.get_from_qp_in(i-1, 'Lr')
                    R_ric = Lr @ Lr.T
                    hess_block = R_ric + B_mat.T @ P_mat @ B_mat

                    # P
                    eigv = np.linalg.eigvals(P_mat)
                    min_eig_P_global = min(min_eig_P_global, np.min(eigv))
                    min_abs_eig_P_global = min(min_abs_eig_P_global, np.min(np.abs(eigv)))

                else:
                    raise ValueError("Wrong input given to function! Possible inputs are FULL_HESSIAN, PROJECTED_HESSIAN")

                eigv = np.linalg.eigvals(hess_block)
                min_eigv = np.min(eigv)
                max_eigv = np.max(eigv)

                min_eigv_global = min(min_eigv, min_eigv_global)
                max_eigv_global = max(max_eigv, max_eigv_global)
                min_abs_eigv = min(min_abs_eigv, np.min(np.abs(eigv)))
                max_abs_eigv = max(max_abs_eigv, np.max(np.abs(eigv)))

                max_eigv_stage.append(max_eigv)
                min_eigv_stage.append(min_eigv)
                condition_number_stage.append(np.max(np.abs(eigv))/np.min(np.abs(eigv)))

            condition_number_global = max_abs_eigv/min_abs_eigv

            qp_diagnostic['max_eigv_global'] = max_eigv_global
            qp_diagnostic['min_eigv_global'] = min_eigv_global
            qp_diagnostic['min_abs_eigv_global'] = min_abs_eigv
            qp_diagnostic['condition_number_global'] = condition_number_global
            qp_diagnostic['max_eigv_stage'] = max_eigv_stage
            qp_diagnostic['min_eigv_stage'] = min_eigv_stage
            qp_diagnostic['condition_number_stage'] = condition_number_stage

            if hessian_type == "PROJECTED_HESSIAN":
                qp_diagnostic['min_eigv_P_global'] = min_eig_P_global
                qp_diagnostic['min_abs_eigv_P_global'] = min_abs_eig_P_global

            return qp_diagnostic


    def dump_last_qp_to_json(self, filename: str = '', overwrite=False):
        """
        Dumps the latest QP data into a json file

            :param filename: if not set, use name + timestamp + '.json'
            :param overwrite: if false and filename exists add timestamp to filename
        """
        if filename == '':
            filename = f'{self.name}_QP.json'

        if not overwrite:
            # append timestamp
            if os.path.isfile(filename):
                filename = filename[:-5]
                filename += datetime.now().strftime('%Y-%m-%d-%H:%M:%S.%f') + '.json'

        # get QP data:
        qp_data = dict()

        lN = len(str(self.N+1))
        for field in self.__qp_dynamics_fields:
            for i in range(self.N):
                qp_data[f'{field}_{i:0{lN}d}'] = self.get_from_qp_in(i,field)

        for field in self.__qp_constraint_fields + self.__qp_cost_fields + self.__qp_constraint_int_fields:
            for i in range(self.N+1):
                qp_data[f'{field}_{i:0{lN}d}'] = self.get_from_qp_in(i,field)

        # remove empty fields
        for k in list(qp_data.keys()):
            if len(qp_data[k]) == 0:
                del qp_data[k]

        # save
        with open(filename, 'w') as f:
            json.dump(qp_data, f, default=make_object_json_dumpable, indent=4, sort_keys=True)
        print("stored qp from solver memory in ", os.path.join(os.getcwd(), filename))



    def load_iterate(self, filename:str, verbose: bool = True):
        """
        Loads the iterate stored in json file with filename into the ocp solver.
        Note: This does not contain the iterate of the integrators, and the parameters.
        """
        if not os.path.isfile(filename):
            raise Exception('load_iterate: failed, file does not exist: ' + os.path.join(os.getcwd(), filename))

        with open(filename, 'r') as f:
            solution = json.load(f)

        if verbose:
            print(f"loading iterate {filename}")
        for key in solution.keys():
            (field, stage) = key.split('_')
            self.set(int(stage), field, np.array(solution[key]))


    def store_iterate_to_obj(self) -> AcadosOcpIterate:
        """
        Returns the current iterate of the OCP solver as an AcadosOcpIterate.
        """
        d = {}
        for field in ["x", "u", "z", "sl", "su", "pi", "lam"]:
            traj = []
            for n in range(self.N+1):
                if n < self.N or not (field in ["u", "pi", "z"]):
                    traj.append(self.get(n, field))

            d[f"{field}_traj"] = traj

        return AcadosOcpIterate(**d)


    def load_iterate_from_obj(self, iterate: AcadosOcpIterate):
        """
        Loads the provided iterate into the OCP solver.
        Note: The iterate object does not contain the the parameters.
        """

        for key, traj in iterate.__dict__.items():
            field = key.replace('_traj', '')

            for n, val in enumerate(traj):
                self.set(n, field, val)


    def store_iterate_to_flat_obj(self) -> AcadosOcpFlattenedIterate:
        """
        Returns the current iterate of the OCP solver as an AcadosOcpFlattenedIterate.
        """
        return AcadosOcpFlattenedIterate(x = self.get_flat("x"),
                                         u = self.get_flat("u"),
                                         z = self.get_flat("z"),
                                         sl = self.get_flat("sl"),
                                         su = self.get_flat("su"),
                                         pi = self.get_flat("pi"),
                                         lam = self.get_flat("lam"))

    def load_iterate_from_flat_obj(self, iterate: AcadosOcpFlattenedIterate) -> None:
        """
        Loads the provided iterate into the OCP solver.
        Note: The iterate object does not contain the the parameters.
        """
        self.set_flat("x", iterate.x)
        self.set_flat("u", iterate.u)
        self.set_flat("z", iterate.z)
        self.set_flat("sl", iterate.sl)
        self.set_flat("su", iterate.su)
        self.set_flat("pi", iterate.pi)
        self.set_flat("lam", iterate.lam)


    # TODO this should be a property
    def get_status(self) -> int:
        """
        Returns the status of the last solver call.

        Status codes:
            - 0: Success (ACADOS_SUCCESS)
            - 1: NaN detected (ACADOS_NAN_DETECTED)
            - 2: Maximum number of iterations reached (ACADOS_MAXITER)
            - 3: Minimum step size reached (ACADOS_MINSTEP)
            - 4: QP solver failed (ACADOS_QP_FAILURE)
            - 5: Solver created (ACADOS_READY)
            - 6: Problem unbounded (ACADOS_UNBOUNDED)
            - 7: Solver timeout (ACADOS_TIMEOUT)

        See `return_values` in https://github.com/acados/acados/blob/main/acados/utils/types.h
        """
        return self.status

    def get_stats(self, field_: str) -> Union[int, float, np.ndarray]:
        """
        Get the information of the last solver call.

            :param field: string in ['statistics', 'time_tot', 'time_lin', 'time_sim', 'time_sim_ad', 'time_sim_la', 'time_qp', 'time_qp_solver_call', 'time_reg', 'nlp_iter', 'sqp_iter', 'residuals', 'qp_iter', 'alpha']

        Available fileds:
            - time_tot: total CPU time previous call
            - time_lin: CPU time for linearization
            - time_sim: CPU time for integrator
            - time_sim_ad: CPU time for integrator contribution of external function calls
            - time_sim_la: CPU time for integrator contribution of linear algebra
            - time_qp: CPU time qp solution
            - time_qp_solver_call: CPU time inside qp solver (without converting the QP)
            - time_qp_xcond: time_glob: CPU time globalization
            - time_solution_sensitivities: CPU time for previous call to eval_param_sens
            - time_solution_sens_lin: CPU time for linearization in eval_param_sens
            - time_solution_sens_solve: CPU time for solving in eval_solution_sensitivity
            - time_reg: CPU time regularization
            - time_preparation: CPU time for last preparation phase, relevant for (AS-)RTI, zero otherwise
            - time_feedback: CPU time for last feedback phase, relevant for (AS-)RTI, otherwise returns total compuation time.
            - sqp_iter: number of SQP iterations
            - nlp_iter: number of NLP solver iterations (DDP or SQP)
            - qp_stat: status of QP solver
            - qp_iter: vector of QP iterations for last SQP call
            - statistics: table with info about last iteration
            - stat_m: number of rows in statistics matrix
            - stat_n: number of columns in statistics matrix
            - residuals: residuals of current iterate
            - alpha: step sizes of SQP iterations
        """

        if field_ == "time_solution_sens_lin":
            return self.time_solution_sens_lin
        elif field_ == "time_solution_sens_solve":
            return self.time_solution_sens_solve

        double_fields = ['time_tot',
                  'time_lin',
                  'time_sim',
                  'time_sim_ad',
                  'time_sim_la',
                  'time_qp',
                  'time_qp_solver_call',
                  'time_qp_xcond',
                  'time_glob',
                  'time_solution_sensitivities',
                  'time_reg',
                  'time_preparation',
                  'time_feedback',
        ]
        fields = double_fields + [
                  'sqp_iter',
                  'ddp_iter',
                  'nlp_iter',
                  'qp_stat',
                  'qp_iter',
                  'statistics',
                  'stat_m',
                  'stat_n',
                  'residuals',
                  'alpha',
                  'res_eq_all',
                  'res_stat_all',
                ]

        field = field_.encode('utf-8')

        if field_ in ['ddp_iter', 'sqp_iter', 'nlp_iter', 'stat_m', 'stat_n']:
            out = c_int(0)
            self.__acados_lib.ocp_nlp_get(self.nlp_solver, field, byref(out))
            return out.value

        elif field_ in double_fields:
            out = c_double(0)
            self.__acados_lib.ocp_nlp_get(self.nlp_solver, field, byref(out))
            return out.value

        elif field_ == 'statistics':
            nlp_iter = self.get_stats("nlp_iter")
            stat_m = self.get_stats("stat_m")
            stat_n = self.get_stats("stat_n")
            min_size = min([stat_m, nlp_iter+1])
            out = np.zeros((stat_n+1, min_size), dtype=np.float64, order="C")
            out_data = cast(out.ctypes.data, POINTER(c_double))
            self.__acados_lib.ocp_nlp_get(self.nlp_solver, field, out_data)
            return out

        elif field_ == 'primal_step_norm':
            nlp_iter = self.get_stats("nlp_iter")
            out = np.zeros((nlp_iter,), dtype=np.float64, order="C")
            out_data = cast(out.ctypes.data, POINTER(c_double))
            self.__acados_lib.ocp_nlp_get(self.nlp_solver, field, out_data)
            return out

        elif field_ == 'qp_stat':
            full_stats = self.get_stats('statistics')
            if self.__solver_options['nlp_solver_type'] == 'SQP':
                return full_stats[5, :]
            elif self.__solver_options['nlp_solver_type'] == 'SQP_RTI':
                return full_stats[1, :]

        elif field_ == 'qp_iter':
            full_stats = self.get_stats('statistics')
            if self.__solver_options['nlp_solver_type'] == 'SQP':
                return full_stats[6, :]
            elif self.__solver_options['nlp_solver_type'] == 'SQP_RTI':
                return full_stats[2, :]

        elif field_ == 'alpha':
            full_stats = self.get_stats('statistics')
            if self.__solver_options['nlp_solver_type'] == 'SQP':
                return full_stats[7, :]
            else: # self.__solver_options['nlp_solver_type'] == 'SQP_RTI':
                raise Exception("alpha values are not available for SQP_RTI")

        elif field_ == 'residuals':
            return self.get_residuals()

        elif field_ == 'res_eq_all':
            full_stats = self.get_stats('statistics')
            if self.__solver_options['nlp_solver_type'] == 'SQP':
                return full_stats[2, :]
            elif self.__solver_options['nlp_solver_type'] == 'SQP_RTI':
                if self.__solver_options['rti_log_residuals'] == 1:
                    return full_stats[4, :]
                else:
                    raise Exception("res_eq_all is not available for SQP_RTI if rti_log_residuals is not enabled.")
            else:
                raise Exception(f"res_eq_all is not available for nlp_solver_type {self.__solver_options['nlp_solver_type']}.")

        elif field_ == 'res_stat_all':
            full_stats = self.get_stats('statistics')
            if self.__solver_options['nlp_solver_type'] == 'SQP':
                return full_stats[1, :]
            elif self.__solver_options['nlp_solver_type'] == 'SQP_RTI':
                if self.__solver_options['rti_log_residuals'] == 1:
                    return full_stats[3, :]
                else:
                    raise Exception("res_stat_all is not available for SQP_RTI if rti_log_residuals is not enabled.")
            else:
                raise Exception(f"res_stat_all is not available for nlp_solver_type {self.__solver_options['nlp_solver_type']}.")

        elif field_ == 'res_ineq_all':
            full_stats = self.get_stats('statistics')
            if self.__solver_options['nlp_solver_type'] == 'SQP':
                return full_stats[3, :]
            elif self.__solver_options['nlp_solver_type'] == 'SQP_RTI':
                if self.__solver_options['rti_log_residuals'] == 1:
                    return full_stats[5, :]
                else:
                    raise Exception("res_ineq_all is not available for SQP_RTI if rti_log_residuals is not enabled.")
            else:
                raise Exception(f"res_ineq_all is not available for nlp_solver_type {self.__solver_options['nlp_solver_type']}.")

        elif field_ == 'res_comp_all':
            full_stats = self.get_stats('statistics')
            if self.__solver_options['nlp_solver_type'] == 'SQP':
                return full_stats[4, :]
            elif self.__solver_options['nlp_solver_type'] == 'SQP_RTI':
                if self.__solver_options['rti_log_residuals'] == 1:
                    return full_stats[6, :]
                else:
                    raise Exception("res_comp_all is not available for SQP_RTI if rti_log_residuals is not enabled.")
            else:
                raise Exception(f"res_comp_all is not available for nlp_solver_type {self.__solver_options['nlp_solver_type']}.")

        else:
            raise Exception(f'AcadosOcpSolver.get_stats(): \'{field}\' is not a valid argument.'
                    + f'\n Possible values are {fields}.')


    def get_cost(self) -> float:
        """
        Returns the cost value of the current solution.
        """
        # compute cost internally
        self.__acados_lib.ocp_nlp_eval_cost(self.nlp_solver, self.nlp_in, self.nlp_out)

        # create output array
        out = np.zeros((1,), dtype=np.float64, order="C")
        out_data = cast(out.ctypes.data, POINTER(c_double))

        # call getter
        field = "cost_value".encode('utf-8')
        self.__acados_lib.ocp_nlp_get(self.nlp_solver, field, out_data)

        return out[0]


    def get_residuals(self, recompute=False):
        """
        Returns an array of the form [res_stat, res_eq, res_ineq, res_comp].
        The residuals has to be computed for SQP_RTI solver, since it is not available by default.

        - res_stat: stationarity residual
        - res_eq: residual wrt equality constraints (dynamics)
        - res_ineq: residual wrt inequality constraints (constraints)
        - res_comp: residual wrt complementarity conditions
        """
        # compute residuals if RTI
        if self.__solver_options['nlp_solver_type'] == 'SQP_RTI' or recompute:
            self.__acados_lib.ocp_nlp_eval_residuals(self.nlp_solver, self.nlp_in, self.nlp_out)

        # create output array
        out = np.zeros((4, 1), dtype=np.float64, order="C")
        out_data = cast(out.ctypes.data, POINTER(c_double))

        # call getters
        field = "res_stat".encode('utf-8')
        self.__acados_lib.ocp_nlp_get(self.nlp_solver, field, out_data)

        out_data = cast(out[1].ctypes.data, POINTER(c_double))
        field = "res_eq".encode('utf-8')
        self.__acados_lib.ocp_nlp_get(self.nlp_solver, field, out_data)

        out_data = cast(out[2].ctypes.data, POINTER(c_double))
        field = "res_ineq".encode('utf-8')
        self.__acados_lib.ocp_nlp_get(self.nlp_solver, field, out_data)

        out_data = cast(out[3].ctypes.data, POINTER(c_double))
        field = "res_comp".encode('utf-8')
        self.__acados_lib.ocp_nlp_get(self.nlp_solver, field, out_data)
        return out.flatten()


    def get_initial_residuals(self) -> np.ndarray:
        """
        Returns an array of the form [res_stat, res_eq, res_ineq, res_comp].
        Residuals: residuals of initial iterate in previous solver call
        """
        full_stats = self.get_stats('statistics')
        if self.__solver_options['nlp_solver_type'] == 'SQP':
            return full_stats[1:5, 0]
        elif self.__solver_options['nlp_solver_type'] == 'SQP_RTI':
            if self.__solver_options['rti_log_residuals'] == 1:
                return full_stats[3:7, 0]
            else:
                raise Exception("initial_residuals is only available for SQP_RTI if rti_log_residuals is enabled, for efficiency the rti_log_only_available_residuals option is recommended.")
        else:
            raise Exception(f"initial_residuals is not available for nlp_solver_type {self.__solver_options['nlp_solver_type']}.")

    # Note: this function should not be used anymore, better use cost_set, constraints_set
    def set(self, stage_: int, field_: str, value_: np.ndarray):
        """
        Set numerical data inside the solver.

            :param stage: integer corresponding to shooting node
            :param field: string in ['x', 'u', 'pi', 'lam', 'p', 'xdot_guess', 'z_guess', 'sens_x', 'sens_u']

            .. note:: regarding lam: \n
                    the inequalities are internally organized in the following order: \n
                    [ lbu lbx lg lh lphi ubu ubx ug uh uphi; \n
                      lsbu lsbx lsg lsh lsphi usbu usbx usg ush usphi]

            .. note:: pi: multipliers for dynamics equality constraints \n
                      lam: multipliers for inequalities \n
                      t: slack variables corresponding to evaluation of all inequalities (at the solution) \n
                      sl: slack variables of soft lower inequality constraints \n
                      su: slack variables of soft upper inequality constraints \n
        """
        cost_fields = ['y_ref', 'yref']
        constraints_fields = ['lbx', 'ubx', 'lbu', 'ubu']
        out_fields = ['x', 'u', 'pi', 'lam', 'z', 'sl', 'su']
        mem_fields = ['xdot_guess', 'z_guess']
        sens_fields = ['sens_x', 'sens_u']

        if not isinstance(stage_, int):
            raise Exception('stage should be integer.')
        elif stage_ < 0 or stage_ > self.N:
            raise Exception(f'stage should be in [0, N], got {stage_}')

        # cast value_ to avoid conversion issues
        if isinstance(value_, (float, int)):
            value_ = np.array([value_])
        value_ = value_.astype(float)

        field = field_.replace("sens_", "").encode('utf-8')

        stage = c_int(stage_)

        # treat parameters separately
        if field_ == 'p':
            value_data = cast(value_.ctypes.data, POINTER(c_double))
            assert getattr(self.shared_lib, f"{self.name}_acados_update_params")(self.capsule, stage, value_data, value_.shape[0])==0
        else:
            if field_ not in constraints_fields + cost_fields + out_fields + mem_fields + sens_fields:
                raise Exception(f"AcadosOcpSolver.set(): '{field}' is not a valid argument.\n"
                    f" Possible values are {constraints_fields + cost_fields + out_fields + mem_fields + sens_fields + ['p']}.")

            dims = self.__acados_lib.ocp_nlp_dims_get_from_attr(self.nlp_config, \
                self.nlp_dims, self.nlp_out, stage_, field)

            if value_.shape[0] != dims:
                msg = f'AcadosOcpSolver.set(): mismatching dimension for field "{field_}" '
                msg += f'with dimension {dims} (you have {value_.shape[0]})'
                raise Exception(msg)

            value_data = cast(value_.ctypes.data, POINTER(c_double))
            value_data_p = cast((value_data), c_void_p)

            if field_ in constraints_fields:
                self.__acados_lib.ocp_nlp_constraints_model_set(self.nlp_config, \
                    self.nlp_dims, self.nlp_in, stage, field, value_data_p)
            elif field_ in cost_fields:
                self.__acados_lib.ocp_nlp_cost_model_set(self.nlp_config, \
                    self.nlp_dims, self.nlp_in, stage, field, value_data_p)
            elif field_ in out_fields:
                self.__acados_lib.ocp_nlp_out_set(self.nlp_config, \
                    self.nlp_dims, self.nlp_out, stage, field, value_data_p)
            elif field_ in mem_fields:
                self.__acados_lib.ocp_nlp_set(self.nlp_solver, stage, field, value_data_p)
            elif field_ in sens_fields:
                self.__acados_lib.ocp_nlp_out_set.argtypes = \
                    [c_void_p, c_void_p, c_void_p, c_int, c_char_p, c_void_p]
                self.__acados_lib.ocp_nlp_out_set(self.nlp_config, \
                    self.nlp_dims, self.sens_out, stage, field, value_data_p)
            # also set z_guess, when setting z.
            if field_ == 'z':
                field = 'z_guess'.encode('utf-8')
                self.__acados_lib.ocp_nlp_set(self.nlp_solver, stage, field, value_data_p)
        return


    def reset_sens_out(self):
        self.__acados_lib.ocp_nlp_out_set_values_to_zero(self.nlp_config, self.nlp_dims, self.sens_out)


    def cost_get(self, stage_: int, field_: str) -> np.ndarray:
        """
        Get numerical data in the cost module of the solver.

            :param stage: integer corresponding to shooting node
            :param field: string in ['yref', 'W', 'ext_cost_num_hess', 'zl', 'zu', 'Zl', 'Zu', 'scaling']
        """

        if not isinstance(stage_, int):
            raise Exception('stage should be integer.')
        elif stage_ < 0 or stage_ > self.N:
            raise Exception(f'stage should be in [0, N], got {stage_}')

        field = field_.encode('utf-8')
        stage = c_int(stage_)

        dims_ = np.zeros((2,), dtype=np.intc, order="C")
        dims_data = cast(dims_.ctypes.data, POINTER(c_int))

        self.__acados_lib.ocp_nlp_cost_dims_get_from_attr(self.nlp_config, \
            self.nlp_dims, self.nlp_out, stage_, field, dims_data)

        # vector-valued fields
        if field_ in ['yref', 'zl', 'zu', 'Zl', 'Zu', 'scaling']:
            dims = (dims_[0],)
        else:
            dims = tuple(dims_)

        out = np.zeros(dims, dtype=np.float64, order="F")
        out_data = cast(out.ctypes.data, POINTER(c_double))

        self.__acados_lib.ocp_nlp_cost_model_get(self.nlp_config, \
            self.nlp_dims, self.nlp_in, stage, field, out_data)

        return out


    def cost_set(self, stage_: int, field_: str, value_, api='warn'):
        """
        Set numerical data in the cost module of the solver.

            :param stage: integer corresponding to shooting node
            :param field: string, e.g. 'yref', 'W', 'ext_cost_num_hess', 'zl', 'zu', 'Zl', 'Zu', 'scaling'
            :param value: of appropriate size

        Note: by default the cost is scaled with the time step, and the terminal cost term scaled with 1.
        This can be overwritten by setting the 'scaling' field.
        """
        # cast value_ to avoid conversion issues
        if isinstance(value_, (float, int)):
            value_ = np.array([value_])

        if not isinstance(stage_, int):
            raise Exception('stage should be integer.')
        elif stage_ < 0 or stage_ > self.N:
            raise Exception(f'stage should be in [0, N], got {stage_}')

        value_ = value_.astype(float)
        field = field_.encode('utf-8')
        stage = c_int(stage_)

        dims = np.zeros((2,), dtype=np.intc, order="C")
        dims_data = cast(dims.ctypes.data, POINTER(c_int))

        self.__acados_lib.ocp_nlp_cost_dims_get_from_attr(self.nlp_config, \
            self.nlp_dims, self.nlp_out, stage_, field, dims_data)

        value_shape = value_.shape
        if len(value_shape) == 1:
            value_shape = (value_shape[0], 0)

        elif len(value_shape) == 2:
            if api=='old':
                pass
            elif api=='warn':
                if not np.all(np.ravel(value_, order='F')==np.ravel(value_, order='K')):
                    raise Exception("Ambiguity in API detected.\n"
                                    "Are you making an acados model from scratch? Add api='new' to cost_set and carry on.\n"
                                    "Are you seeing this error suddenly in previously running code? Read on.\n"
                                    f"  You are relying on a now-fixed bug in cost_set for field '{field_}'.\n" +
                                    "  acados_template now correctly passes on any matrices to acados in column major format.\n" +
                                    "  Two options to fix this error: \n" +
                                    "   * Add api='old' to cost_set to restore old incorrect behaviour\n" +
                                    "   * Add api='new' to cost_set and remove any unnatural manipulation of the value argument " +
                                    "such as non-mathematical transposes, reshaping, casting to fortran order, etc... " +
                                    "If there is no such manipulation, then you have probably been getting an incorrect solution before.")
                # Get elements in column major order
                value_ = np.ravel(value_, order='F')
            elif api=='new':
                # Get elements in column major order
                value_ = np.ravel(value_, order='F')
            else:
                raise Exception("Unknown api: '{}'".format(api))

        if value_shape != tuple(dims):
            raise Exception('AcadosOcpSolver.cost_set(): mismatching dimension' +
                f' for field "{field_}" at stage {stage} with dimension {tuple(dims)} (you have {value_shape})')

        value_data = cast(value_.ctypes.data, POINTER(c_double))
        value_data_p = cast((value_data), c_void_p)

        self.__acados_lib.ocp_nlp_cost_model_set(self.nlp_config, \
            self.nlp_dims, self.nlp_in, stage, field, value_data_p)

        return


    def constraints_get(self, stage_: int, field_: str) -> np.ndarray:
        """
        Get numerical data in the constraint module of the solver.

            :param stage: integer corresponding to shooting node
            :param field: string in ['lbx', 'ubx', 'lbu', 'ubu', 'lg', 'ug', 'lh', 'uh', 'uphi', 'C', 'D']
        """

        if not isinstance(stage_, int):
            raise Exception('stage should be integer.')
        elif stage_ < 0 or stage_ > self.N:
            raise Exception(f'stage should be in [0, N], got {stage_}')

        field = field_.encode('utf-8')
        stage = c_int(stage_)

        dims_ = np.zeros((2,), dtype=np.intc, order="C")
        dims_data = cast(dims_.ctypes.data, POINTER(c_int))

        self.__acados_lib.ocp_nlp_cost_dims_get_from_attr(self.nlp_config, \
            self.nlp_dims, self.nlp_out, stage_, field, dims_data)

        # check whether field is vector-valued
        if field_ not in ['C', 'D']:
            dims = (dims_[0],)
        else:
            dims = tuple(dims_)

        out = np.zeros(dims, dtype=np.float64, order="F")
        out_data = cast(out.ctypes.data, POINTER(c_double))

        self.__acados_lib.ocp_nlp_constraints_model_get(self.nlp_config, \
            self.nlp_dims, self.nlp_in, stage, field, out_data)

        return out


    def constraints_set(self, stage_: int, field_: str, value_: np.ndarray, api='warn'):
        """
        Set numerical data in the constraint module of the solver.

            :param stage: integer corresponding to shooting node
            :param field: string in ['lbx', 'ubx', 'lbu', 'ubu', 'lg', 'ug', 'lh', 'uh', 'uphi', 'C', 'D']
            :param value: of appropriate size
        """
        # cast value_ to avoid conversion issues
        if isinstance(value_, (float, int)):
            value_ = np.array([value_])
        value_ = value_.astype(float)

        if not isinstance(stage_, int):
            raise Exception('stage should be integer.')
        elif stage_ < 0 or stage_ > self.N:
            raise Exception(f'stage should be in [0, N], got {stage_}')

        field = field_.encode('utf-8')
        stage = c_int(stage_)

        dims = np.zeros((2,), dtype=np.intc, order="C")
        dims_data = cast(dims.ctypes.data, POINTER(c_int))

        self.__acados_lib.ocp_nlp_constraint_dims_get_from_attr(self.nlp_config, \
            self.nlp_dims, self.nlp_out, stage_, field, dims_data)

        value_shape = value_.shape
        if len(value_shape) == 1:
            value_shape = (value_shape[0], 0)
        elif len(value_shape) == 2:
            if api=='old':
                pass
            elif api=='warn':
                if not np.all(np.ravel(value_, order='F')==np.ravel(value_, order='K')):
                    raise Exception("Ambiguity in API detected.\n"
                                    "Are you making an acados model from scrach? Add api='new' to constraints_set and carry on.\n"
                                    "Are you seeing this error suddenly in previously running code? Read on.\n"
                                    f"  You are relying on a now-fixed bug in constraints_set for field '{field}'.\n" +
                                    "  acados_template now correctly passes on any matrices to acados in column major format.\n" +
                                    "  Two options to fix this error: \n" +
                                    "   * Add api='old' to constraints_set to restore old incorrect behaviour\n" +
                                    "   * Add api='new' to constraints_set and remove any unnatural manipulation of the value argument " +
                                    "such as non-mathematical transposes, reshaping, casting to fortran order, etc... " +
                                    "If there is no such manipulation, then you have probably been getting an incorrect solution before.")
                # Get elements in column major order
                value_ = np.ravel(value_, order='F')
            elif api=='new':
                # Get elements in column major order
                value_ = np.ravel(value_, order='F')
            else:
                raise Exception(f"Unknown api: '{api}'")

        if value_shape != tuple(dims):
            raise Exception(f'AcadosOcpSolver.constraints_set(): mismatching dimension' +
                f' for field "{field_}" at stage {stage} with dimension {tuple(dims)} (you have {value_shape})')

        value_data = cast(value_.ctypes.data, POINTER(c_double))
        value_data_p = cast((value_data), c_void_p)

        self.__acados_lib.ocp_nlp_constraints_model_set(self.nlp_config, \
            self.nlp_dims, self.nlp_in, stage, field, value_data_p)

        return


    def get_hessian_block(self, stage: int) -> np.ndarray:
        """
        Get Hessian block from last QP at stage i
        In HPIPM form [[R, S^T], [S, Q]]
        """
        Q_mat = self.get_from_qp_in(stage, 'Q')
        R_mat = self.get_from_qp_in(stage, 'R')
        S_mat = self.get_from_qp_in(stage, 'S')
        hess_block = scipy.linalg.block_diag(R_mat, Q_mat)
        nu = R_mat.shape[0]
        hess_block[nu:, :nu] = S_mat.T
        hess_block[:nu, nu:] = S_mat
        return hess_block


    def get_from_qp_in(self, stage_: int, field_: str):
        """
        Get numerical data from the current QP.

            :param stage: integer corresponding to shooting node
            :param field: string in ['A', 'B', 'b', 'Q', 'R', 'S', 'q', 'r', 'C', 'D', 'lg', 'ug', 'lbx', 'ubx', 'lbu', 'ubu']

        Note:
        - additional supported fields are ['P', 'K', 'Lr'], which can be extracted form QP solver PARTIAL_CONDENSING_HPIPM.
        - for PARTIAL_CONDENSING_* QP solvers, the following additional fields are available: ['pcond_Q', 'pcond_R', 'pcond_S']
        """
        # idx* should be added too..
        if not isinstance(stage_, int):
            raise TypeError("stage should be int")
        if stage_ > self.N:
            raise Exception("stage should be <= self.N")
        if field_ in self.__qp_dynamics_fields and stage_ >= self.N:
            raise ValueError(f"dynamics field {field_} not available at terminal stage")
        if field_ not in self.__qp_dynamics_fields + self.__qp_cost_fields + self.__qp_constraint_fields + self.__qp_pc_hpipm_fields + self.__qp_pc_fields + self.__qp_constraint_int_fields:
            raise Exception(f"field {field_} not supported.")
        if field_ in self.__qp_pc_hpipm_fields:
            if self.acados_ocp.solver_options.qp_solver != "PARTIAL_CONDENSING_HPIPM" or self.acados_ocp.solver_options.qp_solver_cond_N != self.acados_ocp.solver_options.N_horizon:
                raise Exception(f"field {field_} only works for PARTIAL_CONDENSING_HPIPM QP solver with qp_solver_cond_N == N.")
            if field_ in ["P", "K", "p"] and stage_ == 0 and self.acados_ocp.dims.nbxe_0 > 0:
                raise Exception(f"getting field {field_} at stage 0 only works without x0 elimination (see nbxe_0).")
        if field_ in self.__qp_pc_fields and not self.acados_ocp.solver_options.qp_solver.startswith("PARTIAL_CONDENSING"):
            raise Exception(f"field {field_} only works for PARTIAL_CONDENSING QP solvers.")

        field = field_.encode('utf-8')
        stage = c_int(stage_)

        # get dims
        dims = np.zeros((2,), dtype=np.intc, order="C")
        dims_data = cast(dims.ctypes.data, POINTER(c_int))

        self.__acados_lib.ocp_nlp_qp_dims_get_from_attr(self.nlp_config, \
            self.nlp_dims, self.nlp_out, stage_, field, dims_data)

        # create output data
        if field_ in self.__qp_constraint_int_fields:
            out =np.zeros((np.prod(dims),), dtype=np.int32, order="C")
        else:
            out = np.zeros((np.prod(dims),), dtype=np.float64, order="C")
        out = out.reshape(dims[0], dims[1], order='F')

        out_data = cast(out.ctypes.data, POINTER(c_double))
        out_data_p = cast((out_data), c_void_p)

        # call getter
        self.__acados_lib.ocp_nlp_get_at_stage(self.nlp_solver, stage, field, out_data_p)

        if field_ in ["Q", "R"]:
            # make symmetric: copy lower triangular part to upper triangular part
            out = np.tril(out) + np.tril(out, -1).T

        return out


    def __ocp_nlp_get_from_iterate(self, iteration_, stage_, field_):
        stage = c_int(stage_)
        field = field_.encode('utf-8')
        iteration = c_int(iteration_)
        dim = self.__acados_lib.ocp_nlp_dims_get_from_attr(self.nlp_config, self.nlp_dims, self.nlp_out,
                    stage, field)

        out = np.zeros((dim,), dtype=np.float64, order="C")
        out_data = cast(out.ctypes.data, POINTER(c_double))
        out_data_p = cast((out_data), c_void_p)
        self.__acados_lib.ocp_nlp_get_from_iterate(self.nlp_solver, iteration, stage, field, out_data_p)
        return out

    def get_iterate(self, iteration: int) -> AcadosOcpIterate:

        nlp_iter = self.get_stats('nlp_iter')
        if iteration < -1 or iteration > nlp_iter:
            raise Exception("get_iterate: iteration needs to be nonnegative and <= nlp_iter or -1.")

        if not self.acados_ocp.solver_options.store_iterates:
            raise Exception("get_iterate: the solver option store_iterates needs to be true in order to get iterates.")

        if self.acados_ocp.solver_options.nlp_solver_type == "SQP_RTI":
            raise Exception("get_iterate: SQP_RTI not supported.")

        # set to nlp_iter if -1
        iteration = nlp_iter if iteration == -1 else iteration

        x_traj = []
        u_traj = []
        z_traj = []
        sl_traj = []
        su_traj = []
        pi_traj = []
        lam_traj = []

        for n in range(self.acados_ocp.solver_options.N_horizon):
            x_traj.append(self.__ocp_nlp_get_from_iterate(iteration, n, "x"))
            u_traj.append(self.__ocp_nlp_get_from_iterate(iteration, n, "u"))
            z_traj.append(self.__ocp_nlp_get_from_iterate(iteration, n, "z"))
            sl_traj.append(self.__ocp_nlp_get_from_iterate(iteration, n, "sl"))
            su_traj.append(self.__ocp_nlp_get_from_iterate(iteration, n, "su"))
            pi_traj.append(self.__ocp_nlp_get_from_iterate(iteration, n, "pi"))
            lam_traj.append(self.__ocp_nlp_get_from_iterate(iteration, n, "lam"))

        n = self.acados_ocp.solver_options.N_horizon
        x_traj.append(self.__ocp_nlp_get_from_iterate(iteration, n, "x"))
        sl_traj.append(self.__ocp_nlp_get_from_iterate(iteration, n, "sl"))
        su_traj.append(self.__ocp_nlp_get_from_iterate(iteration, n, "su"))
        lam_traj.append(self.__ocp_nlp_get_from_iterate(iteration, n, "lam"))

        iterate = AcadosOcpIterate(x_traj=tuple(x_traj),
                                   u_traj=tuple(u_traj),
                                   z_traj=tuple(z_traj),
                                   sl_traj=tuple(sl_traj),
                                   su_traj=tuple(su_traj),
                                   pi_traj=tuple(pi_traj),
                                   lam_traj=tuple(lam_traj))

        return iterate


    def get_iterates(self) -> AcadosOcpIterates:
        return AcadosOcpIterates(iterate_list=[self.get_iterate(n) for n in range(self.get_stats('nlp_iter')+1)])


    def dims_get(self, field_, stage_):
        return self.__acados_lib.ocp_nlp_dims_get_from_attr(self.nlp_config, self.nlp_dims, self.nlp_out,
                    c_int(stage_), field_.encode('utf-8'))


    def options_set(self, field_, value_):
        """
        Set options of the solver.

            :param field: string, possible values are:
                'print_level', 'rti_phase', 'nlp_solver_max_iter, 'as_rti_level',
                'tol_eq', 'tol_stat', 'tol_ineq', 'tol_comp',
                'qp_tol_stat', 'qp_tol_eq', 'qp_tol_ineq', 'qp_tol_comp', 'qp_tau_min',
                'qp_warm_start', 'qp_mu0', 'qp_print_level', 'warm_start_first_qp',
                'globalization_fixed_step_length', 'globalization_alpha_min', 'globalization_alpha_reduction',
                'globalization_line_search_use_sufficient_descent', 'globalization_full_step_dual', 'globalization_use_SOC',
                'globalization_funnel_init_upper_bound', 'globalization_funnel_sufficient_decrease_factor',
                'globalization_funnel_kappa', 'globalization_funnel_fraction_switching_condition',
                'globalization_funnel_initial_penalty_parameter', 'globalization_funnel_init_increase_factor',
                'levenberg_marquardt',
                'adaptive_levenberg_marquardt_lam', 'adaptive_levenberg_marquardt_mu_min', 'adaptive_levenberg_marquardt_mu0',

            :param value: of type int, float, string, bool

            - qp_tol_stat: QP solver tolerance stationarity
            - qp_tol_eq: QP solver tolerance equalities
            - qp_tol_ineq: QP solver tolerance inequalities
            - qp_tol_comp: QP solver tolerance complementarity
            - qp_tau_min: for HPIPM QP solvers: minimum value of barrier parameter in HPIPM
            - qp_mu0: for HPIPM QP solvers: initial value for complementarity slackness
            - warm_start_first_qp: indicates if first QP in SQP is warm_started
            - rti_phase: 0: PREPARATION_AND_FEEDBACK, 1: PREPARATION, 2: FEEDBACK; only support for nlp_solver = 'SQP_RTI'
        """
        int_fields = ['print_level',
                      'rti_phase',
                      'globalization_line_search_use_sufficient_descent',
                      'globalization_full_step_dual',
                      'globalization_use_SOC',
                      'warm_start_first_qp',
                      'as_rti_level',
                      'max_iter',
                      'nlp_solver_max_iter',
                      'qp_warm_start',
                      'qp_print_level']
        double_fields = ['globalization_fixed_step_length',
                         'globalization_alpha_min',
                         'globalization_alpha_reduction',
                         'globalization_eps_sufficient_descent',
                         'globalization_funnel_init_increase_factor',
                         'globalization_funnel_init_upper_bound',
                         'globalization_funnel_sufficient_decrease_factor',
                         'globalization_funnel_kappa',
                         'globalization_funnel_fraction_switching_condition',
                         'globalization_funnel_initial_penalty_parameter',
                         'levenberg_marquardt',
                         'adaptive_levenberg_marquardt_lam',
                         'adaptive_levenberg_marquardt_mu_min',
                         'adaptive_levenberg_marquardt_mu0',
                         'tol_eq',
                         'tol_stat',
                         'tol_ineq',
                         'tol_comp',
                         'qp_tol_stat',
                         'qp_tol_eq',
                         'qp_tol_ineq',
                         'qp_tol_comp',
                         'qp_tau_min',
                         'qp_mu0']
        string_fields = []
        bool_fields = ['with_adaptive_levenberg_marquardt']

        # check field availability and type
        if field_ in int_fields:
            if not isinstance(value_, int):
                raise Exception(f'solver option \'{field_}\' must be of type int. You have {type(value_)}.')
            else:
                value_ctypes = c_int(value_)
        elif field_ in double_fields:
            if not isinstance(value_, float):
                raise Exception(f'solver option \'{field_}\' must be of type float. You have {type(value_)}.')
            else:
                value_ctypes = c_double(value_)
        elif field_ in bool_fields:
            if not isinstance(value_, bool):
                raise Exception(f'solver option \'{field_}\' must be of type bool. You have {type(value_)}.')
            else:
                value_ctypes = c_bool(value_)
        elif field_ in string_fields:
            if not isinstance(value_, str):
                raise Exception(f'solver option \'{field_}\' must be of type str. You have {type(value_)}.')
            else:
                value_ctypes = value_.encode('utf-8')
        else:
            fields = ', '.join(int_fields + double_fields + string_fields)
            raise Exception(f'AcadosOcpSolver.options_set() does not support field \'{field_}\'.\n'\
                f' Possible values are {fields}.')


        if (field_ == 'max_iter' or field_ == 'nlp_solver_max_iter') and value_ > self.__solver_options['nlp_solver_max_iter']:
            raise Exception('AcadosOcpSolver.options_set() cannot increase nlp_solver_max_iter' \
                    f' above initial value {self.__nlp_solver_max_iter} (you have {value_})')
            return

        if field_ == 'rti_phase':
            if value_ < 0 or value_ > 2:
                raise Exception('AcadosOcpSolver.options_set(): argument \'rti_phase\' can '
                    'take only values 0, 1, 2 for SQP-RTI-type solvers')
            if self.__solver_options['nlp_solver_type'] != 'SQP_RTI' and value_ > 0:
                raise Exception('AcadosOcpSolver.options_set(): argument \'rti_phase\' can '
                    'take only value 0 for SQP-type solvers')

        # encode
        field = field_.encode('utf-8')

        # call C interface
        if field_ in string_fields:
            self.__acados_lib.ocp_nlp_solver_opts_set(self.nlp_config, \
                self.nlp_opts, field, value_ctypes)
        else:
            self.__acados_lib.ocp_nlp_solver_opts_set(self.nlp_config, \
                self.nlp_opts, field, byref(value_ctypes))
        return


    def set_params_sparse(self, stage_: int, idx_values_: np.ndarray, param_values_):
        """
        set parameters of the solvers external function partially:
        Pseudo: solver.param[idx_values] = param_values;
        Parameters:

            :param stage: integer corresponding to shooting node
            :param idx_values: 0 based np array (or iterable) of integers: indices of parameter to be set
            :param param_values: new parameter values as numpy array
        """

        if not isinstance(stage_, int):
            raise Exception('stage should be integer.')
        elif stage_ < 0 or stage_ > self.N:
            raise Exception(f'stage should be in [0, N], got {stage_}')

        # if not isinstance(idx_values_, np.ndarray) or not issubclass(type(idx_values_[0]), np.integer):
        #     raise Exception('idx_values_ must be np.array of integers.')

        if not isinstance(param_values_, np.ndarray):
            raise Exception('param_values_ must be np.array.')
        elif np.float64 != param_values_.dtype:
            raise TypeError('param_values_ must be np.array of float64.')

        if param_values_.shape[0] != len(idx_values_):
            raise Exception(f'param_values_ and idx_values_ must be of the same size.' +
                 f' Got sizes idx {param_values_.shape[0]}, param_values {len(idx_values_)}.')

        p_dimension = self.__acados_lib.ocp_nlp_dims_get_from_attr(self.nlp_config, self.nlp_dims, self.nlp_out, stage_, "p".encode('utf-8'))
        if any(idx_values_ >= p_dimension):
            raise Exception(f'idx_values_ contains value >= np = {p_dimension} for stage {stage_}.')

        stage = c_int(stage_)
        n_update = c_int(len(param_values_))

        param_data = cast(param_values_.ctypes.data, POINTER(c_double))
        c_idx_values = np.ascontiguousarray(idx_values_, dtype=np.intc)
        idx_data = cast(c_idx_values.ctypes.data, POINTER(c_int))

        getattr(self.shared_lib, f"{self.name}_acados_update_params_sparse") \
                                    (self.capsule, stage, idx_data, param_data, n_update)

    def set_p_global_and_precompute_dependencies(self, data_: np.ndarray):
        """
        Sets values of p_global and precomputes all parts of the CasADi graphs of all other functions that only depend on p_global.
        """

        # checks
        np_global = self.__acados_lib.ocp_nlp_dims_get_from_attr(self.nlp_config, \
                self.nlp_dims, self.nlp_out, 0, "p_global".encode('utf-8'))
        if not isinstance(data_, np.ndarray):
            raise Exception('data must be np.array.')
        if np.float64 != data_.dtype:
            raise TypeError('data must be np.array of float64.')
        if data_.ndim != 1:
            raise Exception('data must be one-dimensional np array.')

        data = np.ascontiguousarray(data_, dtype=np.float64)
        c_data = cast(data.ctypes.data, POINTER(c_double))
        data_len = len(data)
        if data_len != np_global:
            raise Exception(f'data must have length {np_global}, got {data_len}.')

        status = getattr(self.shared_lib, f"{self.name}_acados_set_p_global_and_precompute_dependencies")(self.capsule, c_data, data_len)
        if self.__save_p_global:
            self.__p_global_values = data_

        return status


    def __del__(self):
        if self.solver_created:
            getattr(self.shared_lib, f"{self.name}_acados_free")(self.capsule)
            getattr(self.shared_lib, f"{self.name}_acados_free_capsule")(self.capsule)

            try:
                self.dlclose(self.shared_lib._handle)
            except:
                print(f"WARNING: acados Python interface could not close shared_lib handle of AcadosOcpSolver {self.name}.\n",
                     "Attempting to create a new one with the same name will likely result in the old one being used!")
                pass
