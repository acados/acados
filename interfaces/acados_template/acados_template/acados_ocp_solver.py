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
from copy import deepcopy
from ctypes import (POINTER, byref, c_char_p, c_double, c_int, c_int64,
                    c_void_p, cast)
if os.name == 'nt':
    from ctypes import wintypes
    from ctypes import WinDLL as DllLoader
else:
    from ctypes import CDLL as DllLoader
from datetime import datetime
from typing import Union, Optional

import numpy as np
import scipy.linalg

from .acados_ocp import AcadosOcp
from .acados_model import AcadosModel
from .acados_ocp_cost import AcadosOcpCost
from .acados_ocp_constraints import AcadosOcpConstraints
from .acados_dims import AcadosOcpDims
from .acados_ocp_options import AcadosOcpOptions

from .builders import CMakeBuilder

from .casadi_function_generation import ocp_generate_external_functions, mocp_generate_external_functions

from .gnsf.detect_gnsf_structure import detect_gnsf_structure
from .utils import (format_class_dict,
                    get_shared_lib_ext, get_shared_lib_prefix, get_shared_lib_dir, get_python_interface_path,
                    make_object_json_dumpable, render_template,
                    set_up_imported_gnsf_model, verbose_system_call)
from .zoro_description import ZoroDescription

from .acados_multiphase_ocp import AcadosMultiphaseOcp


def get_simulink_default_opts():
    python_interface_path = get_python_interface_path()
    abs_path = os.path.join(python_interface_path, 'simulink_default_opts.json')
    with open(abs_path , 'r') as f:
        simulink_default_opts = json.load(f)
    return simulink_default_opts


def ocp_formulation_json_dump(ocp: AcadosOcp, json_file: str, simulink_opts: Optional[dict]=None) -> None:
    ocp_nlp_dict = ocp.to_dict()
    # TODO: move simulink_opts to class AcadosOcp
    if simulink_opts is not None:
        ocp_nlp_dict['simulink_opts'] = simulink_opts

    with open(json_file, 'w') as f:
        json.dump(ocp_nlp_dict, f, default=make_object_json_dumpable, indent=4, sort_keys=True)
    return

def mocp_formulation_json_dump(ocp: AcadosMultiphaseOcp, json_file: str) -> None:
    ocp_nlp_dict = ocp.to_dict()
    with open(json_file, 'w') as f:
        json.dump(ocp_nlp_dict, f, default=make_object_json_dumpable, indent=4, sort_keys=True)
    return


def ocp_get_default_cmake_builder() -> CMakeBuilder:
    """
    If :py:class:`~acados_template.acados_ocp_solver.AcadosOcpSolver` is used with `CMake` this function returns a good first setting.
    :return: default :py:class:`~acados_template.builders.CMakeBuilder`
    """
    cmake_builder = CMakeBuilder()
    cmake_builder.options_on = ['BUILD_ACADOS_OCP_SOLVER_LIB']
    return cmake_builder



def ocp_render_templates(ocp: AcadosOcp, json_file: str, cmake_builder=None, simulink_opts=None):

    # check json file
    json_path = os.path.abspath(json_file)
    if not os.path.exists(json_path):
        raise Exception(f'Path "{json_path}" not found!')

    template_list = __ocp_get_template_list(ocp, cmake_builder=cmake_builder, simulink_opts=simulink_opts)

    # Render templates
    for tup in template_list:
        if len(tup) > 2:
            output_dir = tup[2]
        else:
            output_dir = ocp.code_export_directory
        render_template(tup[0], tup[1], output_dir, json_path)

    # Custom templates
    acados_template_path = os.path.dirname(os.path.abspath(__file__))
    custom_template_glob = os.path.join(acados_template_path, 'custom_update_templates', '*')
    for tup in ocp.solver_options.custom_templates:
        render_template(tup[0], tup[1], ocp.code_export_directory, json_path, template_glob=custom_template_glob)

    return



def __mocp_get_template_list(ocp: AcadosMultiphaseOcp, cmake_builder=None, simulink_opts=None) -> list:
    """
    returns a list of tuples in the form:
    (input_filename, output_filname)
    or
    (input_filename, output_filname, output_directory)
    """
    name = ocp.name
    template_list = []

    template_list.append(('main_multi.in.c', f'main_{name}.c'))
    template_list.append(('acados_multi_solver.in.h', f'acados_solver_{name}.h'))
    template_list.append(('acados_multi_solver.in.c', f'acados_solver_{name}.c'))
    # template_list.append(('acados_solver.in.pxd', f'acados_solver.pxd'))
    if cmake_builder is not None:
        raise NotImplementedError('CMake not yet supported for multiphase OCPs.')
        template_list.append(('CMakeLists.in.txt', 'CMakeLists.txt'))
    else:
        template_list.append(('multi_Makefile.in', 'Makefile'))

    # Simulink
    if simulink_opts is not None:
        raise NotImplementedError('Simulink not yet supported for multiphase OCPs.')

    return template_list


def mocp_render_templates(ocp: AcadosMultiphaseOcp, json_file: str, cmake_builder=None, simulink_opts=None):

    # model templates
    for i, dummy_ocp in enumerate(ocp.dummy_ocp_list):
        # this is the only option that can vary and influence external functions to be generated
        dummy_ocp.solver_options.integrator_type = ocp.mocp_opts.integrator_type[i]
        template_list = __ocp_get_external_function_header_templates(dummy_ocp)
        # dump dummy_ocp
        tmp_json_file = 'tmp_ocp.json'
        ocp_formulation_json_dump(dummy_ocp, json_file=tmp_json_file)
        tmp_json_path = os.path.abspath(tmp_json_file)

        # renter templates
        for tup in template_list:
            if len(tup) > 2:
                output_dir = tup[2]
            else:
                output_dir = ocp.code_export_directory
            render_template(tup[0], tup[1], output_dir, tmp_json_path)

    print("rendered model templates successfully")

    # check json file
    json_path = os.path.abspath(json_file)
    if not os.path.exists(json_path):
        raise Exception(f'Path "{json_path}" not found!')

    # solver templates
    template_list = __mocp_get_template_list(ocp, cmake_builder=cmake_builder, simulink_opts=simulink_opts)

    # Render templates
    for tup in template_list:
        if len(tup) > 2:
            output_dir = tup[2]
        else:
            output_dir = ocp.code_export_directory
        render_template(tup[0], tup[1], output_dir, json_path)

    # # Custom templates
    # acados_template_path = os.path.dirname(os.path.abspath(__file__))
    # custom_template_glob = os.path.join(acados_template_path, 'custom_update_templates', '*')
    # for tup in ocp.solver_options.custom_templates:
    #     render_template(tup[0], tup[1], ocp.code_export_directory, json_path, template_glob=custom_template_glob)
    print("\nmocp_render_templates: rendered solver templates successfully!\n")


    return


def __ocp_get_template_list(ocp: AcadosOcp, cmake_builder=None, simulink_opts=None) -> list:
    """
    returns a list of tuples in the form:
    (input_filename, output_filname)
    or
    (input_filename, output_filname, output_directory)
    """
    name = ocp.model.name
    template_list = []

    template_list.append(('main.in.c', f'main_{name}.c'))
    template_list.append(('acados_solver.in.c', f'acados_solver_{name}.c'))
    template_list.append(('acados_solver.in.h', f'acados_solver_{name}.h'))
    template_list.append(('acados_solver.in.pxd', f'acados_solver.pxd'))
    if cmake_builder is not None:
        template_list.append(('CMakeLists.in.txt', 'CMakeLists.txt'))
    else:
        template_list.append(('Makefile.in', 'Makefile'))

    # sim
    template_list.append(('acados_sim_solver.in.c', f'acados_sim_solver_{name}.c'))
    template_list.append(('acados_sim_solver.in.h', f'acados_sim_solver_{name}.h'))
    template_list.append(('main_sim.in.c', f'main_sim_{name}.c'))

    # model
    template_list += __ocp_get_external_function_header_templates(ocp)

    # Simulink
    if simulink_opts is not None:
        template_file = os.path.join('matlab_templates', 'acados_solver_sfun.in.c')
        template_list.append((template_file, f'acados_solver_sfunction_{name}.c'))
        template_file = os.path.join('matlab_templates', 'make_sfun.in.m')
        template_list.append((template_file, f'make_sfun_{name}.m'))
        template_file = os.path.join('matlab_templates', 'acados_sim_solver_sfun.in.c')
        template_list.append((template_file, f'acados_sim_solver_sfunction_{name}.c'))
        template_file = os.path.join('matlab_templates', 'make_sfun_sim.in.m')
        template_list.append((template_file, f'make_sfun_sim_{name}.m'))

    return template_list


def __ocp_get_external_function_header_templates(ocp: AcadosOcp) -> list:
    dims = ocp.dims
    cost = ocp.cost
    code_export_directory = ocp.code_export_directory
    name = ocp.model.name
    template_list = []

    # dynamics
    model_dir = os.path.join(code_export_directory, f'{name}_model')
    template_list.append(('model.in.h', f'{name}_model.h', model_dir))
    # constraints
    if any(np.array([dims.nh, dims.nh_e, dims.nh_0, dims.nphi, dims.nphi_e, dims.nphi_0]) > 0):
        constraints_dir = os.path.join(code_export_directory, f'{name}_constraints')
        template_list.append(('constraints.in.h', f'{name}_constraints.h', constraints_dir))
    # cost
    if any([cost.cost_type != 'LINEAR_LS', cost.cost_type_0 != 'LINEAR_LS', cost.cost_type_e != 'LINEAR_LS']):
        cost_dir = os.path.join(code_export_directory, f'{name}_cost')
        template_list.append(('cost.in.h', f'{name}_cost.h', cost_dir))

    return template_list




class AcadosOcpSolver:
    """
    Class to interact with the acados ocp solver C object.

        :param acados_ocp: type :py:class:`~acados_template.acados_ocp.AcadosOcp` or :py:class:`~acados_template.acados_multiphase_ocp.AcadosMultiphaseOcp` - description of the OCP for acados
        :param json_file: name for the json file used to render the templated code - default: acados_ocp_nlp.json
        :param simulink_opts: Options to configure Simulink S-function blocks, mainly to activate possible Inputs and Outputs
    """
    if os.name == 'nt':
        dlclose = DllLoader('kernel32', use_last_error=True).FreeLibrary
        dlclose.argtypes = [wintypes.HMODULE]
        winmode = 8 # why 8? what does that mean?
    else:
        dlclose = DllLoader(None).dlclose
        dlclose.argtypes = [c_void_p]
        winmode = None

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

        # set integrator time automatically
        acados_ocp.solver_options.Tsim = acados_ocp.solver_options.time_steps[0]

        # generate code (external functions and templated code)
        if isinstance(acados_ocp, AcadosOcp):
            ocp_generate_external_functions(acados_ocp)
            ocp_formulation_json_dump(acados_ocp, json_file, simulink_opts=simulink_opts)
            ocp_render_templates(acados_ocp, json_file, cmake_builder=cmake_builder, simulink_opts=simulink_opts)
        elif isinstance(acados_ocp, AcadosMultiphaseOcp):
            mocp_generate_external_functions(acados_ocp)
            mocp_formulation_json_dump(acados_ocp, json_file)
            mocp_render_templates(acados_ocp, json_file, cmake_builder=cmake_builder, simulink_opts=simulink_opts)

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
        cwd = os.getcwd()
        os.chdir(code_export_dir)

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
        os.chdir(cwd)


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


    def __init__(self, acados_ocp: Union[AcadosOcp, AcadosMultiphaseOcp], json_file=None, simulink_opts=None, build=True, generate=True, cmake_builder: CMakeBuilder = None, verbose=True):

        self.solver_created = False

        if isinstance(acados_ocp, AcadosOcp):
            if json_file is None:
                json_file = 'acados_ocp_nlp.json'
        elif isinstance(acados_ocp, AcadosMultiphaseOcp):
            if json_file is None:
                json_file = 'mocp.json'
        else:
            raise Exception('acados_ocp should be of type AcadosOcp or AcadosMultiphaseOcp.')

        if generate:
            self.generate(acados_ocp, json_file=json_file, simulink_opts=simulink_opts, cmake_builder=cmake_builder)
        else:
            acados_ocp.make_consistent()

        # load json, store options in object
        with open(json_file, 'r') as f:
            acados_ocp_json = json.load(f)
        if isinstance(acados_ocp, AcadosOcp):
            self.N = acados_ocp_json['dims']['N']
        elif isinstance(acados_ocp, AcadosMultiphaseOcp):
            self.N = acados_ocp_json['N_horizon']
        self.solver_options = acados_ocp_json['solver_options']
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
        if self.winmode is not None:
            self.__acados_lib = DllLoader(libacados_filepath, winmode=self.winmode)
        else:
            # for compatibility with older python versions
            self.__acados_lib = DllLoader(libacados_filepath)

        # find out if acados was compiled with OpenMP
        try:
            self.__acados_lib_uses_omp = getattr(self.__acados_lib, 'omp_get_thread_num') is not None
        except AttributeError as e:
            self.__acados_lib_uses_omp = False
        if verbose:
            if self.__acados_lib_uses_omp:
                print('acados was compiled with OpenMP.')
            else:
                print('acados was compiled without OpenMP.')
        libacados_ocp_solver_name = f'{lib_prefix}acados_ocp_solver_{self.name}{lib_ext}'
        self.shared_lib_name = os.path.join(code_export_directory, libacados_ocp_solver_name)

        # get shared_lib
        if self.winmode is not None:
            self.shared_lib = DllLoader(self.shared_lib_name, winmode=self.winmode)
        else:
            # for compatibility with older python versions
            self.shared_lib = DllLoader(self.shared_lib_name)

        # create capsule
        getattr(self.shared_lib, f"{self.name}_acados_create_capsule").restype = c_void_p
        self.capsule = getattr(self.shared_lib, f"{self.name}_acados_create_capsule")()

        # create solver
        getattr(self.shared_lib, f"{self.name}_acados_create").argtypes = [c_void_p]
        getattr(self.shared_lib, f"{self.name}_acados_create").restype = c_int
        assert getattr(self.shared_lib, f"{self.name}_acados_create")(self.capsule)==0
        self.solver_created = True

        self.acados_ocp = acados_ocp

        # get pointers solver
        self.__get_pointers_solver()

        self.status = 0

        # gettable fields
        self.__qp_dynamics_fields = ['A', 'B', 'b']
        self.__qp_cost_fields = ['Q', 'R', 'S', 'q', 'r']
        self.__qp_constraint_fields = ['C', 'D', 'lg', 'ug', 'lbx', 'ubx', 'lbu', 'ubu']
        self.__qp_pc_hpipm_fields = ['P', 'K', 'Lr', 'p']

        return

    def __get_pointers_solver(self):
        # """
        # Private function to get the pointers for solver
        # """
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


    def solve(self):
        """
        Solve the ocp with current input.
        """
        getattr(self.shared_lib, f"{self.name}_acados_solve").argtypes = [c_void_p]
        getattr(self.shared_lib, f"{self.name}_acados_solve").restype = c_int
        self.status = getattr(self.shared_lib, f"{self.name}_acados_solve")(self.capsule)

        return self.status


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

        getattr(self.shared_lib, f"{self.name}_acados_custom_update").argtypes = [c_void_p, POINTER(c_double), c_int]
        getattr(self.shared_lib, f"{self.name}_acados_custom_update").restype = c_int
        status = getattr(self.shared_lib, f"{self.name}_acados_custom_update")(self.capsule, c_data, data_len)

        return status


    def reset(self, reset_qp_solver_mem=1):
        """
        Sets current iterate to all zeros.
        """
        getattr(self.shared_lib, f"{self.name}_acados_reset").argtypes = [c_void_p, c_int]
        getattr(self.shared_lib, f"{self.name}_acados_reset").restype = c_int
        getattr(self.shared_lib, f"{self.name}_acados_reset")(self.capsule, reset_qp_solver_mem)

        return


    def set_new_time_steps(self, new_time_steps):
        """
        Set new time steps.
        Recreates the solver if N changes.

            :param new_time_steps: 1 dimensional np array of new time steps for the solver

            .. note:: This allows for different use-cases: either set a new size of time_steps or a new distribution of
                      the shooting nodes without changing the number, e.g., to reach a different final time. Both cases
                      do not require a new code export and compilation.
        """

        # unlikely but still possible
        if not self.solver_created:
            raise Exception('Solver was not yet created!')

        # check if time steps really changed in value
        if np.array_equal(self.solver_options['time_steps'], new_time_steps):
            return

        N = new_time_steps.size
        new_time_steps_data = cast(new_time_steps.ctypes.data, POINTER(c_double))

        # check if recreation of acados is necessary (no need to recreate acados if sizes are identical)
        if len(self.solver_options['time_steps']) == N:
            getattr(self.shared_lib, f"{self.name}_acados_update_time_steps").argtypes = [c_void_p, c_int, c_void_p]
            getattr(self.shared_lib, f"{self.name}_acados_update_time_steps").restype = c_int
            assert getattr(self.shared_lib, f"{self.name}_acados_update_time_steps")(self.capsule, N, new_time_steps_data) == 0
        else:  # recreate the solver with the new time steps
            self.solver_created = False

            # delete old memory (analog to __del__)
            getattr(self.shared_lib, f"{self.name}_acados_free").argtypes = [c_void_p]
            getattr(self.shared_lib, f"{self.name}_acados_free").restype = c_int
            getattr(self.shared_lib, f"{self.name}_acados_free")(self.capsule)

            # create solver with new time steps
            getattr(self.shared_lib, f"{self.name}_acados_create_with_discretization").argtypes = [c_void_p, c_int, c_void_p]
            getattr(self.shared_lib, f"{self.name}_acados_create_with_discretization").restype = c_int
            assert getattr(self.shared_lib, f"{self.name}_acados_create_with_discretization")(self.capsule, N, new_time_steps_data) == 0

            self.solver_created = True

            # get pointers solver
            self.__get_pointers_solver()

        # store time_steps, N
        self.solver_options['time_steps'] = new_time_steps
        self.N = N
        self.solver_options['Tsim'] = self.solver_options['time_steps'][0]


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
        if self.solver_options['qp_solver_cond_N'] != qp_solver_cond_N:
            self.solver_created = False

            # recreate the solver
            fun_name = f'{self.name}_acados_update_qp_solver_cond_N'
            getattr(self.shared_lib, fun_name).argtypes = [c_void_p, c_int]
            getattr(self.shared_lib, fun_name).restype = c_int
            assert getattr(self.shared_lib, fun_name)(self.capsule, qp_solver_cond_N) == 0

            # store the new value
            self.solver_options['qp_solver_cond_N'] = qp_solver_cond_N
            self.solver_created = True

            # get pointers solver
            self.__get_pointers_solver()


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

        if not (self.acados_ocp.solver_options.qp_solver == 'FULL_CONDENSING_HPIPM' or
                self.acados_ocp.solver_options.qp_solver == 'PARTIAL_CONDENSING_HPIPM'):
            raise Exception("Parametric sensitivities are only available with HPIPM as QP solver.")


        if not (
           (self.acados_ocp.solver_options.hessian_approx == 'EXACT' or
           (self.acados_ocp.cost.cost_type == 'LINEAR_LS' and
            self.acados_ocp.cost.cost_type_0 == 'LINEAR_LS' and
            self.acados_ocp.cost.cost_type_e == 'LINEAR_LS'))
            and
            self.acados_ocp.solver_options.regularize_method == 'NO_REGULARIZE' and
            self.acados_ocp.solver_options.levenberg_marquardt == 0
        ):
            raise Exception("Parametric sensitivities are only correct if an exact Hessian is used!")

        field_ = field
        field = field_.encode('utf-8')

        # checks
        if not isinstance(index, int):
            raise Exception('AcadosOcpSolver.eval_param_sens(): index must be Integer.')

        self.__acados_lib.ocp_nlp_dims_get_from_attr.argtypes = [c_void_p, c_void_p, c_void_p, c_int, c_char_p]
        self.__acados_lib.ocp_nlp_dims_get_from_attr.restype = c_int
        nx = self.__acados_lib.ocp_nlp_dims_get_from_attr(self.nlp_config, self.nlp_dims, self.nlp_out, 0, "x".encode('utf-8'))

        if index < 0 or index > nx:
            raise Exception(f'AcadosOcpSolver.eval_param_sens(): index must be in [0, nx-1], got: {index}.')

        # actual eval_param
        self.__acados_lib.ocp_nlp_eval_param_sens.argtypes = [c_void_p, c_char_p, c_int, c_int, c_void_p]
        self.__acados_lib.ocp_nlp_eval_param_sens.restype = None
        self.__acados_lib.ocp_nlp_eval_param_sens(self.nlp_solver, field, stage, index, self.sens_out)

        return


    def get(self, stage_: int, field_: str):
        """
        Get the last solution of the solver:

            :param stage: integer corresponding to shooting node
            :param field: string in ['x', 'u', 'z', 'pi', 'lam', 't', 'sl', 'su',]

            .. note:: regarding lam, t: \n
                    the inequalities are internally organized in the following order: \n
                    [ lbu lbx lg lh lphi ubu ubx ug uh uphi; \n
                      lsbu lsbx lsg lsh lsphi usbu usbx usg ush usphi]

            .. note:: pi: multipliers for dynamics equality constraints \n
                      lam: multipliers for inequalities \n
                      t: slack variables corresponding to evaluation of all inequalities (at the solution) \n
                      sl: slack variables of soft lower inequality constraints \n
                      su: slack variables of soft upper inequality constraints \n
        """

        out_fields = ['x', 'u', 'z', 'pi', 'lam', 't', 'sl', 'su']
        # mem_fields = ['sl', 'su']
        sens_fields = ['sens_u', "sens_x"]
        all_fields = out_fields + sens_fields

        field = field_

        if (field_ not in all_fields):
            raise Exception(f'AcadosOcpSolver.get(stage={stage_}, field={field_}): \'{field_}\' is an invalid argument.\
                    \n Possible values are {all_fields}.')

        if not isinstance(stage_, int):
            raise Exception(f'AcadosOcpSolver.get(stage={stage_}, field={field_}): stage index must be an integer, got type {type(stage_)}.')

        if stage_ < 0 or stage_ > self.N:
            raise Exception(f'AcadosOcpSolver.get(stage={stage_}, field={field_}): stage index must be in [0, {self.N}], got: {stage_}.')

        if stage_ == self.N and field_ == 'pi':
            raise Exception(f'AcadosOcpSolver.get(stage={stage_}, field={field_}): field \'{field_}\' does not exist at final stage {stage_}.')

        if field_ in sens_fields:
            field = field_.replace('sens_', '')

        field = field.encode('utf-8')

        self.__acados_lib.ocp_nlp_dims_get_from_attr.argtypes = \
            [c_void_p, c_void_p, c_void_p, c_int, c_char_p]
        self.__acados_lib.ocp_nlp_dims_get_from_attr.restype = c_int

        dims = self.__acados_lib.ocp_nlp_dims_get_from_attr(self.nlp_config, \
            self.nlp_dims, self.nlp_out, stage_, field)

        out = np.ascontiguousarray(np.zeros((dims,)), dtype=np.float64)
        out_data = cast(out.ctypes.data, POINTER(c_double))

        if (field_ in out_fields):
            self.__acados_lib.ocp_nlp_out_get.argtypes = \
                [c_void_p, c_void_p, c_void_p, c_int, c_char_p, c_void_p]
            self.__acados_lib.ocp_nlp_out_get(self.nlp_config, \
                self.nlp_dims, self.nlp_out, stage_, field, out_data)
        # elif field_ in mem_fields:
        #     self.__acados_lib.ocp_nlp_get_at_stage.argtypes = \
        #         [c_void_p, c_void_p, c_void_p, c_int, c_char_p, c_void_p]
        #     self.__acados_lib.ocp_nlp_get_at_stage(self.nlp_config, \
        #         self.nlp_dims, self.nlp_solver, stage_, field, out_data)
        elif field_ in sens_fields:
            self.__acados_lib.ocp_nlp_out_get.argtypes = \
                [c_void_p, c_void_p, c_void_p, c_int, c_char_p, c_void_p]
            self.__acados_lib.ocp_nlp_out_get(self.nlp_config, \
                self.nlp_dims, self.sens_out, stage_, field, out_data)

        return out


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

        if self.solver_options['nlp_solver_type'] == 'SQP':
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
        elif self.solver_options['nlp_solver_type'] == 'SQP_RTI':
            print('\niter\tqp_stat\tqp_iter')
            if stat.shape[0]>3:
                print('\tqp_res_stat\tqp_res_eq\tqp_res_ineq\tqp_res_comp')
            for jj in range(stat.shape[1]):
                print('{:d}\t{:d}\t{:d}'.format( int(stat[0][jj]), int(stat[1][jj]), int(stat[2][jj])))
                if stat.shape[0]>3:
                    print('\t{:e}\t{:e}\t{:e}\t{:e}'.format( \
                         stat[3][jj], stat[4][jj], stat[5][jj], stat[6][jj]))
            print('\n')

        return


    def store_iterate(self, filename: str = '', overwrite=False, verbose=True):
        """
        Stores the current iterate of the ocp solver in a json file.
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
                filename += datetime.utcnow().strftime('%Y-%m-%d-%H:%M:%S.%f') + '.json'

        # get iterate:
        solution = dict()

        lN = len(str(self.N+1))
        for i in range(self.N+1):
            i_string = f'{i:0{lN}d}'
            solution['x_'+i_string] = self.get(i,'x')
            solution['u_'+i_string] = self.get(i,'u')
            solution['z_'+i_string] = self.get(i,'z')
            solution['lam_'+i_string] = self.get(i,'lam')
            solution['t_'+i_string] = self.get(i, 't')
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
                filename += datetime.utcnow().strftime('%Y-%m-%d-%H:%M:%S.%f') + '.json'

        # get QP data:
        qp_data = dict()

        lN = len(str(self.N+1))
        for field in self.__qp_dynamics_fields:
            for i in range(self.N):
                qp_data[f'{field}_{i:0{lN}d}'] = self.get_from_qp_in(i,field)

        for field in self.__qp_constraint_fields + self.__qp_cost_fields:
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



    def load_iterate(self, filename:str, verbose: bool=True):
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


    def get_stats(self, field_: str) -> Union[int, float, np.ndarray]:
        """
        Get the information of the last solver call.

            :param field: string in ['statistics', 'time_tot', 'time_lin', 'time_sim', 'time_sim_ad', 'time_sim_la', 'time_qp', 'time_qp_solver_call', 'time_reg', 'sqp_iter', 'residuals', 'qp_iter', 'alpha']

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
            - time_reg: CPU time regularization
            - sqp_iter: number of SQP iterations
            - qp_stat: status of QP solver
            - qp_iter: vector of QP iterations for last SQP call
            - statistics: table with info about last iteration
            - stat_m: number of rows in statistics matrix
            - stat_n: number of columns in statistics matrix
            - residuals: residuals of last iterate
            - alpha: step sizes of SQP iterations
        """

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
                  'time_reg'
        ]
        fields = double_fields + [
                  'sqp_iter',
                  'qp_stat',
                  'qp_iter',
                  'statistics',
                  'stat_m',
                  'stat_n',
                  'residuals',
                  'alpha',
                ]
        field = field_.encode('utf-8')


        if field_ in ['sqp_iter', 'stat_m', 'stat_n']:
            out = c_int(0)
            self.__acados_lib.ocp_nlp_get.argtypes = [c_void_p, c_void_p, c_char_p, c_void_p]
            self.__acados_lib.ocp_nlp_get(self.nlp_config, self.nlp_solver, field, byref(out))
            return out.value

        elif field_ in double_fields:
            out = c_double(0)
            self.__acados_lib.ocp_nlp_get.argtypes = [c_void_p, c_void_p, c_char_p, c_void_p]
            self.__acados_lib.ocp_nlp_get(self.nlp_config, self.nlp_solver, field, byref(out))
            return out.value

        elif field_ == 'statistics':
            sqp_iter = self.get_stats("sqp_iter")
            stat_m = self.get_stats("stat_m")
            stat_n = self.get_stats("stat_n")
            min_size = min([stat_m, sqp_iter+1])
            out = np.ascontiguousarray(
                        np.zeros((stat_n+1, min_size)), dtype=np.float64)
            out_data = cast(out.ctypes.data, POINTER(c_double))
            self.__acados_lib.ocp_nlp_get.argtypes = [c_void_p, c_void_p, c_char_p, c_void_p]
            self.__acados_lib.ocp_nlp_get(self.nlp_config, self.nlp_solver, field, out_data)
            return out

        elif field_ == 'qp_stat':
            full_stats = self.get_stats('statistics')
            if self.solver_options['nlp_solver_type'] == 'SQP':
                return full_stats[5, :]
            elif self.solver_options['nlp_solver_type'] == 'SQP_RTI':
                return full_stats[1, :]

        elif field_ == 'qp_iter':
            full_stats = self.get_stats('statistics')
            if self.solver_options['nlp_solver_type'] == 'SQP':
                return full_stats[6, :]
            elif self.solver_options['nlp_solver_type'] == 'SQP_RTI':
                return full_stats[2, :]

        elif field_ == 'alpha':
            full_stats = self.get_stats('statistics')
            if self.solver_options['nlp_solver_type'] == 'SQP':
                return full_stats[7, :]
            else: # self.solver_options['nlp_solver_type'] == 'SQP_RTI':
                raise Exception("alpha values are not available for SQP_RTI")

        elif field_ == 'residuals':
            return self.get_residuals()

        elif field_ == 'residuals':
            if self.acados_ocp.solver_options.nlp_solver_type == 'SQP':
                full_stats = self.get_stats('statistics')
                if self.status != 2:
                    out = (full_stats.T)[-1][1:5]
                else: # when exiting with max_iter, residuals are computed for second last iterate only
                    out = (full_stats.T)[-2][1:5]
            else:
                Exception("residuals are not computed for SQP_RTI")


        else:
            raise Exception(f'AcadosOcpSolver.get_stats(): \'{field}\' is not a valid argument.'
                    + f'\n Possible values are {fields}.')


    def get_cost(self):
        """
        Returns the cost value of the current solution.
        """
        # compute cost internally
        self.__acados_lib.ocp_nlp_eval_cost.argtypes = [c_void_p, c_void_p, c_void_p]
        self.__acados_lib.ocp_nlp_eval_cost(self.nlp_solver, self.nlp_in, self.nlp_out)

        # create output array
        out = np.ascontiguousarray(np.zeros((1,)), dtype=np.float64)
        out_data = cast(out.ctypes.data, POINTER(c_double))

        # call getter
        self.__acados_lib.ocp_nlp_get.argtypes = [c_void_p, c_void_p, c_char_p, c_void_p]

        field = "cost_value".encode('utf-8')
        self.__acados_lib.ocp_nlp_get(self.nlp_config, self.nlp_solver, field, out_data)

        return out[0]


    def get_optimal_value_gradient(self):
        """
        Returns the gradient of the optimal value function w.r.t. the current initial state.
        Disclaimer: This function only returns reasonable values if the solver has converged for the current problem instance.
        """
        if not self.acados_ocp.constraints.has_x0:
            raise Exception("OCP does not have an initial state constraint.")

        nx = self.acados_ocp.dims.nx
        nbu = self.acados_ocp.dims.nbu
        lam = self.get(0, 'lam')
        nlam_non_slack = lam.shape[0]//2 - self.acados_ocp.dims.ns_0
        grad = lam[nbu:nbu+nx] - lam[nlam_non_slack+nbu : nlam_non_slack+nbu+nx]

        return grad


    def get_residuals(self, recompute=False):
        """
        Returns an array of the form [res_stat, res_eq, res_ineq, res_comp].
        This residual has to be computed for SQP_RTI solver, since it is not available by default.

        - res_stat: stationarity residual
        - res_eq: residual wrt equality constraints (dynamics)
        - res_ineq: residual wrt inequality constraints (constraints)
        - res_comp: residual wrt complementarity conditions
        """
        # compute residuals if RTI
        if self.solver_options['nlp_solver_type'] == 'SQP_RTI' or recompute:
            self.__acados_lib.ocp_nlp_eval_residuals.argtypes = [c_void_p, c_void_p, c_void_p]
            self.__acados_lib.ocp_nlp_eval_residuals(self.nlp_solver, self.nlp_in, self.nlp_out)

        # create output array
        out = np.ascontiguousarray(np.zeros((4, 1)), dtype=np.float64)
        out_data = cast(out.ctypes.data, POINTER(c_double))

        # call getters
        self.__acados_lib.ocp_nlp_get.argtypes = [c_void_p, c_void_p, c_char_p, c_void_p]

        field = "res_stat".encode('utf-8')
        self.__acados_lib.ocp_nlp_get(self.nlp_config, self.nlp_solver, field, out_data)

        out_data = cast(out[1].ctypes.data, POINTER(c_double))
        field = "res_eq".encode('utf-8')
        self.__acados_lib.ocp_nlp_get(self.nlp_config, self.nlp_solver, field, out_data)

        out_data = cast(out[2].ctypes.data, POINTER(c_double))
        field = "res_ineq".encode('utf-8')
        self.__acados_lib.ocp_nlp_get(self.nlp_config, self.nlp_solver, field, out_data)

        out_data = cast(out[3].ctypes.data, POINTER(c_double))
        field = "res_comp".encode('utf-8')
        self.__acados_lib.ocp_nlp_get(self.nlp_config, self.nlp_solver, field, out_data)
        return out.flatten()


    # Note: this function should not be used anymore, better use cost_set, constraints_set
    def set(self, stage_: int, field_: str, value_):
        """
        Set numerical data inside the solver.

            :param stage: integer corresponding to shooting node
            :param field: string in ['x', 'u', 'pi', 'lam', 't', 'p', 'xdot_guess', 'z_guess']

            .. note:: regarding lam, t: \n
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
        out_fields = ['x', 'u', 'pi', 'lam', 't', 'z', 'sl', 'su']
        mem_fields = ['xdot_guess', 'z_guess']

        if not isinstance(stage_, int):
            raise Exception('stage should be integer.')
        elif stage_ < 0 or stage_ > self.N:
            raise Exception(f'stage should be in [0, N], got {stage_}')

        # cast value_ to avoid conversion issues
        if isinstance(value_, (float, int)):
            value_ = np.array([value_])
        value_ = value_.astype(float)

        field = field_.encode('utf-8')

        stage = c_int(stage_)

        # treat parameters separately
        if field_ == 'p':
            getattr(self.shared_lib, f"{self.name}_acados_update_params").argtypes = [c_void_p, c_int, POINTER(c_double), c_int]
            getattr(self.shared_lib, f"{self.name}_acados_update_params").restype = c_int

            value_data = cast(value_.ctypes.data, POINTER(c_double))

            assert getattr(self.shared_lib, f"{self.name}_acados_update_params")(self.capsule, stage, value_data, value_.shape[0])==0
        else:
            if field_ not in constraints_fields + cost_fields + out_fields + mem_fields:
                raise Exception(f"AcadosOcpSolver.set(): '{field}' is not a valid argument.\n"
                    f" Possible values are {constraints_fields + cost_fields + out_fields + mem_fields + ['p']}.")

            self.__acados_lib.ocp_nlp_dims_get_from_attr.argtypes = \
                [c_void_p, c_void_p, c_void_p, c_int, c_char_p]
            self.__acados_lib.ocp_nlp_dims_get_from_attr.restype = c_int

            dims = self.__acados_lib.ocp_nlp_dims_get_from_attr(self.nlp_config, \
                self.nlp_dims, self.nlp_out, stage_, field)

            if value_.shape[0] != dims:
                msg = f'AcadosOcpSolver.set(): mismatching dimension for field "{field_}" '
                msg += f'with dimension {dims} (you have {value_.shape[0]})'
                raise Exception(msg)

            value_data = cast(value_.ctypes.data, POINTER(c_double))
            value_data_p = cast((value_data), c_void_p)

            if field_ in constraints_fields:
                self.__acados_lib.ocp_nlp_constraints_model_set.argtypes = \
                    [c_void_p, c_void_p, c_void_p, c_int, c_char_p, c_void_p]
                self.__acados_lib.ocp_nlp_constraints_model_set(self.nlp_config, \
                    self.nlp_dims, self.nlp_in, stage, field, value_data_p)
            elif field_ in cost_fields:
                self.__acados_lib.ocp_nlp_cost_model_set.argtypes = \
                    [c_void_p, c_void_p, c_void_p, c_int, c_char_p, c_void_p]
                self.__acados_lib.ocp_nlp_cost_model_set(self.nlp_config, \
                    self.nlp_dims, self.nlp_in, stage, field, value_data_p)
            elif field_ in out_fields:
                self.__acados_lib.ocp_nlp_out_set.argtypes = \
                    [c_void_p, c_void_p, c_void_p, c_int, c_char_p, c_void_p]
                self.__acados_lib.ocp_nlp_out_set(self.nlp_config, \
                    self.nlp_dims, self.nlp_out, stage, field, value_data_p)
            elif field_ in mem_fields:
                self.__acados_lib.ocp_nlp_set.argtypes = \
                    [c_void_p, c_void_p, c_int, c_char_p, c_void_p]
                self.__acados_lib.ocp_nlp_set(self.nlp_config, \
                    self.nlp_solver, stage, field, value_data_p)
            # also set z_guess, when setting z.
            if field_ == 'z':
                field = 'z_guess'.encode('utf-8')
                self.__acados_lib.ocp_nlp_set.argtypes = \
                    [c_void_p, c_void_p, c_int, c_char_p, c_void_p]
                self.__acados_lib.ocp_nlp_set(self.nlp_config, \
                    self.nlp_solver, stage, field, value_data_p)
        return


    def cost_set(self, stage_: int, field_: str, value_, api='warn'):
        """
        Set numerical data in the cost module of the solver.

            :param stage: integer corresponding to shooting node
            :param field: string, e.g. 'yref', 'W', 'ext_cost_num_hess', 'zl', 'zu', 'Zl', 'Zu'
            :param value: of appropriate size
        """
        # cast value_ to avoid conversion issues
        if isinstance(value_, (float, int)):
            value_ = np.array([value_])

        if not isinstance(stage_, int):
            raise Exception('stage should be integer.')
        elif stage_ < 0 or stage_ > self.N:
            raise Exception(f'stage should be in [0, N], got {stage_}')

        value_ = value_.astype(float)

        field = field_
        field = field.encode('utf-8')

        stage = c_int(stage_)
        self.__acados_lib.ocp_nlp_cost_dims_get_from_attr.argtypes = \
            [c_void_p, c_void_p, c_void_p, c_int, c_char_p, POINTER(c_int)]
        self.__acados_lib.ocp_nlp_cost_dims_get_from_attr.restype = c_int

        dims = np.ascontiguousarray(np.zeros((2,)), dtype=np.intc)
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
                                    "Are you making an acados model from scrach? Add api='new' to cost_set and carry on.\n"
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

        self.__acados_lib.ocp_nlp_cost_model_set.argtypes = \
            [c_void_p, c_void_p, c_void_p, c_int, c_char_p, c_void_p]
        self.__acados_lib.ocp_nlp_cost_model_set(self.nlp_config, \
            self.nlp_dims, self.nlp_in, stage, field, value_data_p)

        return


    def constraints_set(self, stage_, field_, value_, api='warn'):
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

        field = field_
        field = field.encode('utf-8')

        stage = c_int(stage_)
        self.__acados_lib.ocp_nlp_constraint_dims_get_from_attr.argtypes = \
            [c_void_p, c_void_p, c_void_p, c_int, c_char_p, POINTER(c_int)]
        self.__acados_lib.ocp_nlp_constraint_dims_get_from_attr.restype = c_int

        dims = np.ascontiguousarray(np.zeros((2,)), dtype=np.intc)
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

        self.__acados_lib.ocp_nlp_constraints_model_set.argtypes = \
            [c_void_p, c_void_p, c_void_p, c_int, c_char_p, c_void_p]
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
        hess_block = scipy.linalg.block_diag(
            R_mat, Q_mat
        )
        nu = R_mat.shape[0]
        hess_block[nu:, :nu] = S_mat.T
        hess_block[:nu, nu:] = S_mat
        return hess_block


    def get_from_qp_in(self, stage_: int, field_: str):
        """
        Get numerical data from the current QP.

            :param stage: integer corresponding to shooting node
            :param field: string in ['A', 'B', 'b', 'Q', 'R', 'S', 'q', 'r', 'C', 'D', 'lg', 'ug', 'lbx', 'ubx', 'lbu', 'ubu']

        Note: additional supported fields are ['P', 'K', 'Lr'], which can be extracted form QP solver PARTIAL_CONDENSING_HPIPM.
        """
        # idx* should be added too..
        if not isinstance(stage_, int):
            raise TypeError("stage should be int")
        if stage_ > self.N:
            raise Exception("stage should be <= self.N")
        if field_ in self.__qp_dynamics_fields and stage_ >= self.N:
            raise ValueError(f"dynamics field {field_} not available at terminal stage")
        if field_ not in self.__qp_dynamics_fields + self.__qp_cost_fields + self.__qp_constraint_fields + self.__qp_pc_hpipm_fields:
            raise Exception(f"field {field_} not supported.")
        if field_ in self.__qp_pc_hpipm_fields:
            if self.acados_ocp.solver_options.qp_solver != "PARTIAL_CONDENSING_HPIPM" or self.acados_ocp.solver_options.qp_solver_cond_N != self.acados_ocp.dims.N:
                raise Exception(f"field {field_} only works for PARTIAL_CONDENSING_HPIPM QP solver with qp_solver_cond_N == N.")
            if field_ in ["P", "K", "p"] and stage_ == 0 and self.acados_ocp.dims.nbxe_0 > 0:
                raise Exception(f"getting field {field_} at stage 0 only works without x0 elimination (see nbxe_0).")

        field = field_.encode('utf-8')
        stage = c_int(stage_)

        # get dims
        self.__acados_lib.ocp_nlp_qp_dims_get_from_attr.argtypes = \
            [c_void_p, c_void_p, c_void_p, c_int, c_char_p, POINTER(c_int)]
        self.__acados_lib.ocp_nlp_qp_dims_get_from_attr.restype = c_int

        dims = np.ascontiguousarray(np.zeros((2,)), dtype=np.intc)
        dims_data = cast(dims.ctypes.data, POINTER(c_int))

        self.__acados_lib.ocp_nlp_qp_dims_get_from_attr(self.nlp_config, \
            self.nlp_dims, self.nlp_out, stage_, field, dims_data)

        # create output data
        out = np.ascontiguousarray(np.zeros((np.prod(dims),)), dtype=np.float64)
        out = out.reshape(dims[0], dims[1], order='F')

        out_data = cast(out.ctypes.data, POINTER(c_double))
        out_data_p = cast((out_data), c_void_p)

        # call getter
        self.__acados_lib.ocp_nlp_get_at_stage.argtypes = \
            [c_void_p, c_void_p, c_void_p, c_int, c_char_p, c_void_p]
        self.__acados_lib.ocp_nlp_get_at_stage(self.nlp_config, \
            self.nlp_dims, self.nlp_solver, stage, field, out_data_p)

        return out


    def options_set(self, field_, value_):
        """
        Set options of the solver.

            :param field: string, e.g. 'print_level', 'rti_phase', 'initialize_t_slacks', 'step_length', 'alpha_min', 'alpha_reduction', 'qp_warm_start', 'line_search_use_sufficient_descent', 'full_step_dual', 'globalization_use_SOC', 'qp_tol_stat', 'qp_tol_eq', 'qp_tol_ineq', 'qp_tol_comp', 'qp_tau_min', 'qp_mu0'

            :param value: of type int, float, string

            - qp_tol_stat: QP solver tolerance stationarity
            - qp_tol_eq: QP solver tolerance equalities
            - qp_tol_ineq: QP solver tolerance inequalities
            - qp_tol_comp: QP solver tolerance complementarity
            - qp_tau_min: for HPIPM QP solvers: minimum value of barrier parameter in HPIPM
            - qp_mu0: for HPIPM QP solvers: initial value for complementarity slackness
            - warm_start_first_qp: indicates if first QP in SQP is warm_started
        """
        int_fields = ['print_level', 'rti_phase', 'initialize_t_slacks', 'qp_warm_start',
                      'line_search_use_sufficient_descent', 'full_step_dual', 'globalization_use_SOC', 'warm_start_first_qp']
        double_fields = ['step_length', 'tol_eq', 'tol_stat', 'tol_ineq', 'tol_comp', 'alpha_min', 'alpha_reduction',
                         'eps_sufficient_descent', 'qp_tol_stat', 'qp_tol_eq', 'qp_tol_ineq', 'qp_tol_comp', 'qp_tau_min', 'qp_mu0']
        string_fields = ['globalization']

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

        elif field_ in string_fields:
            if not isinstance(value_, str):
                raise Exception(f'solver option \'{field_}\' must be of type str. You have {type(value_)}.')
            else:
                value_ctypes = value_.encode('utf-8')
        else:
            fields = ', '.join(int_fields + double_fields + string_fields)
            raise Exception(f'AcadosOcpSolver.options_set() does not support field \'{field_}\'.\n'\
                f' Possible values are {fields}.')


        if field_ == 'rti_phase':
            if value_ < 0 or value_ > 2:
                raise Exception('AcadosOcpSolver.options_set(): argument \'rti_phase\' can '
                    'take only values 0, 1, 2 for SQP-RTI-type solvers')
            if self.solver_options['nlp_solver_type'] != 'SQP_RTI' and value_ > 0:
                raise Exception('AcadosOcpSolver.options_set(): argument \'rti_phase\' can '
                    'take only value 0 for SQP-type solvers')

        # encode
        field = field_
        field = field.encode('utf-8')

        # call C interface
        if field_ in string_fields:
            self.__acados_lib.ocp_nlp_solver_opts_set.argtypes = \
                [c_void_p, c_void_p, c_char_p, c_char_p]
            self.__acados_lib.ocp_nlp_solver_opts_set(self.nlp_config, \
                self.nlp_opts, field, value_ctypes)
        else:
            self.__acados_lib.ocp_nlp_solver_opts_set.argtypes = \
                [c_void_p, c_void_p, c_char_p, c_void_p]
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

        if any(idx_values_ >= self.acados_ocp.dims.np):
            raise Exception(f'idx_values_ contains value >= np = {self.acados_ocp.dims.np}')

        stage = c_int(stage_)
        n_update = c_int(len(param_values_))

        param_data = cast(param_values_.ctypes.data, POINTER(c_double))
        c_idx_values = np.ascontiguousarray(idx_values_, dtype=np.intc)
        idx_data = cast(c_idx_values.ctypes.data, POINTER(c_int))

        getattr(self.shared_lib, f"{self.name}_acados_update_params_sparse").argtypes = \
                        [c_void_p, c_int, POINTER(c_int), POINTER(c_double), c_int]
        getattr(self.shared_lib, f"{self.name}_acados_update_params_sparse").restype = c_int
        getattr(self.shared_lib, f"{self.name}_acados_update_params_sparse") \
                                    (self.capsule, stage, idx_data, param_data, n_update)

    def __del__(self):
        if self.solver_created:
            getattr(self.shared_lib, f"{self.name}_acados_free").argtypes = [c_void_p]
            getattr(self.shared_lib, f"{self.name}_acados_free").restype = c_int
            getattr(self.shared_lib, f"{self.name}_acados_free")(self.capsule)

            getattr(self.shared_lib, f"{self.name}_acados_free_capsule").argtypes = [c_void_p]
            getattr(self.shared_lib, f"{self.name}_acados_free_capsule").restype = c_int
            getattr(self.shared_lib, f"{self.name}_acados_free_capsule")(self.capsule)

            try:
                self.dlclose(self.shared_lib._handle)
            except:
                print(f"WARNING: acados Python interface could not close shared_lib handle of AcadosOcpSolver {self.name}.\n",
                     "Attempting to create a new one with the same name will likely result in the old one being used!")
                pass
