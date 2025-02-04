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

from typing import Union, List
import numpy as np
import casadi as ca
from copy import deepcopy

import os, json

from .acados_model import AcadosModel
from .acados_dims import AcadosOcpDims
from .acados_ocp_cost import AcadosOcpCost
from .acados_ocp_constraints import AcadosOcpConstraints
from .acados_ocp_options import AcadosOcpOptions, INTEGRATOR_TYPES, COLLOCATION_TYPES, COST_DISCRETIZATION_TYPES
from .acados_ocp import AcadosOcp
from .casadi_function_generation import GenerateContext
from .utils import make_object_json_dumpable, get_acados_path, format_class_dict, get_shared_lib_ext, render_template, is_empty


def find_non_default_fields_of_obj(obj: Union[AcadosOcpCost, AcadosOcpConstraints, AcadosOcpOptions], stage_type='all') -> list:

    all_fields = [field for field in dir(obj) if not field.startswith("_")]

    # remove special properties which are translated to other fields
    if isinstance(obj, AcadosOcpConstraints):
        all_fields.remove('x0')
        all_fields = [field for field in all_fields if not field.startswith("J")]

    if isinstance(obj, AcadosOcpOptions):
        all_fields.remove('qp_tol')
        all_fields.remove('tol')

    all_fields = [field for field in all_fields if not callable(getattr(obj, field))]

    if stage_type == 'all':
        pass
    elif stage_type == 'initial':
        all_fields = [field for field in all_fields if field.endswith("_0")]
    elif stage_type == 'terminal':
        all_fields = [field for field in all_fields if field.endswith("_e")]
    else:
        raise Exception(f"stage_type {stage_type} not supported.")

    obj_type = type(obj)
    dummy_obj = obj_type()
    nondefault_fields = []
    for field in all_fields:
        val = getattr(obj, field)
        default_val = getattr(dummy_obj, field)
        if not isinstance(val, type(default_val)):
            nondefault_fields.append(field)

        elif isinstance(val, np.ndarray):
            if not np.array_equal(val, default_val):
                nondefault_fields.append(field)

        elif val != default_val:
            nondefault_fields.append(field)

    return nondefault_fields


class AcadosMultiphaseOptions:
    """
    Class containing options that might be different for each phase.

    All of the fields can be either None, then the corresponding value from ocp.solver_options is used,
    or a list of length n_phases describing the value for this option at each phase.

    - integrator_type: list of strings, must be in ["ERK", "IRK", "GNSF", "DISCRETE", "LIFTED_IRK"]
    - collocation_type: list of strings, must be in ["GAUSS_RADAU_IIA", "GAUSS_LEGENDRE", "EXPLICIT_RUNGE_KUTTA"]
    - cost_discretization: list of strings, must be in ["EULER", "INTEGRATOR"]
    """
    def __init__(self):
        self.integrator_type = None
        self.collocation_type = None
        self.cost_discretization = None

    def make_consistent(self, opts: AcadosOcpOptions, n_phases: int) -> None:
        for field, variants in zip(['integrator_type', 'collocation_type', 'cost_discretization'],
                                [INTEGRATOR_TYPES, COLLOCATION_TYPES, COST_DISCRETIZATION_TYPES]
                                ):
            if getattr(self, field) is None:
                # non varying field, use value from ocp opts
                setattr(self, field, [getattr(opts, field) for _ in range(n_phases)])
            if not isinstance(getattr(self, field), list):
                raise Exception(f'AcadosMultiphaseOptions.{field} must be a list, got {getattr(self, field)}.')
            if len(getattr(self, field)) != n_phases:
                raise Exception(f'AcadosMultiphaseOptions.{field} must be a list of length n_phases, got {getattr(self, field)}.')
            if not all([item in variants for item in getattr(self, field)]):
                raise Exception(f'AcadosMultiphaseOptions.{field} must be a list of strings in {variants}, got {getattr(self, field)}.')


class AcadosMultiphaseOcp:
    """
    Class containing the description of an optimal control problem with multiple phases.
    This object can be used to create an :py:class:`acados_template.acados_ocp_solver.AcadosOcpSolver`.

    Initial cost and constraints are defined by the first phase, terminal cost and constraints by the last phase.
    All other phases are treated as intermediate phases, where only dynamics and path cost and constraints are used.

    Solver options are shared between all phases. Options that can vary phase-wise must be set via self.mocp_opts of type :py:class:`acados_template.acados_multiphase_ocp.AcadosMultiphaseOptions`.

    :param N_list: list containing the number of shooting intervals for each phase
    """
    def __init__(self, N_list: list):

        if not isinstance(N_list, list) or len(N_list) < 1:
            raise Exception("N_list must be a list of integers.")
        if any([not isinstance(N, int) for N in N_list]):
            raise Exception("N_list must be a list of integers.")

        n_phases = len(N_list)

        self.n_phases = n_phases
        self.N_list = N_list

        self.name = 'multiphase_ocp'
        self.model = [AcadosModel() for _ in range(n_phases)]
        """Model definitions, type :py:class:`acados_template.acados_model.AcadosModel`"""
        self.cost = [AcadosOcpCost() for _ in range(n_phases)]
        """Cost definitions, type :py:class:`acados_template.acados_ocp_cost.AcadosOcpCost`"""
        self.constraints = [AcadosOcpConstraints() for _ in range(n_phases)]
        """Constraints definitions, type :py:class:`acados_template.acados_ocp_constraints.AcadosOcpConstraints`"""

        self.phases_dims = [AcadosOcpDims() for _ in range(n_phases)]

        self.dummy_ocp_list: List[AcadosOcp] = []

        # NOTE: this is the same for AcadosOcp
        self.solver_options = AcadosOcpOptions()
        """Solver Options, type :py:class:`acados_template.acados_ocp_options.AcadosOcpOptions`"""
        self.mocp_opts = AcadosMultiphaseOptions()
        """Phase-wise varying solver Options, type :py:class:`acados_template.acados_multiphase_ocp.AcadosMultiphaseOptions`"""

        acados_path = get_acados_path()

        self.acados_include_path = os.path.join(acados_path, 'include').replace(os.sep, '/') # the replace part is important on Windows for CMake
        """Path to acados include directory (set automatically), type: `string`"""
        self.acados_lib_path = os.path.join(acados_path, 'lib').replace(os.sep, '/') # the replace part is important on Windows for CMake
        """Path to where acados library is located, type: `string`"""
        self.shared_lib_ext = get_shared_lib_ext()

        # get cython paths
        from sysconfig import get_paths
        self.cython_include_dirs = [np.get_include(), get_paths()['include']]

        self.__parameter_values = [np.array([]) for _ in range(n_phases)]
        self.__p_global_values = np.array([])
        self.__problem_class = "MOCP"
        self.__json_file = 'mocp.json'

        self.code_export_directory = 'c_generated_code'
        """Path to where code will be exported. Default: `c_generated_code`."""

        self.simulink_opts = None
        """Options to configure Simulink S-function blocks, mainly to activate possible Inputs and Outputs."""


    @property
    def parameter_values(self):
        """:math:`p` - list of initial values for parameter vector.
        Type: `list` of `numpy.ndarray` of shape `(np_i, )`.
        - can be updated stagewise."""
        return self.__parameter_values

    @parameter_values.setter
    def parameter_values(self, parameter_values):
        if not isinstance(parameter_values, list):
            raise Exception('parameter_values must be a list of numpy.ndarrays.')
        elif len(parameter_values) != self.n_phases:
            raise Exception('parameter_values must be a list of length n_phases.')
        self.__parameter_values = parameter_values


    @property
    def p_global_values(self):
        r"""initial values for :math:`p_\text{global}` vector, see `AcadosModel.p_global` - can be updated.
        NOTE: `p_global` is shared between all phases.
        Type: `numpy.ndarray` of shape `(np_global, )`.
        """
        return self.__p_global_values

    @p_global_values.setter
    def p_global_values(self, p_global_values):
        if not isinstance(p_global_values, np.ndarray):
            raise Exception('p_global_values must be a single numpy.ndarrays.')
        self.__p_global_values = p_global_values


    @property
    def json_file(self):
        """Name of the json file where the problem description is stored."""
        return self.__json_file

    @json_file.setter
    def json_file(self, json_file):
        self.__json_file = json_file

    def set_phase(self, ocp: AcadosOcp, phase_idx: int) -> None:
        """
        Set phase of the multiphase OCP to match the given OCP.

        NOTE: model, cost, constraints and parameter_values are taken from phase OCP,
              all other fields, especially options are ignored.

        :param ocp: OCP to be set as phase
        :param phase_idx: index of the phase, must be in [0, n_phases-1]
        """
        if phase_idx >= self.n_phases:
            raise Exception(f"phase_idx {phase_idx} out of bounds, must be in [0, {self.n_phases-1}].")

        # check options
        non_default_opts = find_non_default_fields_of_obj(ocp.solver_options)
        if len(non_default_opts) > 0:
            print(f"WARNING: set_phase: Phase {phase_idx} contains non-default solver options: {non_default_opts}, which will be ignored.\n",
                   "Solver options need to be set via AcadosMultiphaseOcp.solver_options or mocp_opts instead.")

        # set phase
        self.model[phase_idx] = ocp.model
        self.cost[phase_idx] = ocp.cost
        self.constraints[phase_idx] = ocp.constraints
        self.parameter_values[phase_idx] = ocp.parameter_values

        if ocp.p_global_values.size > 0:
            print(f"WARNING: set_phase: Phase {phase_idx} contains p_global_values which will be ignored.")

        return

    def make_consistent(self) -> None:

        self.N_horizon = sum(self.N_list)
        self.solver_options.N_horizon = self.N_horizon # NOTE: to not change options when making ocp consistent

        # check options
        self.mocp_opts.make_consistent(self.solver_options, n_phases=self.n_phases)

        # check phases formulation objects are distinct
        warning = "\nNOTE: this can happen if set_phase() is called with the same ocp object for multiple phases."
        for field in ['model', 'cost', 'constraints']:
            if len(set(getattr(self, field))) != self.n_phases:
                raise Exception(f"AcadosMultiphaseOcp: make_consistent: {field} objects are not distinct.{warning}")

        # p_global check:
        p_global = self.model[0].p_global
        for i in range(self.n_phases):
            if is_empty(p_global) and not is_empty(self.model[i].p_global):
                raise Exception(f"p_global is empty for phase 0, but not for phase {i}. Should be the same for all phases.")
            if not is_empty(p_global) and not ca.is_equal(p_global, self.model[i].p_global):
                raise Exception(f"p_global is different for phase 0 and phase {i}. Should be the same for all phases.")

        # compute phase indices
        phase_idx = np.cumsum([0] + self.N_list).tolist()

        self.start_idx = phase_idx[:-1]
        self.end_idx = phase_idx[1:]

        self.cost_start_idx = phase_idx.copy()
        self.cost_start_idx[0] += 1

        # make model names unique if necessary
        model_name_list = [self.model[i].name for i in range(self.n_phases)]
        n_names = len(set(model_name_list))
        if n_names != self.n_phases:
            print(f"model names are not unique: got {model_name_list}")
            print("adding _i to model names")
            for i in range(self.n_phases):
                self.model[i].name = f"{self.model[i].name}_{i}"
            model_name_list = [self.model[i].name for i in range(self.n_phases)]
            print(f"new model names are {model_name_list}")

        # make phase OCPs consistent, warn about unused fields
        for i in range(self.n_phases):
            ocp = AcadosOcp()
            ocp.dims = self.phases_dims[i]
            ocp.model = self.model[i]
            ocp.constraints = self.constraints[i]
            ocp.cost = self.cost[i]
            ocp.parameter_values = self.parameter_values[i]
            ocp.p_global_values = self.p_global_values
            ocp.solver_options = self.solver_options

            # set phase dependent options
            ocp.solver_options.integrator_type = self.mocp_opts.integrator_type[i]
            ocp.solver_options.collocation_type = self.mocp_opts.collocation_type[i]
            ocp.solver_options.cost_discretization = self.mocp_opts.cost_discretization[i]

            if i != self.n_phases - 1:
                nondefault_fields = []
                nondefault_fields += find_non_default_fields_of_obj(ocp.cost, stage_type='terminal')
                nondefault_fields += find_non_default_fields_of_obj(ocp.constraints, stage_type='terminal')
                nondefault_fields += find_non_default_fields_of_obj(ocp.model, stage_type='terminal')
                if len(nondefault_fields) > 0:
                    print(f"Phase {i} contains non-default terminal fields: {nondefault_fields}, which will be ignored.")
            elif i != 0:
                nondefault_fields = []
                nondefault_fields += find_non_default_fields_of_obj(ocp.cost, stage_type='initial')
                nondefault_fields += find_non_default_fields_of_obj(ocp.constraints, stage_type='initial')
                nondefault_fields += find_non_default_fields_of_obj(ocp.model, stage_type='initial')
                if len(nondefault_fields) > 0:
                    print(f"Phase {i} contains non-default initial fields: {nondefault_fields}, which will be ignored.")

            print(f"Calling make_consistent for phase {i}.")
            ocp.make_consistent(is_mocp_phase=True)

            self.dummy_ocp_list.append(ocp)

        # check for transition consistency
        nx_list = [self.phases_dims[i].nx for i in range(self.n_phases)]
        for i in range(1, self.n_phases):
            if nx_list[i] != nx_list[i-1]:
                if self.phases_dims[i].nx != self.phases_dims[i-1].nx_next:
                    raise Exception(f"detected stage transition with different nx from phase {i-1} to {i}, nx_next at phase {i-1} = {self.phases_dims[i-1].nx_next} should match nx at phase {i} = {nx_list[i]}.")
                if self.N_list[i-1] != 1 or self.mocp_opts.integrator_type[i-1] != 'DISCRETE':
                    raise Exception(f"detected stage transition with different nx from phase {i-1} to {i}, which is only supported for integrator_type='DISCRETE' and N_list[i] == 1.")
        return


    def to_dict(self) -> dict:
        # Copy ocp object dictionary
        ocp_dict = dict(deepcopy(self).__dict__)
        del ocp_dict['dummy_ocp_list']

        # convert acados classes to dicts
        for key, v in ocp_dict.items():
            if isinstance(v, (AcadosOcpOptions, AcadosMultiphaseOptions)):
                ocp_dict[key]=dict(getattr(self, key).__dict__)
            if isinstance(v, list):
                for i, item in enumerate(v):
                    if isinstance(item, (AcadosModel, AcadosOcpDims, AcadosOcpConstraints, AcadosOcpCost)):
                        ocp_dict[key][i] = format_class_dict(dict(item.__dict__))

        ocp_dict = format_class_dict(ocp_dict)

        # delete keys that should not be used
        del ocp_dict['solver_options']['integrator_type']
        del ocp_dict['solver_options']['collocation_type']
        del ocp_dict['solver_options']['cost_discretization']

        return ocp_dict


    def dump_to_json(self) -> None:
        ocp_nlp_dict = self.to_dict()
        with open(self.json_file, 'w') as f:
            json.dump(ocp_nlp_dict, f, default=make_object_json_dumpable, indent=4, sort_keys=True)
        return


    def __get_template_list(self, cmake_builder=None) -> list:
        """
        returns a list of tuples in the form:
        (input_filename, output_filname)
        or
        (input_filename, output_filname, output_directory)
        """
        name = self.name
        template_list = []

        template_list.append(('main_multi.in.c', f'main_{name}.c'))
        template_list.append(('acados_multi_solver.in.h', f'acados_solver_{name}.h'))
        template_list.append(('acados_multi_solver.in.c', f'acados_solver_{name}.c'))
        # template_list.append(('acados_solver.in.pxd', f'acados_solver.pxd'))
        if cmake_builder is not None:
            template_list.append(('multi_CMakeLists.in.txt', 'CMakeLists.txt'))
        else:
            template_list.append(('multi_Makefile.in', 'Makefile'))

        if self.phases_dims[0].np_global > 0:
            template_list.append(('p_global_precompute_fun.in.h', f'{self.name}_p_global_precompute_fun.h'))

        # Simulink
        if self.simulink_opts is not None:
            raise NotImplementedError('Simulink not yet supported for multiphase OCPs.')

        return template_list


    def render_templates(self, cmake_builder=None):

        # model templates
        for i, dummy_ocp in enumerate(self.dummy_ocp_list):
            # this is the only option that can vary and influence external functions to be generated
            dummy_ocp.solver_options.integrator_type = self.mocp_opts.integrator_type[i]

            template_list = dummy_ocp._get_external_function_header_templates()
            # dump dummy_ocp
            dummy_ocp.json_file = 'tmp_ocp.json'
            dummy_ocp.dump_to_json()
            tmp_json_path = os.path.abspath(dummy_ocp.json_file)

            # render templates
            for tup in template_list:
                output_dir = self.code_export_directory if len(tup) <= 2 else tup[2]
                render_template(tup[0], tup[1], output_dir, tmp_json_path)

        print("rendered model templates successfully")

        # check json file
        json_path = os.path.abspath(self.json_file)
        if not os.path.exists(json_path):
            raise Exception(f'Path "{json_path}" not found!')

        # solver templates
        template_list = self.__get_template_list(cmake_builder=cmake_builder)

        # Render templates
        for tup in template_list:
            output_dir = self.code_export_directory if len(tup) <= 2 else tup[2]
            render_template(tup[0], tup[1], output_dir, json_path)

        # # Custom templates
        # acados_template_path = os.path.dirname(os.path.abspath(__file__))
        # custom_template_glob = os.path.join(acados_template_path, 'custom_update_templates', '*')
        # for tup in ocp.solver_options.custom_templates:
        #     render_template(tup[0], tup[1], ocp.code_export_directory, json_path, template_glob=custom_template_glob)
        print("\nmocp_render_templates: rendered solver templates successfully!\n")

        return


    def generate_external_functions(self) -> GenerateContext:

        # options for code generation
        code_gen_opts = dict()
        code_gen_opts['generate_hess'] = self.solver_options.hessian_approx == 'EXACT'
        code_gen_opts['with_solution_sens_wrt_params'] = self.solver_options.with_solution_sens_wrt_params
        code_gen_opts['with_value_sens_wrt_params'] = self.solver_options.with_value_sens_wrt_params
        code_gen_opts['code_export_directory'] = self.code_export_directory
        code_gen_opts['ext_fun_expand'] = self.solver_options.ext_fun_expand

        context = GenerateContext(self.model[0].p_global, self.name, code_gen_opts)

        for i in range(self.n_phases):
            ignore_initial = True if i != 0 else False
            ignore_terminal = True if i != self.n_phases-1 else False
            # this is the only option that can vary and influence external functions to be generated
            self.dummy_ocp_list[i].solver_options.integrator_type = self.mocp_opts.integrator_type[i]
            context = self.dummy_ocp_list[i]._setup_code_generation_context(context, ignore_initial, ignore_terminal)
            self.dummy_ocp_list[i].code_export_directory = self.code_export_directory

        context.finalize()
        self.__external_function_files_model = context.get_external_function_file_list(ocp_specific=False)
        self.__external_function_files_ocp = context.get_external_function_file_list(ocp_specific=True)
        for i in range(self.n_phases):
            self.phases_dims[i].n_global_data = context.get_n_global_data()

        return context
