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

import os, json
import numpy as np
from copy import deepcopy
from .acados_model import AcadosModel
from .acados_dims import AcadosSimDims
from .builders import CMakeBuilder
from .utils import (get_acados_path, get_shared_lib_ext, format_class_dict, check_casadi_version,
                    make_object_json_dumpable, render_template)
from .casadi_function_generation import (
                    GenerateContext,
                    AcadosCodegenOptions,
                    generate_c_code_explicit_ode,
                    generate_c_code_gnsf,
                    generate_c_code_implicit_ode)

class AcadosSimOptions:
    """
    class containing the solver options
    """
    def __init__(self):
        self.__integrator_type = 'ERK'
        self.__collocation_type = 'GAUSS_LEGENDRE'
        self.__Tsim = None
        # ints
        self.__sim_method_num_stages = 4
        self.__sim_method_num_steps = 1
        self.__sim_method_newton_iter = 3
        # doubles
        self.__sim_method_newton_tol = 0.0
        # bools
        self.__sens_forw = True
        self.__sens_adj = False
        self.__sens_algebraic = False
        self.__sens_hess = False
        self.__output_z = True
        self.__sim_method_jac_reuse = 0
        env = os.environ
        self.__ext_fun_compile_flags = '-O2' if 'ACADOS_EXT_FUN_COMPILE_FLAGS' not in env else env['ACADOS_EXT_FUN_COMPILE_FLAGS']
        self.__ext_fun_expand_dyn = False
        self.__num_threads_in_batch_solve: int = 1
        self.__with_batch_functionality: bool = False

    @property
    def integrator_type(self):
        """Integrator type. Default: 'ERK'."""
        return self.__integrator_type

    @property
    def num_stages(self):
        """Number of stages in the integrator. Default: 4"""
        return self.__sim_method_num_stages

    @property
    def num_steps(self):
        """Number of steps in the integrator. Default: 1"""
        return self.__sim_method_num_steps

    @property
    def newton_iter(self):
        """Number of Newton iterations in simulation method. Default: 3"""
        return self.__sim_method_newton_iter

    @property
    def newton_tol(self):
        """
        Tolerance for Newton system solved in implicit integrator (IRK, GNSF).
        0.0 means this is not used and exactly newton_iter iterations are carried out.
        Default: 0.0
        """
        return self.__sim_method_newton_tol

    @property
    def sens_forw(self):
        """Boolean determining if forward sensitivities are computed. Default: True"""
        return self.__sens_forw

    @property
    def sens_adj(self):
        """Boolean determining if adjoint sensitivities are computed. Default: False"""
        return self.__sens_adj

    @property
    def sens_algebraic(self):
        """Boolean determining if sensitivities wrt algebraic variables are computed. Default: False"""
        return self.__sens_algebraic

    @property
    def sens_hess(self):
        """Boolean determining if hessians are computed. Default: False"""
        return self.__sens_hess

    @property
    def output_z(self):
        """Boolean determining if values for algebraic variables (corresponding to start of simulation interval) are computed. Default: True"""
        return self.__output_z

    @property
    def sim_method_jac_reuse(self):
        """Integer determining if jacobians are reused (0 or 1). Default: 0"""
        return self.__sim_method_jac_reuse

    @property
    def T(self):
        """Time horizon"""
        return self.__Tsim

    @property
    def collocation_type(self):
        """Collocation type: relevant for implicit integrators
        -- string in {'GAUSS_RADAU_IIA', 'GAUSS_LEGENDRE', 'EXPLICIT_RUNGE_KUTTA'}.

        Default: GAUSS_LEGENDRE.

        .. note:: GAUSS_LEGENDRE tableaus yield integration methods that are A-stable, but not L-stable and have order `2 * num_stages`,
        .. note:: GAUSS_RADAU_IIA tableaus yield integration methods that are L-stable and have order `2 * num_stages - 1`.
        .. note:: EXPLICIT_RUNGE_KUTTA tableaus can be used for comparisons of ERK and IRK to ensure correctness, but are only recommended with ERK for users.
        """
        return self.__collocation_type

    @property
    def ext_fun_compile_flags(self):
        """
        String with compiler flags for external function compilation.
        Default: '-O2' if environment variable ACADOS_EXT_FUN_COMPILE_FLAGS is not set, else ACADOS_EXT_FUN_COMPILE_FLAGS is used as default.
        """
        return self.__ext_fun_compile_flags


    @property
    def ext_fun_expand_dyn(self):
        """
        Flag indicating whether CasADi.MX should be expanded to CasADi.SX before code generation.
        Default: False
        """
        return self.__ext_fun_expand_dyn


    @property
    def num_threads_in_batch_solve(self):
        """
        DEPRECATED, use the flag with_batch_functionality instead and pass the number of threads directly to the BatchSolver.
        Integer indicating how many threads should be used within the batch solve.
        If more than one thread should be used, the sim solver is compiled with openmp.
        Default: 1.
        """
        return self.__num_threads_in_batch_solve

    @property
    def with_batch_functionality(self):
        """
        Whether the AcadosSimBatchSolver can be used.
        In this case, the sim solver is compiled with openmp.
        Default: False.
        """
        return self.__with_batch_functionality

    @ext_fun_compile_flags.setter
    def ext_fun_compile_flags(self, ext_fun_compile_flags):
        if isinstance(ext_fun_compile_flags, str):
            self.__ext_fun_compile_flags = ext_fun_compile_flags
        else:
            raise TypeError('Invalid ext_fun_compile_flags value, expected a string.\n')

    @ext_fun_expand_dyn.setter
    def ext_fun_expand_dyn(self, ext_fun_expand_dyn):
        if isinstance(ext_fun_expand_dyn, bool):
            self.__ext_fun_expand_dyn = ext_fun_expand_dyn
        else:
            raise TypeError('Invalid ext_fun_expand_dyn value, expected bool.\n')

    @integrator_type.setter
    def integrator_type(self, integrator_type):
        integrator_types = ('ERK', 'IRK', 'GNSF')
        if integrator_type in integrator_types:
            self.__integrator_type = integrator_type
        else:
            raise ValueError('Invalid integrator_type value. Possible values are:\n\n' \
                    + ',\n'.join(integrator_types) + '.\n\nYou have: ' + integrator_type + '.\n\n')

    @collocation_type.setter
    def collocation_type(self, collocation_type):
        collocation_types = ('GAUSS_RADAU_IIA', 'GAUSS_LEGENDRE')
        if collocation_type in collocation_types:
            self.__collocation_type = collocation_type
        else:
            raise ValueError('Invalid collocation_type value. Possible values are:\n\n' \
                    + ',\n'.join(collocation_types) + '.\n\nYou have: ' + collocation_type + '.\n\n')

    @T.setter
    def T(self, T):
        self.__Tsim = T

    @num_stages.setter
    def num_stages(self, num_stages):
        if isinstance(num_stages, int):
            self.__sim_method_num_stages = num_stages
        else:
            raise ValueError('Invalid num_stages value. num_stages must be an integer.')

    @num_steps.setter
    def num_steps(self, num_steps):
        if isinstance(num_steps, int):
            self.__sim_method_num_steps = num_steps
        else:
            raise TypeError('Invalid num_steps value. num_steps must be an integer.')

    @newton_iter.setter
    def newton_iter(self, newton_iter):
        if isinstance(newton_iter, int):
            self.__sim_method_newton_iter = newton_iter
        else:
            raise TypeError('Invalid newton_iter value. newton_iter must be an integer.')

    @newton_tol.setter
    def newton_tol(self, newton_tol):
        if isinstance(newton_tol, float):
            self.__sim_method_newton_tol = newton_tol
        else:
            raise TypeError('Invalid newton_tol value. newton_tol must be a float.')

    @sens_forw.setter
    def sens_forw(self, sens_forw):
        if sens_forw in (True, False):
            self.__sens_forw = sens_forw
        else:
            raise ValueError('Invalid sens_forw value. sens_forw must be a Boolean.')

    @sens_adj.setter
    def sens_adj(self, sens_adj):
        if sens_adj in (True, False):
            self.__sens_adj = sens_adj
        else:
            raise ValueError('Invalid sens_adj value. sens_adj must be a Boolean.')

    @sens_hess.setter
    def sens_hess(self, sens_hess):
        if sens_hess in (True, False):
            self.__sens_hess = sens_hess
        else:
            raise ValueError('Invalid sens_hess value. sens_hess must be a Boolean.')

    @sens_algebraic.setter
    def sens_algebraic(self, sens_algebraic):
        if sens_algebraic in (True, False):
            self.__sens_algebraic = sens_algebraic
        else:
            raise ValueError('Invalid sens_algebraic value. sens_algebraic must be a Boolean.')

    @output_z.setter
    def output_z(self, output_z):
        if output_z in (True, False):
            self.__output_z = output_z
        else:
            raise ValueError('Invalid output_z value. output_z must be a Boolean.')

    @sim_method_jac_reuse.setter
    def sim_method_jac_reuse(self, sim_method_jac_reuse):
        if sim_method_jac_reuse in (0, 1):
            self.__sim_method_jac_reuse = sim_method_jac_reuse
        else:
            raise ValueError('Invalid sim_method_jac_reuse value. sim_method_jac_reuse must be 0 or 1.')

    @num_threads_in_batch_solve.setter
    def num_threads_in_batch_solve(self, num_threads_in_batch_solve):
        print("Warning: num_threads_in_batch_solve is deprecated, set the flag with_batch_functionality instead and pass the number of threads directly to the BatchSolver.")
        if isinstance(num_threads_in_batch_solve, int) and num_threads_in_batch_solve > 0:
            self.__num_threads_in_batch_solve = num_threads_in_batch_solve
        else:
            raise ValueError('Invalid num_threads_in_batch_solve value. num_threads_in_batch_solve must be a positive integer.')

    @with_batch_functionality.setter
    def with_batch_functionality(self, with_batch_functionality):
        if isinstance(with_batch_functionality, bool):
            self.__with_batch_functionality = with_batch_functionality
        else:
            raise Exception('Invalid with_batch_functionality value. Expected bool.')

class AcadosSim:
    """
    The class has the following properties that can be modified to formulate a specific simulation problem, see below:

    :param acados_path: string with the path to acados. It is used to generate the include and lib paths.

    - :py:attr:`dims` of type :py:class:`acados_template.acados_dims.AcadosSimDims` - are automatically detected from model
    - :py:attr:`model` of type :py:class:`acados_template.acados_model.AcadosModel`
    - :py:attr:`solver_options` of type :py:class:`acados_template.acados_sim.AcadosSimOptions`

    - :py:attr:`acados_include_path` (set automatically)
    - :py:attr:`shared_lib_ext` (set automatically)
    - :py:attr:`acados_lib_path` (set automatically)
    - :py:attr:`parameter_values` - used to initialize the parameters (can be changed)

    """
    def __init__(self, acados_path=''):
        if acados_path == '':
            acados_path = get_acados_path()
        self.dims = AcadosSimDims()
        """Dimension definitions, automatically detected from :py:attr:`model`. Type :py:class:`acados_template.acados_dims.AcadosSimDims`"""
        self.model = AcadosModel()
        """Model definitions, type :py:class:`acados_template.acados_model.AcadosModel`"""
        self.solver_options = AcadosSimOptions()
        """Solver Options, type :py:class:`acados_template.acados_sim.AcadosSimOptions`"""

        self.acados_include_path = os.path.join(acados_path, 'include').replace(os.sep, '/') # the replace part is important on Windows for CMake
        """Path to acados include directory (set automatically), type: `string`"""
        self.acados_lib_path = os.path.join(acados_path, 'lib').replace(os.sep, '/') # the replace part is important on Windows for CMake
        """Path to where acados library is located (set automatically), type: `string`"""

        self.code_export_directory = 'c_generated_code'
        """Path to where code will be exported. Default: `c_generated_code`."""
        self.shared_lib_ext = get_shared_lib_ext()

        self.simulink_opts = None
        """Options to configure Simulink S-function blocks, if not None, MATLAB related files will be generated. More options may be added in the future, similar to OCP interface"""

        # get cython paths
        from sysconfig import get_paths
        self.cython_include_dirs = [np.get_include(), get_paths()['include']]

        self.__parameter_values = np.array([])
        self.__problem_class = 'SIM'

    @property
    def parameter_values(self):
        """:math:`p` - initial values for parameter - can be updated"""
        return self.__parameter_values

    @parameter_values.setter
    def parameter_values(self, parameter_values):
        if isinstance(parameter_values, np.ndarray):
            self.__parameter_values = parameter_values
        else:
            raise ValueError('Invalid parameter_values value. ' +
                            f'Expected numpy array, got {type(parameter_values)}.')

    def make_consistent(self):
        self.model.make_consistent(self.dims)

        if self.parameter_values.shape[0] != self.dims.np:
            raise ValueError('inconsistent dimension np, regarding model.p and parameter_values.' + \
                f'\nGot np = {self.dims.np}, acados_sim.parameter_values.shape = {self.parameter_values.shape[0]}\n')

        # check required arguments are given
        if self.solver_options.T is None:
            raise ValueError('acados_sim.solver_options.T is None, should be provided.')


    def to_dict(self) -> dict:
        # Copy input sim object dictionary
        sim_dict = dict(deepcopy(self).__dict__)

        # convert acados classes to dicts
        for key, v in sim_dict.items():
            # skip non dict attributes
            if isinstance(v, (AcadosSim, AcadosSimDims, AcadosSimOptions, AcadosModel)):
                sim_dict[key]=dict(getattr(self, key).__dict__)

        return format_class_dict(sim_dict)


    def dump_to_json(self, json_file='acados_sim.json') -> None:
        with open(json_file, 'w') as f:
            json.dump(self.to_dict(), f, default=make_object_json_dumpable, indent=4, sort_keys=True)


    def render_templates(self, json_file, cmake_options: CMakeBuilder = None):
        # setting up loader and environment
        json_path = os.path.join(os.getcwd(), json_file)
        name = self.model.name

        if not os.path.exists(json_path):
            raise FileNotFoundError(f"{json_path} not found!")

        template_list = [
            ('acados_sim_solver.in.c', f'acados_sim_solver_{name}.c'),
            ('acados_sim_solver.in.h', f'acados_sim_solver_{name}.h'),
            ('acados_sim_solver.in.pxd', 'acados_sim_solver.pxd'),
            ('main_sim.in.c', f'main_sim_{name}.c'),
        ]
        if self.simulink_opts is not None:
            template_file = os.path.join('matlab_templates', 'mex_sim_solver.in.m')
            template_list.append((template_file, f'{name}_mex_sim_solver.m'))
            template_file = os.path.join('matlab_templates', 'make_mex_sim.in.m')
            template_list.append((template_file, f'make_mex_sim_{name}.m'))
            template_file = os.path.join('matlab_templates', 'acados_sim_create.in.c')
            template_list.append((template_file, f'acados_sim_create_{name}.c'))
            template_file = os.path.join('matlab_templates', 'acados_sim_free.in.c')
            template_list.append((template_file, f'acados_sim_free_{name}.c'))
            template_file = os.path.join('matlab_templates', 'acados_sim_set.in.c')
            template_list.append((template_file, f'acados_sim_set_{name}.c'))
            template_file = os.path.join('matlab_templates', 'acados_sim_solver_sfun.in.c')
            template_list.append((template_file, f'acados_sim_solver_sfunction_{name}.c'))
            template_file = os.path.join('matlab_templates', 'make_sfun_sim.in.m')
            template_list.append((template_file, f'make_sfun_sim_{name}.m'))

        # Builder
        if cmake_options is not None:
            template_list.append(('CMakeLists.in.txt', 'CMakeLists.txt'))
        else:
            template_list.append(('Makefile.in', 'Makefile'))

        # Render templates
        for (in_file, out_file) in template_list:
            render_template(in_file, out_file, self.code_export_directory, json_path)

        # folder model
        model_dir = os.path.join(self.code_export_directory, self.model.name + '_model')

        in_file = 'model.in.h'
        out_file = f'{self.model.name}_model.h'
        render_template(in_file, out_file, model_dir, json_path)


    def generate_external_functions(self, ):

        integrator_type = self.solver_options.integrator_type
        code_export_dir = self.code_export_directory

        opts = AcadosCodegenOptions(generate_hess = self.solver_options.sens_hess,
                    code_export_directory = self.code_export_directory,
                    ext_fun_expand_dyn = self.solver_options.ext_fun_expand_dyn,
                    ext_fun_expand_cost = False,
                    ext_fun_expand_constr = False,
                    ext_fun_expand_precompute = False,
                    )

        # create code_export_dir, model_dir
        model_dir = os.path.join(code_export_dir, self.model.name + '_model')
        if not os.path.exists(model_dir):
            os.makedirs(model_dir)

        context = GenerateContext(self.model.p_global, self.model.name, opts)

        # generate external functions
        check_casadi_version()
        if integrator_type == 'ERK':
            generate_c_code_explicit_ode(context, self.model, model_dir)
        elif integrator_type == 'IRK':
            generate_c_code_implicit_ode(context, self.model, model_dir)
        elif integrator_type == 'GNSF':
            generate_c_code_gnsf(context, self.model, model_dir)
        else:
            raise ValueError('Invalid integrator_type value. Possible values are:\n\n' \
                    + ',\n'.join(['ERK', 'IRK', 'GNSF']) + '.\n\nYou have: ' + integrator_type + '.\n\n')

        context.finalize()
        self.__external_function_files_model = context.get_external_function_file_list(ocp_specific=False)
