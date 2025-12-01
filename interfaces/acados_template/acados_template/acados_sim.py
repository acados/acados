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
from typing import Optional
from copy import deepcopy
from deprecated.sphinx import deprecated
import warnings
import casadi as ca
from .acados_model import AcadosModel
from .acados_dims import AcadosSimDims
from .acados_ocp import AcadosOcp
from .builders import CMakeBuilder
from .ros2.sim_node import AcadosSimRosOptions
from .utils import (get_acados_path, get_shared_lib_ext, format_class_dict, check_casadi_version,
                    make_object_json_dumpable, render_template, is_scalar_integer, is_empty)
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

        # NOTE: internal names have sim_method_ prefix for intercompatibility of sim and ocp templates
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

    @integrator_type.setter
    def integrator_type(self, integrator_type):
        integrator_types = ('ERK', 'IRK', 'GNSF')
        if integrator_type in integrator_types:
            self.__integrator_type = integrator_type
        else:
            raise ValueError('Invalid integrator_type value. Possible values are:\n\n' \
                    + ',\n'.join(integrator_types) + '.\n\nYou have: ' + integrator_type + '.\n\n')

    @property
    def num_stages(self):
        """Number of stages in the integrator. Default: 4"""
        return self.__sim_method_num_stages

    @num_stages.setter
    def num_stages(self, num_stages):
        if is_scalar_integer(num_stages):
            self.__sim_method_num_stages = num_stages
        else:
            raise ValueError('Invalid num_stages value. num_stages must be an integer.')

    @property
    def num_steps(self):
        """Number of steps in the integrator. Default: 1"""
        return self.__sim_method_num_steps

    @num_steps.setter
    def num_steps(self, num_steps):
        if is_scalar_integer(num_steps):
            self.__sim_method_num_steps = num_steps
        else:
            raise TypeError('Invalid num_steps value. num_steps must be an integer.')

    @property
    def newton_iter(self):
        """Number of Newton iterations in simulation method. Default: 3"""
        return self.__sim_method_newton_iter

    @newton_iter.setter
    def newton_iter(self, newton_iter):
        if is_scalar_integer(newton_iter):
            self.__sim_method_newton_iter = newton_iter
        else:
            raise TypeError('Invalid newton_iter value. newton_iter must be an integer.')

    @property
    def newton_tol(self):
        """
        Tolerance for Newton system solved in implicit integrator (IRK, GNSF).
        0.0 means this is not used and exactly newton_iter iterations are carried out.
        Default: 0.0
        """
        return self.__sim_method_newton_tol

    @newton_tol.setter
    def newton_tol(self, newton_tol):
        if isinstance(newton_tol, float):
            self.__sim_method_newton_tol = newton_tol
        else:
            raise TypeError('Invalid newton_tol value. newton_tol must be a float.')

    @property
    def sens_forw(self):
        """Boolean determining if forward sensitivities are computed. Default: True"""
        return self.__sens_forw

    @sens_forw.setter
    def sens_forw(self, sens_forw):
        if isinstance(sens_forw, bool):
            self.__sens_forw = sens_forw
        else:
            raise ValueError('Invalid sens_forw value. sens_forw must be a Boolean.')

    @property
    def sens_adj(self):
        """Boolean determining if adjoint sensitivities are computed. Default: False"""
        return self.__sens_adj

    @sens_adj.setter
    def sens_adj(self, sens_adj):
        if isinstance(sens_adj, bool):
            self.__sens_adj = sens_adj
        else:
            raise ValueError('Invalid sens_adj value. sens_adj must be a Boolean.')

    @property
    def sens_algebraic(self):
        """Boolean determining if sensitivities wrt algebraic variables are computed. Default: False"""
        return self.__sens_algebraic

    @sens_algebraic.setter
    def sens_algebraic(self, sens_algebraic):
        if isinstance(sens_algebraic, bool):
            self.__sens_algebraic = sens_algebraic
        else:
            raise ValueError('Invalid sens_algebraic value. sens_algebraic must be a Boolean.')

    @property
    def sens_hess(self):
        """Boolean determining if hessians are computed. Default: False"""
        return self.__sens_hess

    @sens_hess.setter
    def sens_hess(self, sens_hess):
        if isinstance(sens_hess, bool):
            self.__sens_hess = sens_hess
        else:
            raise ValueError('Invalid sens_hess value. sens_hess must be a Boolean.')

    @property
    def output_z(self):
        """Boolean determining if values for algebraic variables (corresponding to start of simulation interval) are computed. Default: True"""
        return self.__output_z

    @output_z.setter
    def output_z(self, output_z):
        if isinstance(output_z, bool):
            self.__output_z = output_z
        else:
            raise ValueError('Invalid output_z value. output_z must be a Boolean.')

    @property
    # TODO: rename to jac_reuse
    def sim_method_jac_reuse(self):
        """Integer determining if jacobians are reused (0 or 1). Default: 0"""
        return self.__sim_method_jac_reuse

    @sim_method_jac_reuse.setter
    def sim_method_jac_reuse(self, sim_method_jac_reuse):
        # TODO: use bool
        if sim_method_jac_reuse in (0, 1):
            self.__sim_method_jac_reuse = sim_method_jac_reuse
        else:
            raise ValueError('Invalid sim_method_jac_reuse value. sim_method_jac_reuse must be 0 or 1.')

    @property
    def T(self):
        """Time horizon"""
        return self.__Tsim

    @T.setter
    def T(self, T):
        self.__Tsim = T

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

    @collocation_type.setter
    def collocation_type(self, collocation_type):
        collocation_types = ('GAUSS_RADAU_IIA', 'GAUSS_LEGENDRE')
        if collocation_type in collocation_types:
            self.__collocation_type = collocation_type
        else:
            raise ValueError('Invalid collocation_type value. Possible values are:\n\n' \
                    + ',\n'.join(collocation_types) + '.\n\nYou have: ' + collocation_type + '.\n\n')

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
    def ext_fun_expand_dyn(self):
        """
        Flag indicating whether CasADi.MX should be expanded to CasADi.SX before code generation.
        Default: False
        """
        return self.__ext_fun_expand_dyn


    @ext_fun_expand_dyn.setter
    def ext_fun_expand_dyn(self, ext_fun_expand_dyn):
        if isinstance(ext_fun_expand_dyn, bool):
            self.__ext_fun_expand_dyn = ext_fun_expand_dyn
        else:
            raise TypeError('Invalid ext_fun_expand_dyn value, expected bool.\n')

    @property
    @deprecated(version="0.4.0", reason="Set the flag with_batch_functionality instead and pass the number of threads directly to the BatchSolver.")
    def num_threads_in_batch_solve(self):
        """
        Integer indicating how many threads should be used within the batch solve.
        If more than one thread should be used, the sim solver is compiled with openmp.
        Default: 1.
        """
        return self.__num_threads_in_batch_solve

    @num_threads_in_batch_solve.setter
    @deprecated(version="0.4.0", reason="Set the flag with_batch_functionality instead and pass the number of threads directly to the BatchSolver.")
    def num_threads_in_batch_solve(self, num_threads_in_batch_solve):
        if isinstance(num_threads_in_batch_solve, int) and num_threads_in_batch_solve > 0:
            self.__num_threads_in_batch_solve = num_threads_in_batch_solve
        else:
            raise ValueError('Invalid num_threads_in_batch_solve value. num_threads_in_batch_solve must be a positive integer.')

    @property
    def with_batch_functionality(self):
        """
        Whether the AcadosSimBatchSolver can be used.
        In this case, the sim solver is compiled with openmp.
        Default: False.
        """
        return self.__with_batch_functionality

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
        self.__json_file = "acados_sim.json"

        self.__ros_opts: Optional[AcadosSimRosOptions] = None

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

    @property
    def json_file(self):
        """Name of the json file where the problem description is stored."""
        return self.__json_file

    @json_file.setter
    def json_file(self, json_file):
        self.__json_file = json_file

    @property
    def ros_opts(self) -> Optional[AcadosSimRosOptions]:
        """Options to configure ROS 2 nodes and topics."""
        return self.__ros_opts

    @ros_opts.setter
    def ros_opts(self, ros_opts: AcadosSimRosOptions):
        if not isinstance(ros_opts, AcadosSimRosOptions):
            raise TypeError('Invalid ros_opts value, expected AcadosOcpRos.\n')
        self.__ros_opts = ros_opts

    def make_consistent(self):
        self.model.make_consistent(self.dims)
        self.name = self.model.name

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
            if isinstance(v, (AcadosSim, AcadosSimDims, AcadosSimOptions)):
                sim_dict[key] = dict(getattr(self, key).__dict__)
            if isinstance(v, (AcadosModel, AcadosSimRosOptions)):
                sim_dict[key] = v.to_dict()

        return format_class_dict(sim_dict)


    def dump_to_json(self) -> None:
        dir_name = os.path.dirname(self.json_file)
        if dir_name:
            os.makedirs(dir_name, exist_ok=True)

        with open(self.json_file, 'w') as f:
            json.dump(self.to_dict(), f, default=make_object_json_dumpable, indent=4, sort_keys=True)


    def _get_ros_template_list(self) -> list:
        template_list = []
        acados_template_path = os.path.dirname(os.path.abspath(__file__))
        ros_template_glob = os.path.join(acados_template_path, 'ros2_templates', '**', '*')

        # --- Interface Package ---
        ros_interface_dir = os.path.join('sim_interface_templates')
        interface_dir = os.path.join(self.ros_opts.generated_code_dir, f'{self.ros_opts.package_name}_interface')
        template_file = os.path.join(ros_interface_dir, 'README.in.md')
        template_list.append((template_file, 'README.md', interface_dir, ros_template_glob))
        template_file = os.path.join(ros_interface_dir, 'CMakeLists.in.txt')
        template_list.append((template_file, 'CMakeLists.txt', interface_dir, ros_template_glob))
        template_file = os.path.join(ros_interface_dir, 'package.in.xml')
        template_list.append((template_file, 'package.xml', interface_dir, ros_template_glob))

        # Messages
        msg_dir = os.path.join(interface_dir, 'msg')
        template_file = os.path.join(ros_interface_dir, 'State.in.msg')
        template_list.append((template_file, 'State.msg', msg_dir, ros_template_glob))
        template_file = os.path.join(ros_interface_dir, 'ControlInput.in.msg')
        template_list.append((template_file, 'ControlInput.msg', msg_dir, ros_template_glob))

        # --- Simulator Package ---
        ros_pkg_dir = os.path.join('sim_node_templates')
        package_dir = os.path.join(self.ros_opts.generated_code_dir, self.ros_opts.package_name)
        template_file = os.path.join(ros_pkg_dir, 'README.in.md')
        template_list.append((template_file, 'README.md', package_dir, ros_template_glob))
        template_file = os.path.join(ros_pkg_dir, 'CMakeLists.in.txt')
        template_list.append((template_file, 'CMakeLists.txt', package_dir, ros_template_glob))
        template_file = os.path.join(ros_pkg_dir, 'package.in.xml')
        template_list.append((template_file, 'package.xml', package_dir, ros_template_glob))

        # Header
        include_dir = os.path.join(package_dir, 'include', self.ros_opts.package_name)
        template_file = os.path.join(ros_pkg_dir, 'config.in.hpp')
        template_list.append((template_file, 'config.hpp', include_dir, ros_template_glob))
        template_file = os.path.join(ros_pkg_dir, 'utils.in.hpp')
        template_list.append((template_file, 'utils.hpp', include_dir, ros_template_glob))
        template_file = os.path.join(ros_pkg_dir, 'node.in.h')
        template_list.append((template_file, 'node.h', include_dir, ros_template_glob))

        # Source
        src_dir = os.path.join(package_dir, 'src')
        template_file = os.path.join(ros_pkg_dir, 'node.in.cpp')
        template_list.append((template_file, 'node.cpp', src_dir, ros_template_glob))

        # Hooks
        hooks_dir = os.path.join(package_dir, 'env-hooks')
        template_file = os.path.join(ros_pkg_dir, 'export_acados_path.sh.in')
        # Note: still a ".in" file, because it needs to be rendered by colcon build
        template_list.append((template_file, 'export_acados_path.sh.in', hooks_dir, ros_template_glob))

        # Test
        test_dir = os.path.join(package_dir, 'test')
        template_file = os.path.join(ros_pkg_dir, 'test.launch.in.py')
        template_list.append((template_file, f'test_{self.ros_opts.package_name}.launch.py', test_dir, ros_template_glob))
        return template_list


    def _get_simulink_template_list(self, name: str) -> list:
        template_list = []
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
        return template_list


    def render_templates(self, cmake_options: CMakeBuilder = None):
        # setting up loader and environment
        json_path = os.path.abspath(self.json_file)
        name = self.model.name

        if not os.path.exists(json_path):
            raise FileNotFoundError(f"{json_path} not found!")

        template_list = [
            ('acados_sim_solver.in.c', f'acados_sim_solver_{name}.c'),
            ('acados_sim_solver.in.h', f'acados_sim_solver_{name}.h'),
            ('acados_sim_solver.in.pxd', 'acados_sim_solver.pxd'),
            ('main_sim.in.c', f'main_sim_{name}.c'),
        ]

        # Model
        model_dir = os.path.join(self.code_export_directory, self.model.name + '_model')
        template_list.append(('model.in.h', f'{self.model.name}_model.h', model_dir))

        # Simulink
        if self.simulink_opts is not None:
            template_list += self._get_simulink_template_list(name)

        # ROS2
        if self.ros_opts is not None:
            template_list += self._get_ros_template_list()

        # Builder
        if cmake_options is not None:
            template_list.append(('CMakeLists.in.txt', 'CMakeLists.txt'))
        else:
            template_list.append(('Makefile.in', 'Makefile'))

        # Render templates
        for tup in template_list:
            output_dir = self.code_export_directory if len(tup) <= 2 else tup[2]
            template_glob = None if len(tup) <= 3 else tup[3]
            render_template(tup[0], tup[1], output_dir, json_path, template_glob=template_glob)


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


    @classmethod
    def from_ocp(cls, ocp: AcadosOcp):
        """
        Create an AcadosSim object from an AcadosOcp object.
        The AcadosSim object matches the integrator of the OCP on the first shooting node.

        :param ocp: AcadosOcp
        :return: AcadosSim
        """
        ocp.make_consistent()

        sim = cls()
        sim.model = deepcopy(ocp.model)

        if not is_empty(sim.model.p_global):
            sim.model.p = ca.vertcat(sim.model.p, sim.model.p_global)
            sim.model.p_global = []

            warnings.warn('Model contained p_global. Appending p_global to p in the sim model.')


        sim.solver_options.integrator_type = ocp.solver_options.integrator_type
        sim.solver_options.collocation_type = ocp.solver_options.collocation_type
        sim.solver_options.T = ocp.solver_options.time_steps[0]

        sim.solver_options.num_stages = ocp.solver_options.sim_method_num_stages[0]
        sim.solver_options.num_steps = ocp.solver_options.sim_method_num_steps[0]
        sim.solver_options.newton_iter = ocp.solver_options.sim_method_newton_iter
        sim.solver_options.newton_tol = ocp.solver_options.sim_method_newton_tol

        sim.solver_options.sim_method_jac_reuse = ocp.solver_options.sim_method_jac_reuse[0]

        return sim
