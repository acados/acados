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
import sys
import warnings
import time
from deprecated.sphinx import deprecated
from ctypes import (POINTER, byref, c_bool, c_char_p, c_double, c_int,
                    c_void_p, cast)
if os.name == 'nt':
    from ctypes import wintypes
    from ctypes import WinDLL as DllLoader
else:
    from ctypes import CDLL as DllLoader

import numpy as np

from .acados_ocp import AcadosOcp
from .acados_sim import AcadosSim

from .builders import CMakeBuilder
from .gnsf import detect_gnsf_structure
from .utils import (get_shared_lib_ext, get_shared_lib_prefix, get_shared_lib_dir, hash_class_instance,
                    status_to_str,
                    verbose_system_call, acados_lib_is_compiled_with_openmp,
                    get_shared_lib, set_directory)


class AcadosSimSolver:
    """
    Class to interact with the acados integrator C object.

    :param sim: type :py:class:`~acados_template.acados_ocp.AcadosOcp` (takes values to generate an instance :py:class:`~acados_template.sim.AcadosSim`) or :py:class:`~acados_template.sim.AcadosSim`
    :param json_file: Default: 'sim.json'
    :param build: Default: True
    :param cmake_builder: type :py:class:`~acados_template.utils.CMakeBuilder` generate a `CMakeLists.txt` and use
        the `CMake` pipeline instead of a `Makefile` (`CMake` seems to be the better option in conjunction with
        `MS Visual Studio`); default: `None`
    """
    if sys.platform=="win32":
        dlclose = DllLoader('kernel32', use_last_error=True).FreeLibrary
        dlclose.argtypes = [wintypes.HMODULE]
        winmode = 8 # why 8? what does that mean?
    else:
        dlclose = DllLoader(None).dlclose
        dlclose.argtypes = [c_void_p]
        winmode = None

    @property
    def acados_lib_uses_omp(self,):
        """``acados_lib_uses_omp`` - flag indicating whether the acados library has been compiled with openMP."""
        return self.__acados_lib_uses_omp

    @property
    def T(self,):
        """``T`` - Simulation time."""
        return self.__T

    @property
    def generated(self) -> bool:
        """Indicates whether code was generated or reused."""
        return self.__generated

    @staticmethod
    def generate(sim: AcadosSim, cmake_builder: CMakeBuilder = None, verbose: bool = True):
        """
        Generates the code for an acados sim solver, given the description in sim
        """

        sim.make_consistent()

        # module dependent post processing
        if sim.solver_options.integrator_type == 'GNSF':
            if sim.solver_options.sens_hess == True:
                raise ValueError("AcadosSimSolver: GNSF does not support sens_hess = True.")
            if 'gnsf_model' in sim.__dict__:
                raise ValueError("AcadosSim should not have gnsf_model, loading GNSF model functions from json is deprecated.")
            elif sim.model.gnsf_model is not None:
                # user provided GNSF model
                pass
            else:
                detect_gnsf_structure(sim.model, sim.dims)

        # generate code for external functions
        t0 = time.time()
        sim.generate_external_functions()
        t1 = time.time()
        if verbose:
            print(f"External functions generated in {1000*(t1-t0):.3f} ms.")

        sim.dump_to_json()

        t0 = time.time()
        sim.render_templates(cmake_builder)
        t1 = time.time()

        if verbose:
            print(f"Templated solver code generated in {1000*(t1-t0):.3f} ms.")


    @staticmethod
    def build(code_export_dir, with_cython=False, cmake_builder: CMakeBuilder = None, verbose: bool = True):

        t0 = time.time()
        with set_directory(code_export_dir):
            if with_cython:
                verbose_system_call(['make', 'clean_sim_cython'], verbose)
                verbose_system_call(['make', 'sim_cython'], verbose)
            else:
                if cmake_builder is not None:
                    cmake_builder.exec(code_export_dir, verbose)
                else:
                    verbose_system_call(['make', 'sim_shared_lib'], verbose)
        t1 = time.time()

        if verbose:
            print(f"Build completed in {1000*(t1-t0):.3f} ms.")


    @staticmethod
    def create_cython_solver(sim: AcadosSim, generate=True, build=True, cmake_builder: CMakeBuilder = None, verbose: bool = True, check_reuse_possible=True):

        sim.make_consistent()

        if generate is False and check_reuse_possible:
            reuse_possible = AcadosSimSolver.is_code_reuse_possible(sim, verbose=verbose)
            if not reuse_possible:
                generate = True
                build = True
                if verbose:
                    print("Code reuse not possible! Setting generate and build to True.")
            elif verbose:
                print("Code reuse possible, skipping code generation.")

        if generate:
            AcadosSimSolver.generate(sim, cmake_builder=cmake_builder, verbose=verbose)

        if build:
            AcadosSimSolver.build(sim.code_gen_options.code_export_directory, with_cython=True, cmake_builder=cmake_builder, verbose=verbose)


        importlib.invalidate_caches()
        sys.path.append(os.path.dirname(sim.code_gen_options.code_export_directory))

        sim_solver_pyx = importlib.import_module(f'{os.path.split(sim.code_gen_options.code_export_directory)[1]}.acados_sim_solver_pyx')

        AcadosSimSolverCython = getattr(sim_solver_pyx, 'AcadosSimSolverCython')
        return AcadosSimSolverCython(sim.name)


    def __init__(self, sim: AcadosSim, json_file=None, generate=True, build=True, cmake_builder: CMakeBuilder = None, verbose: bool = True, check_reuse_possible=True):

        self.solver_created = False

        # reuse existing json and casadi functions, when creating integrator from ocp
        if isinstance(sim, AcadosOcp):
            warnings.warn("An AcadosSimSolver is created from an AcadosOcp description, which is deprecated in 0.5.6. This only works if you created an AcadosOcpSolver before with the same description. Otherwise it leads to undefined behavior. Use AcadosSim.from_ocp() to first obtain a sim obj.",
                          DeprecationWarning,
                          stacklevel=2)
            generate = False

            if sim.dims.np_global > 0:
                raise ValueError("AcadosSimSolver: AcadosOcp with np_global > 0 is not supported.")

            self.__T = sim.solver_options.Tsim
        else:
            # formulation provided
            if json_file is not None:
                warnings.warn(" The `json_file` argument is deprecated in v0.5.6 and will be removed in a future release. Set AcadosSim.code_gen_options.json_file instead.", DeprecationWarning, stacklevel=2)
                sim.code_gen_options.json_file = json_file
            self.__T = sim.solver_options.T
            sim.make_consistent()

        if isinstance(sim, AcadosSim) and generate is False and check_reuse_possible:
            reuse_possible = AcadosSimSolver.is_code_reuse_possible(sim, verbose=verbose)
            if not reuse_possible:
                generate = True
                build = True
                if verbose:
                    print("Code reuse not possible! Setting generate and build to True.")
            elif verbose:
                print("Code reuse possible, skipping code generation.")

        if generate:
            self.generate(sim, cmake_builder=cmake_builder, verbose=verbose)
            self.__generated = True
        else:
            self.__generated = False

        if build:
            self.build(sim.code_gen_options.code_export_directory, cmake_builder=cmake_builder, verbose=verbose)

        self.__sim = sim

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
        libacados_filepath = os.path.join(sim.code_gen_options.acados_lib_path, '..', lib_dir, libacados_name)
        self.__acados_lib = get_shared_lib(libacados_filepath, self.winmode)

        # find out if acados was compiled with OpenMP
        self.__acados_lib_uses_omp = acados_lib_is_compiled_with_openmp(self.__acados_lib, verbose)

        libacados_sim_solver_name = f'{lib_prefix}acados_sim_solver_{self.sim.name}{lib_ext}'
        self.shared_lib_name = os.path.join(sim.code_export_directory, libacados_sim_solver_name)
        self.shared_lib = get_shared_lib(self.shared_lib_name, winmode=self.winmode)

        # create capsule
        getattr(self.shared_lib, f"{self.sim.name}_acados_sim_solver_create_capsule").restype = c_void_p
        self.capsule = getattr(self.shared_lib, f"{self.sim.name}_acados_sim_solver_create_capsule")()

        # create solver
        getattr(self.shared_lib, f"{self.sim.name}_acados_sim_create").argtypes = [c_void_p]
        getattr(self.shared_lib, f"{self.sim.name}_acados_sim_create").restype = c_int
        assert getattr(self.shared_lib, f"{self.sim.name}_acados_sim_create")(self.capsule)==0
        self.solver_created = True

        getattr(self.shared_lib, f"{self.sim.name}_acados_get_sim_opts").argtypes = [c_void_p]
        getattr(self.shared_lib, f"{self.sim.name}_acados_get_sim_opts").restype = c_void_p
        self.sim_opts = getattr(self.shared_lib, f"{self.sim.name}_acados_get_sim_opts")(self.capsule)

        getattr(self.shared_lib, f"{self.sim.name}_acados_get_sim_dims").argtypes = [c_void_p]
        getattr(self.shared_lib, f"{self.sim.name}_acados_get_sim_dims").restype = c_void_p
        self.sim_dims = getattr(self.shared_lib, f"{self.sim.name}_acados_get_sim_dims")(self.capsule)

        getattr(self.shared_lib, f"{self.sim.name}_acados_get_sim_config").argtypes = [c_void_p]
        getattr(self.shared_lib, f"{self.sim.name}_acados_get_sim_config").restype = c_void_p
        self.sim_config = getattr(self.shared_lib, f"{self.sim.name}_acados_get_sim_config")(self.capsule)

        getattr(self.shared_lib, f"{self.sim.name}_acados_get_sim_out").argtypes = [c_void_p]
        getattr(self.shared_lib, f"{self.sim.name}_acados_get_sim_out").restype = c_void_p
        self.sim_out = getattr(self.shared_lib, f"{self.sim.name}_acados_get_sim_out")(self.capsule)

        getattr(self.shared_lib, f"{self.sim.name}_acados_get_sim_in").argtypes = [c_void_p]
        getattr(self.shared_lib, f"{self.sim.name}_acados_get_sim_in").restype = c_void_p
        self.sim_in = getattr(self.shared_lib, f"{self.sim.name}_acados_get_sim_in")(self.capsule)

        getattr(self.shared_lib, f"{self.sim.name}_acados_get_sim_solver").argtypes = [c_void_p]
        getattr(self.shared_lib, f"{self.sim.name}_acados_get_sim_solver").restype = c_void_p
        self.sim_solver = getattr(self.shared_lib, f"{self.sim.name}_acados_get_sim_solver")(self.capsule)

        getattr(self.shared_lib, f"{self.sim.name}_acados_get_sim_mem").argtypes = [c_void_p]
        getattr(self.shared_lib, f"{self.sim.name}_acados_get_sim_mem").restype = c_void_p
        self.sim_mem = getattr(self.shared_lib, f"{self.sim.name}_acados_get_sim_mem")(self.capsule)


        # argtypes and restypes
        self.__acados_lib.sim_out_get.argtypes = [c_void_p, c_void_p, c_void_p, c_char_p, c_void_p]
        self.__acados_lib.sim_dims_get_from_attr.argtypes = [c_void_p, c_void_p, c_char_p, POINTER(c_int)]

        self.__acados_lib.sim_memory_get.argtypes = [c_void_p, c_void_p, c_void_p, c_char_p, c_void_p]
        self.__acados_lib.sim_memory_get.restype = None

        self.__acados_lib.sim_solver_set.argtypes = [c_void_p, c_char_p, c_void_p]
        self.__acados_lib.sim_in_set.argtypes = [c_void_p, c_void_p, c_void_p, c_char_p, c_void_p]
        self.__acados_lib.sim_opts_set.argtypes = [c_void_p, c_void_p, c_char_p, POINTER(c_bool)]

        getattr(self.shared_lib, f"{self.sim.name}_acados_sim_update_params").argtypes = [c_void_p, POINTER(c_double), c_int]

        getattr(self.shared_lib, f"{self.sim.name}_acados_sim_solve").argtypes = [c_void_p]
        getattr(self.shared_lib, f"{self.sim.name}_acados_sim_solve").restype = c_int

        self.gettable_vectors = ['x', 'u', 'z', 'S_adj']
        self.gettable_matrices = ['S_forw', 'Sx', 'Su', 'S_hess', 'S_algebraic', 'S_p']
        self.gettable_scalars = ['CPUtime', 'time_tot', 'ADtime', 'time_ad', 'LAtime', 'time_la']


    @property
    def sim(self):
        return self.__sim

    @property
    @deprecated(version="0.5.6", reason="The property acados_sim is deprecated. Use sim instead")
    def acados_sim(self):
        return self.__sim

    @staticmethod
    def is_code_reuse_possible(sim: AcadosSim, verbose: bool) -> bool:
        try:
            # Check if code_export_dir exists
            if not os.path.exists(sim.code_gen_options.code_export_directory):
                return False

            # Check if JSON file exists
            if not os.path.exists(sim.code_gen_options.json_file):
                return False

            # Load existing JSON and extract hash
            with open(sim.code_gen_options.json_file, 'r') as f:
                existing_data = json.load(f)

            if 'hash' not in existing_data:
                return False

            existing_hash = existing_data['hash']

            # Create hash of current Sim
            current_hash = hash_class_instance(sim)

            # Compare hashes
            reuse_possible = current_hash == existing_hash
            if not reuse_possible and verbose:
                print("Sim formulation has changed, code reuse not possible.")
            return reuse_possible

        except Exception:
            # If any error occurs during comparison, return False to trigger regeneration
            return False

    def simulate(self, x=None, u=None, z=None, xdot=None, p=None):
        """
        Simulate the system forward for the given x, u, p and return x_next.
        The values xdot, z are used as initial guesses for implicit integrators, if provided.
        Wrapper around `solve()` taking care of setting/getting inputs/outputs.
        """
        if x is not None:
            self.set('x', x)
        if u is not None:
            self.set('u', u)
        if self.sim.solver_options.integrator_type == "IRK":
            if z is not None:
                self.set('z', z)
            if xdot is not None:
                self.set('xdot', xdot)
        if p is not None:
            self.set('p', p)

        status = self.solve()

        if status != 0:
            raise RuntimeError(f'AcadosSimSolver {self.sim.name} returned status {status} ({status_to_str(status)}).')

        x_next = self.get('x')
        return x_next


    def solve(self):
        """
        Solve the simulation problem with current input.
        """
        status = getattr(self.shared_lib, f"{self.sim.name}_acados_sim_solve")(self.capsule)
        return status


    def get(self, field_):
        """
        Get the last solution of the solver.

        :param field: string in ['x', 'u', 'z', 'S_forw', 'Sx', 'Su', 'S_adj', 'S_hess', 'S_algebraic', 'CPUtime', 'time_tot', 'ADtime', 'time_ad', 'LAtime', 'time_la']
        """
        field = field_.encode('utf-8')

        if field_ in self.gettable_vectors:
            # get dims
            dims = np.zeros((1,), dtype=np.intc, order='F')
            dims_data = cast(dims.ctypes.data, POINTER(c_int))

            self.__acados_lib.sim_dims_get_from_attr(self.sim_config, self.sim_dims, field, dims_data)

            # allocate array
            out = np.zeros((dims[0],), dtype=np.float64, order='F')
            out_data = cast(out.ctypes.data, POINTER(c_double))

            self.__acados_lib.sim_out_get(self.sim_config, self.sim_dims, self.sim_out, field, out_data)

        elif field_ in self.gettable_matrices:
            # get dims
            dims = np.zeros((2,), dtype=np.intc, order='F')
            dims_data = cast(dims.ctypes.data, POINTER(c_int))

            self.__acados_lib.sim_dims_get_from_attr(self.sim_config, self.sim_dims, field, dims_data)

            out = np.zeros((dims[0], dims[1]), dtype=np.float64, order='F')
            out_data = cast(out.ctypes.data, POINTER(c_double))

            # S_p is stored only in integrator memory (not in sim_out)
            if field_ == 'S_p':
                self.__acados_lib.sim_memory_get(self.sim_config, self.sim_dims, self.sim_mem, field, out_data)
            else:
                self.__acados_lib.sim_out_get(self.sim_config, self.sim_dims, self.sim_out, field, out_data)


        elif field_ in self.gettable_scalars:
            scalar = c_double()
            scalar_data = byref(scalar)
            self.__acados_lib.sim_out_get(self.sim_config, self.sim_dims, self.sim_out, field, scalar_data)

            out = scalar.value
        else:
            raise KeyError(f'AcadosSimSolver.get(): Unknown field {field_},'
                f' available fields are {", ".join(self.gettable_vectors+self.gettable_matrices)}, {", ".join(self.gettable_scalars)}')

        return out



    def set(self, field_: str, value_):
        """
        Set numerical data inside the solver.

        :param field: string in ['x', 'u', 'p', 'xdot', 'z', 'seed_adj', 'T', 't0']
        :param value: the value with appropriate size.
        """
        settable = ['x', 'u', 'p', 'xdot', 'z', 'seed_adj', 'T', 't0'] # S_forw

        # TODO: check and throw error here. then remove checks in Cython for speed
        # cast value_ to avoid conversion issues
        if isinstance(value_, (float, int)):
            value_ = np.array([value_])

        value_ = value_.astype(float)
        value_data = cast(value_.ctypes.data, POINTER(c_double))
        value_data_p = cast((value_data), c_void_p)

        field = field_.encode('utf-8')

        # treat parameters separately
        if field_ == 'p':
            value_data = cast(value_.ctypes.data, POINTER(c_double))
            getattr(self.shared_lib, f"{self.sim.name}_acados_sim_update_params")(self.capsule, value_data, value_.shape[0])
            return
        else:
            # dimension check
            dims = np.zeros((2,), dtype=np.intc, order='F')
            dims_data = cast(dims.ctypes.data, POINTER(c_int))

            self.__acados_lib.sim_dims_get_from_attr(self.sim_config, self.sim_dims, field, dims_data)

            # TODO: isn't the shape check afterwards meaningless if we ravel first?
            value_ = np.ravel(value_, order='F')

            value_shape = value_.shape
            if len(value_shape) == 1:
                value_shape = (value_shape[0], 0)

            if value_shape != tuple(dims):
                raise ValueError(f'AcadosSimSolver.set(): mismatching dimension' \
                    f' for field "{field_}" with dimension {tuple(dims)} (you have {value_shape}).')

            if field_ == 'T':
                self.__T = value_

        # set
        if field_ in ['xdot', 'z']:
            self.__acados_lib.sim_solver_set(self.sim_solver, field, value_data_p)
        elif field_ in settable:
            self.__acados_lib.sim_in_set(self.sim_config, self.sim_dims, self.sim_in, field, value_data_p)
        else:
            raise KeyError(f'AcadosSimSolver.set(): Unknown field {field_},'
                f' available fields are {", ".join(settable)}')

        return


    def options_set(self, field_: str, value_: bool):
        """
        Set solver options

        :param field: string in ['sens_forw', 'sens_adj', 'sens_hess']
        :param value: Boolean
        """
        fields = ['sens_forw', 'sens_adj', 'sens_hess']
        if field_ not in fields:
            raise ValueError(f"field {field_} not supported. Supported values are {', '.join(fields)}.\n")

        field = field_.encode('utf-8')
        value_ctypes = c_bool(value_)

        if not isinstance(value_, bool):
            raise TypeError("options_set: expected boolean for value")

        # only allow setting
        if getattr(self.sim.solver_options, field_) or value_ == False:
            self.__acados_lib.sim_opts_set(self.sim_config, self.sim_opts, field, value_ctypes)
        else:
            raise RuntimeError(f"Cannot set option {field_} to True, because it was False in original solver options.\n")

        return


    def __del__(self):

        if self.solver_created:
            getattr(self.shared_lib, f"{self.sim.name}_acados_sim_free").argtypes = [c_void_p]
            getattr(self.shared_lib, f"{self.sim.name}_acados_sim_free").restype = c_int
            getattr(self.shared_lib, f"{self.sim.name}_acados_sim_free")(self.capsule)

            getattr(self.shared_lib, f"{self.sim.name}_acados_sim_solver_free_capsule").argtypes = [c_void_p]
            getattr(self.shared_lib, f"{self.sim.name}_acados_sim_solver_free_capsule").restype = c_int
            getattr(self.shared_lib, f"{self.sim.name}_acados_sim_solver_free_capsule")(self.capsule)

            try:
                self.dlclose(self.shared_lib._handle)
            except:
                warnings.warn(f"acados Python interface could not close shared_lib handle of AcadosSimSolver {self.sim.name}.\n"
                     "Attempting to create a new one with the same name will likely result in the old one being used!")
                pass
