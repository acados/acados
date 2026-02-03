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
from typing import Optional, Union
import numpy as np

from ctypes import (POINTER, byref, c_char_p, c_double, c_int, c_bool,
                    c_void_p, cast)
if os.name == 'nt':
    from ctypes import wintypes
    from ctypes import WinDLL as DllLoader
else:
    from ctypes import CDLL as DllLoader

from .acados_ocp_qp import AcadosOcpQp
from .acados_ocp_options import AcadosOcpQpOptions

from .utils import (get_acados_path, get_shared_lib_ext, get_shared_lib_prefix, get_shared_lib_dir, get_shared_lib,
                    acados_lib_is_compiled_with_openmp)
from .acados_ocp_iterate import AcadosOcpIterate, AcadosOcpFlattenedIterate


class AcadosOcpQpSolver:
    if os.name == 'nt':
        dlclose = DllLoader('kernel32', use_last_error=True).FreeLibrary
        dlclose.argtypes = [wintypes.HMODULE]
        winmode = 8 # why 8? what does that mean?
    else:
        dlclose = DllLoader(None).dlclose
        dlclose.argtypes = [c_void_p]
        winmode = None

    @property
    def N(self) -> int:
        return self.__N

    @property
    def name(self) -> int:
        return self.__name

    def __init__(self, qp: AcadosOcpQp, opts: Optional[AcadosOcpQpOptions] = None, verbose: bool = False, acados_lib_path: str = None):

        self.__solver_created = False
        self.__N = qp.N
        self.qp = qp
        if opts is None:
            opts = AcadosOcpQpOptions()
        opts.make_consistent(qp.N)

        # prepare library loading
        lib_ext = get_shared_lib_ext()
        lib_prefix = get_shared_lib_prefix()
        lib_dir = get_shared_lib_dir()

        # Load acados library to avoid unloading the library.
        if acados_lib_path is None:
            acados_path = get_acados_path()
            acados_lib_path = os.path.join(acados_path, 'lib')

        libacados_name = f'{lib_prefix}acados{lib_ext}'
        libacados_filepath = os.path.join(acados_lib_path, '..', lib_dir, libacados_name)
        self.__acados_lib = get_shared_lib(libacados_filepath, self.winmode)

        # find out if acados was compiled with OpenMP
        self.__acados_lib_uses_omp = acados_lib_is_compiled_with_openmp(self.__acados_lib, verbose)

        self.__acados_lib.ocp_qp_xcond_solver_config_create_from_name.argtypes = [c_char_p]
        self.__acados_lib.ocp_qp_xcond_solver_config_create_from_name.restype = c_void_p
        self.c_config = self.__acados_lib.ocp_qp_xcond_solver_config_create_from_name(opts.qp_solver.encode('utf-8'))

        # Create dimensions structure
        self.__acados_lib.ocp_qp_xcond_solver_dims_create.argtypes = [c_void_p, c_int]
        self.__acados_lib.ocp_qp_xcond_solver_dims_create.restype = c_void_p
        self.c_dims = self.__acados_lib.ocp_qp_xcond_solver_dims_create(self.c_config, c_int(qp.N))

        # Set dimensions for each stage
        self._set_dimensions_in_c()

        # Opts
        # void *ocp_qp_xcond_solver_opts_create(ocp_qp_xcond_solver_config *config, ocp_qp_xcond_solver_dims *dims)
        self.__acados_lib.ocp_qp_xcond_solver_opts_create.argtypes = [c_void_p, c_void_p]
        self.__acados_lib.ocp_qp_xcond_solver_opts_create.restype = c_void_p
        self.c_opts = self.__acados_lib.ocp_qp_xcond_solver_opts_create(self.c_config, self.c_dims)

        # TODO: set opts!
        # void ocp_qp_xcond_solver_opts_set(ocp_qp_xcond_solver_config *config, ocp_qp_xcond_solver_opts *opts, const char *field, void* value)
        self.__acados_lib.ocp_qp_xcond_solver_opts_set.argtypes = [c_void_p, c_void_p, c_char_p, c_void_p]
        self.__acados_lib.ocp_qp_xcond_solver_opts_set.restype = None
        self._set_opts_from_class(opts)

        # ocp_qp_solver *ocp_qp_create(ocp_qp_xcond_solver_config *config, ocp_qp_xcond_solver_dims *dims, void *opts_)
        self.__acados_lib.ocp_qp_create.argtypes = [c_void_p, c_void_p, c_void_p]
        self.__acados_lib.ocp_qp_create.restype = c_void_p
        self.c_solver = self.__acados_lib.ocp_qp_create(self.c_config, self.c_dims, self.c_opts)

        # Create in and out
        self.__acados_lib.ocp_qp_in_create_from_xcond_dims.argtypes = [c_void_p]
        self.__acados_lib.ocp_qp_in_create_from_xcond_dims.restype = c_void_p
        self.c_in = self.__acados_lib.ocp_qp_in_create_from_xcond_dims(self.c_dims)

        self.__acados_lib.ocp_qp_out_create_from_xcond_dims.argtypes = [c_void_p]
        self.__acados_lib.ocp_qp_out_create_from_xcond_dims.restype = c_void_p
        self.c_out = self.__acados_lib.ocp_qp_out_create_from_xcond_dims(self.c_dims)

        # Set QP data in C structures
        self._set_initial_qp_data_in_c()

        ### initialize function prototypes
        # int ocp_qp_solve(ocp_qp_solver *solver, ocp_qp_in *qp_in, ocp_qp_out *qp_out)
        self.__acados_lib.ocp_qp_solve.argtypes = [c_void_p, c_void_p, c_void_p]
        self.__acados_lib.ocp_qp_solve.restype = c_int

        # void ocp_qp_out_get(ocp_qp_out *out, int stage, const char *field, void *value)
        self.__acados_lib.ocp_qp_out_get.argtypes = [c_void_p, c_int, c_char_p, c_void_p]
        self.__acados_lib.ocp_qp_out_get.restype = None

        # void ocp_qp_xcond_solver_get_scalar(ocp_qp_solver *solver, ocp_qp_out *qp_out, const char *field, void* value)
        self.__acados_lib.ocp_qp_xcond_solver_get_scalar.argtypes = [c_void_p, c_void_p, c_char_p, c_void_p]
        self.__acados_lib.ocp_qp_xcond_solver_get_scalar.restype = None

        self.__solver_created = True
        self._status = 0

        return

    def _set_dimensions_in_c(self):
        """
        Set the dimensions in the C solver from the QP object.
        """
        N = self.qp.N
        dims = self.qp.dims

        self.__acados_lib.ocp_qp_xcond_solver_dims_set.argtypes = [c_void_p, c_void_p, c_int, c_char_p, POINTER(c_int)]
        self.__acados_lib.ocp_qp_xcond_solver_dims_set.restype = None

        # Set dimensions for each stage (0 to N)
        for i in range(N + 1):
            # Number of states
            nx_data = (c_int)(dims.nx[i])
            self.__acados_lib.ocp_qp_xcond_solver_dims_set(self.c_config, self.c_dims, c_int(i), 'nx'.encode('utf-8'), byref(nx_data))

            # Number of inputs
            nu_data = (c_int)(dims.nu[i])
            self.__acados_lib.ocp_qp_xcond_solver_dims_set(self.c_config, self.c_dims, c_int(i), 'nu'.encode('utf-8'), nu_data)

            # Number of box constraints
            nbx_data = (c_int)(dims.nbx[i])
            self.__acados_lib.ocp_qp_xcond_solver_dims_set(self.c_config, self.c_dims, c_int(i), 'nbx'.encode('utf-8'), nbx_data)

            nbu_data = (c_int)(dims.nbu[i])
            self.__acados_lib.ocp_qp_xcond_solver_dims_set(self.c_config, self.c_dims, c_int(i), 'nbu'.encode('utf-8'), nbu_data)

            # Number of general constraints
            ng_data = (c_int)(dims.ng[i])
            self.__acados_lib.ocp_qp_xcond_solver_dims_set(self.c_config, self.c_dims, c_int(i), 'ng'.encode('utf-8'), ng_data)

            # Number of slack variables
            ns_data = (c_int)(dims.ns[i])
            self.__acados_lib.ocp_qp_xcond_solver_dims_set(self.c_config, self.c_dims, c_int(i), 'ns'.encode('utf-8'), ns_data)

        return

    def opts_set(self, field: str, value):
        """
        Set option in the C solver.

        :param field: string - name of the option to set.
        :param value: value of the option to set.
        """
        fields = ['tol_stat',
                'tol_eq',
                'tol_ineq',
                'tol_comp',
                'iter_max',
                'cond_N',
                'cond_block_size',
                'warm_start',
                'cond_ric_alg',
                'ric_alg',
                'mu0',
                't0_init',
                'print_level',
                'hpipm_mode'
                ]
        if field not in fields:
            raise ValueError(f'AcadosOcpQpSolver.opts_set(field={field}, value={value}): \'{field}\' is an invalid argument.'
                             f'\n Possible values are {fields}.')
        if self.__solver_created and field in ['cond_N', 'qp_solver', 'cond_block_size']:
            raise RuntimeError(f'AcadosOcpQpSolver.opts_set(field={field}, value={value}): cannot set option \'{field}\' after solver creation.')
        if field == 'cond_block_size':
            value = np.ascontiguousarray(value, dtype=np.intc)
            value_ptr = cast(value.ctypes.data, POINTER(c_int))
        elif isinstance(value, float):
            value_c = c_double(value)
            value_ptr = byref(value_c)
        elif isinstance(value, int):
            value_c = c_int(value)
            value_ptr = byref(value_c)
        elif isinstance(value, bool):
            value_c = c_bool(value)
            value_ptr = byref(value_c)
        elif isinstance(value, str):
            value_ptr = value.encode('utf-8')
        else:
            raise TypeError(f'AcadosOcpQpSolver.opts_set(field={field}, value={value}): unsupported type {type(value)} for value.')

        self.__acados_lib.ocp_qp_xcond_solver_opts_set(self.c_config, self.c_opts, field.encode('utf-8'), value_ptr)

    def _set_opts_from_class(self, opts: AcadosOcpQpOptions):
        if 'HPIPM' in opts.qp_solver:
            self.opts_set('hpipm_mode', opts.hpipm_mode)
            self.opts_set('t0_init', opts.t0_init)
            if opts.mu0 is not None:
                self.opts_set('mu0', opts.mu0)
        if opts.qp_solver == "PARTIAL_CONDENSING_HPIPM":
            self.opts_set('ric_alg', opts.ric_alg)
        # tols and iter
        self.opts_set('tol_stat', opts.tol_stat)
        self.opts_set('tol_eq', opts.tol_eq)
        self.opts_set('tol_ineq', opts.tol_ineq)
        self.opts_set('tol_comp', opts.tol_comp)
        self.opts_set('iter_max', opts.iter_max)
        # condensing
        if 'PARTIAL_CONDENSING' in opts.qp_solver:
            self.opts_set('cond_N', opts.cond_N)
        if opts.cond_block_size is not None:
            self.opts_set('cond_block_size', opts.cond_block_size)
        self.opts_set('cond_ric_alg', opts.cond_ric_alg)
        # general
        self.opts_set('warm_start', opts.warm_start)
        self.opts_set('print_level', opts.print_level)


    def _set_initial_qp_data_in_c(self):
        # void ocp_qp_in_set(ocp_qp_xcond_solver_config *config, ocp_qp_in *in,
                #    int stage, char *field, void *value)
        self.__acados_lib.ocp_qp_in_set.argtypes = [c_void_p, c_void_p, c_int, c_char_p, c_void_p]
        self.__acados_lib.ocp_qp_in_set.restype = None
        N = self.qp.N
        qp = self.qp

        # Define field-value list pairs outside the loop
        fieldname_list_pairs = [('A', qp.A), ('B', qp.B), ('b', qp.b),
                                ('Q', qp.Q), ('S', qp.S), ('R', qp.R), ('q', qp.q), ('r', qp.r),
                                ('idxb', qp.idxb),
                                ('lbx', qp.lbx), ('ubx', qp.ubx),
                                ('lbu', qp.lbu), ('ubu', qp.ubu),
                                ('lls', qp.lls), ('lus', qp.lus),
                                ('lbx_mask', qp.lbx_mask), ('ubx_mask', qp.ubx_mask),
                                ('lbu_mask', qp.lbu_mask), ('ubu_mask', qp.ubu_mask),
                                ('lls_mask', qp.lls_mask), ('lus_mask', qp.lus_mask),
                                ('C', qp.C), ('D', qp.D), ('lg', qp.lg), ('ug', qp.ug),
                                ('Zl', qp.Zl), ('Zu', qp.Zu), ('zl', qp.zl), ('zu', qp.zu),
                                ]
        int_fields = ['idxb']

        for i in range(N + 1):
            for field_name, value_list in fieldname_list_pairs:
                if i == N and field_name in qp.dynamics_fields:
                    continue
                # Get elements in column major order
                if field_name in int_fields:
                    value_ = np.ravel(value_list[i].astype(np.int32), order='F')
                    value_data = cast(value_.ctypes.data, POINTER(c_int))
                else:
                    value_ = np.ravel(value_list[i].astype(float), order='F')
                    value_data = cast(value_.ctypes.data, POINTER(c_double))
                self.__acados_lib.ocp_qp_in_set(self.c_config, self.c_in, c_int(i), field_name.encode('utf-8'), cast(value_data, c_void_p))


    def solve(self) -> int:
        """
        Solve the ocp with current input.

        :return: status of the solver
        """
        self._status = self.__acados_lib.ocp_qp_solve(self.c_solver, self.c_in, self.c_out)

        return self._status


    # def get_dim_flat(self, field: str):
    #     """
    #     Get dimension of flattened iterate.
    #     """
    #     if field not in ['x', 'u', 'z', 'pi', 'lam', 'sl', 'su', 'p']:
    #         raise ValueError(f'AcadosOcpSolver.get_dim_flat(field={field}): \'{field}\' is an invalid argument.')

    #     return self.__acados_lib.ocp_nlp_dims_get_total_from_attr(self.nlp_config, self.nlp_dims, self.nlp_out, field.encode('utf-8'))


    # def reset(self, reset_qp_solver_mem=1):
    #     """
    #     Sets current iterate to all zeros.
    #     """
    #     getattr(self.shared_lib, f"{self.name}_acados_reset")(self.capsule, reset_qp_solver_mem)


    def get(self, stage_: int, field_: str):
        """
        Get the last solution of the solver.

        :param stage: integer - stage of OCP-QP.
        :param field: string in ['x', 'u', 'pi', 'lam', 'sl', 'su']
        """

        out_fields = ['x', 'u', 'pi', 'lam', 'sl', 'su']


        if (field_ not in out_fields):
            raise ValueError(f'AcadosOcpQpSolver.get(stage={stage_}, field={field_}): \'{field_}\' is an invalid argument.'
                             f'\n Possible values are {out_fields}.')

        if not isinstance(stage_, int):
            raise TypeError(f'AcadosOcpQpSolver.get(stage={stage_}, field={field_}): stage index must be an integer, got type {type(stage_)}.')

        if stage_ < 0 or stage_ > self.N:
            raise ValueError(f'AcadosOcpQpSolver.get(stage={stage_}, field={field_}): stage index must be in [0, {self.N}], got: {stage_}.')

        if stage_ == self.N and field_ == 'pi':
            raise KeyError(f'AcadosOcpQpSolver.get(stage={stage_}, field={field_}): field \'{field_}\' does not exist at final stage {stage_}.')

        field = field_.encode('utf-8')

        if field_ == 'x':
            dim = self.qp.dims.nx[stage_]
        elif field_ == 'u':
            dim = self.qp.dims.nu[stage_]
        elif field_ == 'pi':
            dim = self.qp.dims.nx[stage_+1]
        elif field_ == 'lam':
            dim = 2*(self.qp.dims.ng[stage_] + self.qp.dims.nbx[stage_] + self.qp.dims.nbu[stage_] + self.qp.dims.ns[stage_])
        elif field_ in ['sl', 'su']:
            dim = self.qp.dims.ns[stage_]

        out = np.zeros((dim,), dtype=np.float64, order="C")
        out_data = cast(out.ctypes.data, POINTER(c_double))

        self.__acados_lib.ocp_qp_out_get(self.c_out, c_int(stage_), field, cast(out_data, c_void_p))

        return out


    def get_flat(self, field_: str) -> np.ndarray:
        raise NotImplementedError("get_flat() not implemented yet.")


    def set_flat(self, field_: str, value_: np.ndarray) -> None:
        raise NotImplementedError("set_flat() not implemented yet.")

    def print_statistics(self):
        raise NotImplementedError("print_statistics() not implemented yet.")
        return


    def qp_diagnostics(self, hessian_type: str = 'FULL_HESSIAN'):
        raise NotImplementedError("qp_diagnostics() not implemented yet.")
        # TODO: implement something joint for NLP and QP solver

    def get_stats(self, field_: str) -> Union[int, float, np.ndarray]:
        int_fields = ['iter']
        double_fields = ['tau_iter', 'time_qp_solver_call', 'time_qp_xcond', 'time_tot']

        print(f"Getting stats for field '{field_}'...")

        if field_ in int_fields:
            value = c_int()
            value_data = byref(value)

            self.__acados_lib.ocp_qp_xcond_solver_get_scalar(self.c_solver, self.c_out, field_.encode('utf-8'), cast(value_data, c_void_p))
            return value.value
        elif field_ in double_fields:
            value = c_double(0)
            self.__acados_lib.ocp_qp_xcond_solver_get_scalar(self.c_solver, self.c_out, field_.encode('utf-8'), byref(value))
            return value.value
        else:
            raise NotImplementedError(f"get_stats() does not support field '{field_}' yet.")


    def get_cost(self) -> float:
        """
        Returns the cost value of the current solution.
        """
        raise NotImplementedError("get_cost() not implemented yet.")


    def get_residuals(self, recompute=False):
        raise NotImplementedError("get_residuals() not implemented yet.")


    def set(self, stage_: int, field_: str, value_: np.ndarray):
        raise NotImplementedError("set() not implemented yet.")


    def get_iterate(self,) -> AcadosOcpIterate:
        d = {}
        for field in ["x", "u", "sl", "su", "pi", "lam"]:
            traj = []
            for n in range(self.N+1):
                if n < self.N or not (field in ["pi"]):
                    traj.append(self.get(n, field))
            d[field] = traj
        d["z"] = self.N * [np.array([])]  # z not supported
        return AcadosOcpIterate(**d)
