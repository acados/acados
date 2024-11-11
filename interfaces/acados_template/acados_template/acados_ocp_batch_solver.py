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

from .acados_ocp_solver import AcadosOcpSolver
from .acados_ocp import AcadosOcp
from typing import List
from ctypes import (POINTER, c_int, c_void_p, cast, c_double)
import numpy as np

class AcadosOcpBatchSolver():
    """
    Batch Integrator for parallel integration.

        :param sim: type :py:class:`~acados_template.acados_sim.AcadosOcp`
        :param N_batch: batch size, positive integer
        :param json_file: Default: 'acados_sim.json'
        :verbose: bool, default: True
    """

    __ocp_solvers : List[AcadosOcpSolver]

    def __init__(self, ocp: AcadosOcp, N_batch: int, json_file: str = 'acados_ocp.json', verbose: bool=True):

        if not isinstance(N_batch, int) or N_batch <= 0:
            raise Exception("AcadosOcpBatchSolver: argument N_batch should be a positive integer.")

        self.__N_batch = N_batch
        self.__ocp_solvers = [AcadosOcpSolver(ocp, json_file=json_file, build=n==0, generate=n==0, verbose=verbose) for n in range(self.N_batch)]

        self.__shared_lib = self.ocp_solvers[0].shared_lib
        self.__name = self.ocp_solvers[0].name
        self.__ocp_solvers_pointer = (c_void_p * self.N_batch)()

        for i in range(self.N_batch):
            self.__ocp_solvers_pointer[i] = self.ocp_solvers[i].capsule

        getattr(self.__shared_lib, f"{self.__name}_acados_batch_solve").argtypes = [POINTER(c_void_p), c_int]
        getattr(self.__shared_lib, f"{self.__name}_acados_batch_solve").restype = c_void_p

        if self.ocp_solvers[0].acados_lib_uses_omp:
            msg = "Note: Please make sure that the acados shared library is compiled with the number of threads set to 1,\n"
        else:
            msg = "Warning: Please compile the acados shared library with openmp and the number of threads set to 1,\n"

        msg += "i.e. with the flags -DACADOS_WITH_OPENMP=ON -DACADOS_NUM_THREADS=1.\n" + \
                   "See https://github.com/acados/acados/pull/1089 for more details."
        print(msg)


    @property
    def ocp_solvers(self):
        """List of AcadosOcpSolvers."""
        return self.__ocp_solvers


    @property
    def N_batch(self):
        """Batch size."""
        return self.__N_batch


    def solve(self):
        """
        Call solve for all `N_batch` solvers.
        """
        getattr(self.__shared_lib, f"{self.__name}_acados_batch_solve")(self.__ocp_solvers_pointer, self.__N_batch)


    def set_flat(self, field_: str, value_: np.ndarray) -> None:
        """
        Set concatenation solver initialization for all `N_batch` solvers.

            :param field_: string in ['x', 'u', 'z', 'pi', 'lam', 'sl', 'su', 'p']
            :param value_: np.array of shape (N_batch, n_field)
        """

        field = field_.encode('utf-8')
        if field_ not in ['x', 'u', 'z', 'pi', 'lam', 'sl', 'su', 'p']:
            raise Exception(f'AcadosOcpSolver.get_flat(field={field_}): \'{field_}\' is an invalid argument.')

        dim = self.ocp_solvers[0].get_dim_flat(field_)

        if value_.shape != (self.N_batch, dim):
            raise Exception(f'AcadosOcpBatchSolver.set_flat(field={field_}, value): value has wrong shape, expected ({self.N_batch}, {dim}), got {value_.shape}.')

        value_ = value_.reshape((-1,), order='C')
        N_data = value_.shape[0]

        value_ = value_.astype(float)
        value_data = cast(value_.ctypes.data, POINTER(c_double))
        value_data_p = cast((value_data), c_void_p)

        getattr(self.__shared_lib, f"{self.__name}_acados_batch_set_flat")(self.__ocp_solvers_pointer, field, value_data_p, N_data, self.__N_batch)

