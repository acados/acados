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

from .acados_sim_solver import AcadosSimSolver
from .acados_sim import AcadosSim
from typing import List, Optional
from ctypes import (POINTER, c_int, c_void_p)
import warnings


class AcadosSimBatchSolver():
    """
    Batch Integrator for parallel integration.

        :param sim: type :py:class:`~acados_template.acados_sim.AcadosSim`
        :param N_batch_init: initial batch size, positive integer
        :param num_threads_in_batch_solve: number of threads used for parallelizing the batch methods. Default: 1
        :param build: Flag indicating whether solver should be (re)compiled. If False an attempt is made to load an already compiled shared library for the solver. Default: True
        :param generate: Flag indicating whether problem functions should be code generated. Default: True
        :verbose: bool, default: True
    """

    __sim_solvers : List[AcadosSimSolver]

    def __init__(self, sim: AcadosSim,
                 N_batch_init: int, num_threads_in_batch_solve: int = 1,
                 json_file: Optional[str] = None,
                 build: bool = True,
                 generate: bool = True,
                 verbose: bool=True):

        if not isinstance(N_batch_init, int) or N_batch_init <= 0:
            raise ValueError("AcadosSimBatchSolver: argument N_batch_init should be a positive integer.")
        if not isinstance(num_threads_in_batch_solve, int) or num_threads_in_batch_solve <= 0:
            raise ValueError("AcadosSimBatchSolver: argument num_threads_in_batch_solve should be a positive integer.")
        if not sim.solver_options.with_batch_functionality:
            warnings.warn("Using AcadosSimBatchSolver, but sim.solver_options.with_batch_functionality is False. Attempting to compile with openmp nonetheless.")
            sim.solver_options.with_batch_functionality = True

        if json_file is not None:
            warnings.warn(" The `json_file` argument is deprecated in v0.5.6 and will be removed in a future release. Set AcadosSim.code_gen_options.json_file instead.", DeprecationWarning, stacklevel=2)
            sim.code_gen_options.json_file = json_file

        self.__num_threads_in_batch_solve = num_threads_in_batch_solve
        self.__n_batch_current = N_batch_init
        self.__sim_solvers = [AcadosSimSolver(sim,
                                              build=n==0 if build else False,
                                              generate=n==0 if generate else False,
                                              verbose=verbose if n==0 else False,
                                              )
                              for n in range(self.n_batch_current)]

        self.__shared_lib = self.sim_solvers[0].shared_lib
        self.__name = self.sim_solvers[0].sim.name
        self.__sim_solvers_pointer = (c_void_p * self.n_batch_current)()

        for i in range(self.n_batch_current):
            self.__sim_solvers_pointer[i] = self.sim_solvers[i].capsule

        getattr(self.__shared_lib, f"{self.__name}_acados_sim_batch_solve").argtypes = [POINTER(c_void_p), c_int, c_int]
        getattr(self.__shared_lib, f"{self.__name}_acados_sim_batch_solve").restype = c_void_p

        if not self.sim_solvers[0].acados_lib_uses_omp:
            warnings.warn("Please compile the acados shared library with openmp and the number of threads set to 1, i.e. with the flags -DACADOS_WITH_OPENMP=ON -DACADOS_NUM_THREADS=1.")


    def solve(self, n_batch: int = None) -> None:
        """
        Solve the simulation problem with current input for the first `n_batch` integrators.
        Or `n_batch_current` if `n_batch` is None.
        """
        if n_batch is None:
            n_batch = self.n_batch_current
        elif n_batch <= self.n_batch_current:
            self.__n_batch_current = n_batch
        else:
            raise ValueError("You are attempting to solve more problem instances than what have been initialized so far.")

        getattr(self.__shared_lib, f"{self.__name}_acados_sim_batch_solve")(self.__sim_solvers_pointer, n_batch, self.__num_threads_in_batch_solve)


    @property
    def sim_solvers(self):
        """List of AcadosSimSolvers."""
        return self.__sim_solvers

    @property
    def n_batch_current(self):
        """The current default batch size which depends on the previous batch method call."""
        return self.__n_batch_current

    @property
    def num_threads_in_batch_solve(self):
        """Number of threads used for parallelizing the batch methods."""
        return self.__num_threads_in_batch_solve
    
    @num_threads_in_batch_solve.setter
    def num_threads_in_batch_solve(self, num_threads_in_batch_solve):
        self.__num_threads_in_batch_solve = num_threads_in_batch_solve
