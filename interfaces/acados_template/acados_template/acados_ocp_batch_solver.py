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
from typing import Optional, List, Tuple
from collections.abc import Sequence
from ctypes import (POINTER, c_int, c_void_p, cast, c_double, c_char_p)
import numpy as np
import time

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
        self.__acados_lib = self.ocp_solvers[0].acados_lib
        self.__name = self.ocp_solvers[0].name
        self.__ocp_solvers_pointer = (c_void_p * self.N_batch)()

        for i in range(self.N_batch):
            self.__ocp_solvers_pointer[i] = self.ocp_solvers[i].capsule

        getattr(self.__shared_lib, f"{self.__name}_acados_batch_solve").argtypes = [POINTER(c_void_p), c_int]
        getattr(self.__shared_lib, f"{self.__name}_acados_batch_solve").restype = c_void_p

        getattr(self.__shared_lib, f"{self.__name}_acados_batch_eval_params_jac").argtypes = [POINTER(c_void_p), c_int]
        getattr(self.__shared_lib, f"{self.__name}_acados_batch_eval_params_jac").restype = c_void_p

        getattr(self.__shared_lib, f"{self.__name}_acados_batch_eval_solution_sens_adj_p").argtypes = [POINTER(c_void_p), c_char_p, c_int, POINTER(c_double), c_int, c_int]
        getattr(self.__shared_lib, f"{self.__name}_acados_batch_eval_solution_sens_adj_p").restype = c_void_p

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



    def eval_adjoint_solution_sensitivity(self,
                                          seed_x: Optional[Sequence[Tuple[int, np.ndarray]]],
                                          seed_u: Optional[Sequence[Tuple[int, np.ndarray]]],
                                          with_respect_to: str = "p_global",
                                          sanity_checks: bool = True,
                                          ) -> np.ndarray:
        """
        Evaluate the adjoint sensitivity of the solution with respect to the parameters.
            :param seed_x : Sequence of tuples of the form (stage: int, seed_vec: np.ndarray).
                    The stage is the stage at which the seed_vec is applied, and seed_vec is the seed for the states at that stage with shape (N_batch, nx, n_seeds)
            :param seed_u : Sequence of tuples of the form (stage: int, seed_vec: np.ndarray).
                    The stage is the stage at which the seed_vec is applied, and seed_vec is the seed for the controls at that stage with shape (N_batch, nu, n_seeds)
            :param with_respect_to : string in ["p_global"]
            :param sanity_checks : bool - whether to perform sanity checks, turn off for minimal overhead, default: True
        """

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
            n_seeds = seed_x[0][1].shape[2]

        if len(seed_u) > 0:
            if not isinstance(seed_u[0], tuple) or len(seed_u[0]) != 2:
                raise Exception(f"seed_u[0] should be tuple of length 2, got seed_u[0] {seed_u[0]}")
            s = seed_u[0][1]
            if not isinstance(s, np.ndarray):
                raise Exception(f"seed_u[0][1] should be np.ndarray, got {type(s)}")
            n_seeds = seed_u[0][1].shape[2]

        if sanity_checks:
            N_horizon = self.__ocp_solvers[0].acados_ocp.solver_options.N_horizon

            for n in range(self.N_batch):
                self.__ocp_solvers[n].sanity_check_parametric_sensitivities()

            nx = self.__acados_lib.ocp_nlp_dims_get_from_attr(
                self.__ocp_solvers[0].nlp_config, self.__ocp_solvers[0].nlp_dims, self.__ocp_solvers[0].nlp_out, 0, "x".encode('utf-8'))
            nu = self.__acados_lib.ocp_nlp_dims_get_from_attr(
                self.__ocp_solvers[0].nlp_config, self.__ocp_solvers[0].nlp_dims, self.__ocp_solvers[0].nlp_out, 0, "u".encode('utf-8'))

            # check seeds
            for seed, name, dim in [(seed_x, "seed_x", nx,), (seed_u, "seed_u", nu)]:
                for stage, seed_stage in seed:
                    if not isinstance(stage, int) or stage < 0 or stage > N_horizon:
                        raise Exception(f"AcadosOcpSolver.eval_solution_sensitivity(): stage {stage} for {name} is not valid.")
                    if not isinstance(seed_stage, np.ndarray):
                        raise Exception(f"{name} for stage {stage} should be np.ndarray, got {type(seed_stage)}")
                    if seed_stage.shape != (self.N_batch, dim, n_seeds):
                        raise Exception(f"{name} for stage {stage} should have shape (N_batch, dim, n_seeds) = ({self.N_batch}, {dim}, {n_seeds}), got {seed_stage.shape}.")

        if with_respect_to == "p_global":
            field = "p_global".encode('utf-8')

            np_global = self.__acados_lib.ocp_nlp_dims_get_from_attr(
                self.__ocp_solvers[0].nlp_config, self.__ocp_solvers[0].nlp_dims, self.__ocp_solvers[0].nlp_out, 0, field)

            # compute jacobian wrt params
            t0 = time.time()
            getattr(self.__shared_lib, f"{self.__name}_acados_batch_eval_params_jac")(self.__ocp_solvers_pointer, self.__N_batch)
            self.time_solution_sens_lin = time.time() - t0

            t1 = time.time()

            out = np.zeros((self.N_batch, n_seeds, np_global))
            for i_seed in range(n_seeds):
                # set seed:
                self.reset_sens_out()

                # TODO this could be a batch_set
                for m in range(self.N_batch):
                    for (stage, sx) in seed_x:
                        self.ocp_solvers[m].set(stage, 'sens_x', sx[m, :, i_seed])
                    for (stage, su) in seed_u:
                        self.ocp_solvers[m].set(stage, 'sens_u', su[m, :, i_seed])


                grad = np.zeros((self.N_batch, np_global))
                grad_p = np.ascontiguousarray(grad, dtype=np.float64)

                c_grad_p = cast(grad_p.ctypes.data, POINTER(c_double))

                # solve adjoint sensitivities
                getattr(self.__shared_lib, f"{self.__name}_acados_batch_eval_solution_sens_adj_p")(
                    self.__ocp_solvers_pointer, field, 0, c_grad_p, np_global, self.__N_batch)

                out[:, i_seed, :] = grad_p

            self.time_solution_sens_solve = time.time() - t1

            return out
        else:
            raise NotImplementedError(f"with_respect_to {with_respect_to} not implemented.")


    def reset_sens_out(self, ):
        # TODO batch this
        for n in range(self.N_batch):
            self.__acados_lib.ocp_nlp_out_set_values_to_zero(self.__ocp_solvers[n].nlp_config, self.__ocp_solvers[n].nlp_dims, self.__ocp_solvers[n].sens_out)



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

