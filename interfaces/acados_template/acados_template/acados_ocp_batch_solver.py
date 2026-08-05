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
from .acados_ocp_iterate import AcadosOcpFlattenedBatchIterate
from typing import Optional, List, Tuple, Sequence, Union, Dict
from ctypes import (POINTER, c_int, c_void_p, cast, c_double, c_char_p)
import numpy as np
import time
import warnings
from deprecated.sphinx import deprecated


class AcadosOcpBatchSolver():
    """
    Batch OCP solver for parallel solves.

    :param ocp: type :py:class:`~acados_template.acados_ocp.AcadosOcp`
    :param N_batch_init: initial batch size, batch size can change dynamically, positive integer
    :param num_threads_in_batch_solve: number of threads used for parallelizing the batch methods. Default: 1
    :param json_file: Default: 'acados_ocp.json'
    :param build: Flag indicating whether solver should be (re)compiled. If False, an attempt is made to load an already compiled shared library for the solver. Default: True
    :param generate: Flag indicating whether problem functions should be code generated. Default: True
    :param verbose: bool, default: True
    :param save_p_global: bool, default: False
    :param check_code_reuse_possible: If generate or build is false, compares the data in the json_file to the ocp object and sets generate or build to True if necessary, Default: True

    """

    __ocp_solvers : List[AcadosOcpSolver]

    def __init__(self, ocp: AcadosOcp, N_batch_init: int,
                 num_threads_in_batch_solve: int = 1,
                 json_file: Optional[str] = None,
                 build: bool = True, generate: bool = True,
                 verbose: bool = True, check_code_reuse_possible: bool = True,
                 save_p_global: bool = False):

        if not isinstance(N_batch_init, int) or N_batch_init <= 0:
            raise ValueError("AcadosOcpBatchSolver: argument N_batch_init should be a positive integer.")
        if not isinstance(num_threads_in_batch_solve, int) or num_threads_in_batch_solve <= 0:
            raise ValueError("AcadosOcpBatchSolver: argument num_threads_in_batch_solve should be a positive integer.")
        if not ocp.solver_options.with_batch_functionality:
            warnings.warn("Using AcadosOcpBatchSolver, but ocp.solver_options.with_batch_functionality is False. Setting to True and attempting to compile with openmp nonetheless.")
            ocp.solver_options.with_batch_functionality = True

        if json_file is not None:
            warnings.warn(" The `json_file` argument is deprecated in v0.5.6 and will be removed in a future release. Set AcadosOcp.code_gen_options.json_file instead.", DeprecationWarning, stacklevel=2)
            ocp.code_gen_options.json_file = json_file

        self.__num_threads_in_batch_solve = num_threads_in_batch_solve

        self.__n_batch_current = N_batch_init
        self.__ocp_solvers = [AcadosOcpSolver(ocp,
                                              build=n==0 if build else False,
                                              generate=n==0 if generate else False,
                                              verbose=n==0 if verbose else False,
                                              save_p_global=save_p_global,
                                              check_reuse_possible=n==0 if check_code_reuse_possible else False,
                                              )
                               for n in range(self.n_batch_current)]

        self.__shared_lib = self.ocp_solvers[0].shared_lib
        self.__acados_lib = self.ocp_solvers[0].acados_lib
        self.__name = self.ocp_solvers[0].ocp.name
        self.__ocp_solvers_pointer = (c_void_p * self.n_batch_current)()

        for i in range(self.n_batch_current):
            self.__ocp_solvers_pointer[i] = self.ocp_solvers[i].capsule

        # out data for solve
        self.__status = np.zeros((self.n_batch_current,), dtype=np.intc, order="C")
        self.__status_p = cast(self.__status.ctypes.data, POINTER(c_int))

        getattr(self.__shared_lib, f"{self.__name}_acados_batch_solve").argtypes = [POINTER(c_void_p), POINTER(c_int), c_int, c_int]
        getattr(self.__shared_lib, f"{self.__name}_acados_batch_solve").restype = c_void_p

        getattr(self.__shared_lib, f"{self.__name}_acados_batch_eval_params_jac").argtypes = [POINTER(c_void_p), c_int, c_int]
        getattr(self.__shared_lib, f"{self.__name}_acados_batch_eval_params_jac").restype = c_void_p

        getattr(self.__shared_lib, f"{self.__name}_acados_batch_eval_param_sens").argtypes = [POINTER(c_void_p), c_char_p, c_int, c_int, c_int, c_int]
        getattr(self.__shared_lib, f"{self.__name}_acados_batch_eval_param_sens").restype = c_void_p

        getattr(self.__shared_lib, f"{self.__name}_acados_batch_eval_solution_sens_adj_p").argtypes = [POINTER(c_void_p), c_char_p, c_int, POINTER(c_double), c_int, c_int, c_int]
        getattr(self.__shared_lib, f"{self.__name}_acados_batch_eval_solution_sens_adj_p").restype = c_void_p

        getattr(self.__shared_lib, f"{self.__name}_acados_batch_set_flat").argtypes = [POINTER(c_void_p), c_char_p, POINTER(c_double), c_int, c_int, c_int]
        getattr(self.__shared_lib, f"{self.__name}_acados_batch_set_flat").restype = c_void_p

        getattr(self.__shared_lib, f"{self.__name}_acados_batch_get_flat").argtypes = [POINTER(c_void_p), c_char_p, POINTER(c_double), c_int, c_int, c_int]
        getattr(self.__shared_lib, f"{self.__name}_acados_batch_get_flat").restype = c_void_p

        getattr(self.__shared_lib, f"{self.__name}_acados_batch_set").argtypes = [POINTER(c_void_p), c_char_p, c_int, POINTER(c_double), c_int, c_int, c_int]
        getattr(self.__shared_lib, f"{self.__name}_acados_batch_set").restype = c_void_p

        getattr(self.__shared_lib, f"{self.__name}_acados_batch_get").argtypes = [POINTER(c_void_p), c_char_p, c_int, POINTER(c_double), c_int, c_int, c_int]
        getattr(self.__shared_lib, f"{self.__name}_acados_batch_get").restype = c_void_p

        getattr(self.__shared_lib, f"{self.__name}_acados_batch_constraints_set").argtypes = [POINTER(c_void_p), c_char_p, c_int, POINTER(c_double), c_int, c_int, c_int]
        getattr(self.__shared_lib, f"{self.__name}_acados_batch_constraints_set").restype = c_void_p

        getattr(self.__shared_lib, f"{self.__name}_acados_batch_reset_sens_out").argtypes = [POINTER(c_void_p), c_int, c_int]
        getattr(self.__shared_lib, f"{self.__name}_acados_batch_reset_sens_out").restype = c_void_p

        if self.ocp_solvers[0].acados_lib_uses_omp:
            msg = "Note: Please make sure that the acados shared library is compiled with the number of threads set to 1,\n"
        else:
            msg = "Warning: Please compile the acados shared library with openmp and the number of threads set to 1,\n"

        msg += "i.e. with the flags -DACADOS_WITH_OPENMP=ON -DACADOS_NUM_THREADS=1.\n" + \
                   "See https://github.com/acados/acados/pull/1089 for more details."

        if verbose:
            print(msg)
        self.verbose = verbose
        self.time_solution_sens_lin = 0.0
        self.time_solution_sens_solve = 0.0


    @property
    def ocp_solvers(self):
        """List of AcadosOcpSolvers."""
        return self.__ocp_solvers

    @property
    def num_threads_in_batch_solve(self):
        """Number of threads used for parallelizing the batch methods."""
        return self.__num_threads_in_batch_solve

    @num_threads_in_batch_solve.setter
    def num_threads_in_batch_solve(self, num_threads_in_batch_solve):
        self.__num_threads_in_batch_solve = num_threads_in_batch_solve

    @property
    def status(self):
        """
        Returns the status of the last batch_solver call. Has shape (`n_batch_current`, ).
        NOTE: If you make single solver calls, e.g., via `batch_solver.ocp_solvers[i].solve()`, this status will not reflect them.

        Status codes:
            - 0: Success (ACADOS_SUCCESS)
            - 1: NaN detected (ACADOS_NAN_DETECTED)
            - 2: Maximum number of iterations reached (ACADOS_MAXITER)
            - 3: Minimum step size reached (ACADOS_MINSTEP)
            - 4: QP solver failed (ACADOS_QP_FAILURE)
            - 5: Solver created (ACADOS_READY)
            - 6: Problem unbounded (ACADOS_UNBOUNDED)
            - 7: Solver timeout (ACADOS_TIMEOUT)
            - 8: QP scaling could not satisfy bounds (ACADOS_QPSCALING_BOUNDS_NOT_SATISFIED); NOTE: this status is typically not returned by the solver, but can be checked via `get_stats('qpscaling_status')`

        See `return_values` in https://github.com/acados/acados/blob/main/acados/utils/types.h
        """
        return self.__status[:self.n_batch_current].copy()

    @property
    def n_batch_current(self):
        """The current default batch size which depends on the previous batch method call."""
        return self.__n_batch_current


    def solve(self, n_batch: Optional[int] = None) -> None:
        """
        Call solve for the first `n_batch` solvers. Or `n_batch_current` if `n_batch` is None.
        """
        if n_batch is None:
            n_batch = self.n_batch_current
        elif n_batch <= self.n_batch_current:
            self.__n_batch_current = n_batch
        else:
            raise ValueError("You are attempting to solve more problem instances than what have been initialized so far. "
                             "First initialize enough problem instances by using the setter methods.")

        getattr(self.__shared_lib, f"{self.__name}_acados_batch_solve")(self.__ocp_solvers_pointer, self.__status_p, n_batch, self.__num_threads_in_batch_solve)

        # to be consistent with non-batched solve
        for s, solver in zip(self.__status, self.ocp_solvers):
            solver._status = s


    def setup_qp_matrices_and_factorize(self, n_batch: Optional[int] = None) -> None:
        """
        Call setup_qp_matrices_and_factorize for the first `n_batch` solvers. Or `n_batch_current` if `n_batch` is None.
        """
        if n_batch is None:
            n_batch = self.n_batch_current
        elif n_batch <= self.n_batch_current:
            self.__n_batch_current = n_batch
        else:
            raise ValueError("You are attempting to factorize for more problem instances than what have been initialized so far. "
                             "First initialize enough problem instances by using the setter methods.")

        getattr(self.__shared_lib, f"{self.__name}_acados_batch_setup_qp_matrices_and_factorize")(self.__ocp_solvers_pointer, self.__status_p, n_batch, self.__num_threads_in_batch_solve)

        # to be consistent with non-batched solve
        for s, solver in zip(self.__status, self.ocp_solvers):
            solver._status = s


    def eval_adjoint_solution_sensitivity(self,
                                          seed_x: Optional[Sequence[Tuple[int, np.ndarray]]],
                                          seed_u: Optional[Sequence[Tuple[int, np.ndarray]]],
                                          with_respect_to: str = "p_global",
                                          sanity_checks: bool = True,
                                          ) -> np.ndarray:
        """
        Evaluate the adjoint sensitivity of the solution with respect to the parameters.
            :param seed_x : Sequence of tuples of the form (stage: int, seed_vec: np.ndarray).
                    The stage is the stage at which the seed_vec is applied, and seed_vec is the seed for the states at that stage with shape (n_batch, nx, n_seeds)
            :param seed_u : Sequence of tuples of the form (stage: int, seed_vec: np.ndarray).
                    The stage is the stage at which the seed_vec is applied, and seed_vec is the seed for the controls at that stage with shape (n_batch, nu, n_seeds)
            :param with_respect_to : string in ["p_global"]
            :param sanity_checks : bool - whether to perform sanity checks on solver settings, turn off for minimal overhead, default: True
            :returns : np.ndarray of shape (n_batch, n_seeds, np_global)
        """

        if seed_x is None:
            seed_x = []
        elif not isinstance(seed_x, Sequence):
            raise TypeError(f"seed_x should be a Sequence, got {type(seed_x)}")

        if seed_u is None:
            seed_u = []
        elif not isinstance(seed_u, Sequence):
            raise TypeError(f"seed_u should be a Sequence, got {type(seed_u)}")

        if len(seed_x) == 0 and len(seed_u) == 0:
            raise ValueError("seed_x and seed_u cannot both be empty.")

        if len(seed_x) > 0:
            if not isinstance(seed_x[0], tuple) or len(seed_x[0]) != 2:
                raise TypeError(f"seed_x[0] should be tuple of length 2, got seed_x[0] {seed_x[0]}")
            s = seed_x[0][1]
            if not isinstance(s, np.ndarray):
                raise TypeError(f"seed_x[0][1] should be np.ndarray, got {type(s)}")
            n_seeds = seed_x[0][1].shape[2]
            n_batch = seed_x[0][1].shape[0]

        if len(seed_u) > 0:
            if not isinstance(seed_u[0], tuple) or len(seed_u[0]) != 2:
                raise ValueError(f"seed_u[0] should be tuple of length 2, got seed_u[0] {seed_u[0]}")
            s = seed_u[0][1]
            if not isinstance(s, np.ndarray):
                raise TypeError(f"seed_u[0][1] should be np.ndarray, got {type(s)}")
            n_seeds = seed_u[0][1].shape[2]
            n_batch = seed_u[0][1].shape[0]

        if n_batch <= self.n_batch_current:
            self.__n_batch_current = n_batch
        else:
            raise ValueError("You are attempting to obtain more sensitivities than problem instances initialized so far. "
                             "First initialize enough problem instances by using the setter methods.")

        if sanity_checks:
            N_horizon = self.__ocp_solvers[0].acados_ocp.solver_options.N_horizon

            for n in range(n_batch):
                self.__ocp_solvers[n]._ensure_solution_sensitivities_available()

            nx = self.__acados_lib.ocp_nlp_dims_get_from_attr(
                self.__ocp_solvers[0].nlp_config, self.__ocp_solvers[0].nlp_dims, self.__ocp_solvers[0].nlp_out, 0, "x".encode('utf-8'))
            nu = self.__acados_lib.ocp_nlp_dims_get_from_attr(
                self.__ocp_solvers[0].nlp_config, self.__ocp_solvers[0].nlp_dims, self.__ocp_solvers[0].nlp_out, 0, "u".encode('utf-8'))

            # check seeds
            for seed, name, dim in [(seed_x, "seed_x", nx,), (seed_u, "seed_u", nu)]:
                for stage, seed_stage in seed:
                    if not isinstance(stage, int) or stage < 0 or stage > N_horizon:
                        raise ValueError(f"AcadosOcpBatchSolver.eval_adjoint_solution_sensitivity(): stage {stage} for {name} is not valid.")
                    if not isinstance(seed_stage, np.ndarray):
                        raise TypeError(f"{name} for stage {stage} should be np.ndarray, got {type(seed_stage)}")
                    if seed_stage.shape != (n_batch, dim, n_seeds):
                        raise ValueError(f"{name} for stage {stage} should have shape (n_batch, dim, n_seeds) = ({n_batch}, {dim}, {n_seeds}), got {seed_stage.shape}.")

        if with_respect_to == "p_global":
            field = "p_global".encode('utf-8')

            np_global = self.__acados_lib.ocp_nlp_dims_get_from_attr(
                self.__ocp_solvers[0].nlp_config, self.__ocp_solvers[0].nlp_dims, self.__ocp_solvers[0].nlp_out, 0, field)

            # compute jacobian wrt params
            t1 = time.time()

            grad_p = np.zeros((n_batch, n_seeds, np_global), order="C", dtype=np.float64)
            offset = n_seeds*np_global

            for i_seed in range(n_seeds):
                # set seed:
                self._reset_sens_out(n_batch)

                for (stage, sx) in seed_x:
                    self.set(stage, 'sens_x', sx[:n_batch, :, i_seed])
                for (stage, su) in seed_u:
                    self.set(stage, 'sens_u', su[:n_batch, :, i_seed])

                c_grad_p = cast(grad_p[0, i_seed].ctypes.data, POINTER(c_double))

                # solve adjoint sensitivities
                getattr(self.__shared_lib, f"{self.__name}_acados_batch_eval_solution_sens_adj_p")(
                    self.__ocp_solvers_pointer, field, 0, c_grad_p, offset, n_batch, self.__num_threads_in_batch_solve)

            self.time_solution_sens_solve = time.time() - t1

            return grad_p
        else:
            raise NotImplementedError(f"with_respect_to {with_respect_to} not implemented.")


    def _reset_sens_out(self, n_batch: int) -> None:
        getattr(self.__shared_lib, f"{self.__name}_acados_batch_reset_sens_out")(
            self.__ocp_solvers_pointer, c_int(n_batch), c_int(self.__num_threads_in_batch_solve)
        )


    def eval_solution_sensitivity(self,
                                  stages: Union[int, List[int]],
                                  with_respect_to: str,
                                  return_sens_x: bool = False,
                                  return_sens_u: bool = True,
                                  return_sens_pi: bool = False,
                                  return_sens_lam: bool = False,
                                  return_sens_su: bool = False,
                                  return_sens_sl: bool = False,
                                  sanity_checks: bool = True,
                                  ) -> Dict:
        """
        Evaluate forward sensitivities of the current batched solutions with respect to the initial state or global parameters.

        :param stages: stages for which sensitivities are returned, int or list of int
        :param with_respect_to: string in ["initial_state", "p_global"]
        :param return_sens_x: whether to return sensitivities of x. Default: True.
        :param return_sens_u: whether to return sensitivities of u. Default: True.
        :param return_sens_pi: whether to return sensitivities of pi. Default: False.
        :param return_sens_lam: whether to return sensitivities of lam. Default: False.
        :param return_sens_su: whether to return sensitivities of su. Default: False.
        :param return_sens_sl: whether to return sensitivities of sl. Default: False.
        :param sanity_checks: bool - whether to perform sanity checks, turn off for minimal overhead, default: True

        :returns: dictionary with requested sensitivity fields.
                  For each requested field and stage, entries have shape (n_batch, field_dim, ngrad).
                  If stages is a list, each dictionary value is a list over stages.
                  If stages is a scalar, each dictionary value is a single np.ndarray.
        """

        stages_is_list = isinstance(stages, list)
        stages_ = stages if stages_is_list else [stages]

        n_batch = self.n_batch_current
        solver0 = self.__ocp_solvers[0]
        N_horizon = solver0.N

        sens_x = []
        sens_u = []
        sens_pi = []
        sens_lam = []
        sens_sl = []
        sens_su = []

        if sanity_checks:
            for s in stages_:
                if not isinstance(s, int) or s < 0 or s > N_horizon:
                    raise TypeError(
                        "AcadosOcpBatchSolver.eval_solution_sensitivity(): stages need to be int or list[int] "
                        f"and in [0, N], got stages = {stages_}."
                    )

        if with_respect_to == "initial_state":
            nx = self.__acados_lib.ocp_nlp_dims_get_from_attr(
                solver0.nlp_config, solver0.nlp_dims, solver0.nlp_out, 0, "x".encode('utf-8'))
            ngrad = nx
            field = "ex".encode('utf-8')
            if sanity_checks:
                self.__ocp_solvers[0]._ensure_solution_sensitivities_available(parametric=False)
            self.time_solution_sens_lin = 0.0

        elif with_respect_to == "p_global":
            np_global = self.__acados_lib.ocp_nlp_dims_get_from_attr(
                solver0.nlp_config, solver0.nlp_dims, solver0.nlp_out, 0, "p_global".encode('utf-8'))
            ngrad = np_global
            field = "p_global".encode('utf-8')
            if sanity_checks:
                self.__ocp_solvers[0]._ensure_solution_sensitivities_available(forward=True, parametric=True)

            # Compute jacobians wrt params for all batch instances.
            t0 = time.time()
            getattr(self.__shared_lib, f"{self.__name}_acados_batch_eval_params_jac")(
                self.__ocp_solvers_pointer, c_int(n_batch), c_int(self.__num_threads_in_batch_solve)
            )
            self.time_solution_sens_lin = time.time() - t0

        else:
            raise ValueError(
                f"AcadosOcpBatchSolver.eval_solution_sensitivity(): Unknown field: with_respect_to = {with_respect_to}"
            )

        # Initialize output arrays.
        for s in stages_:
            if return_sens_x:
                nx = self.__acados_lib.ocp_nlp_dims_get_from_attr(
                    solver0.nlp_config, solver0.nlp_dims, solver0.nlp_out, s, "x".encode('utf-8'))
                sens_x.append(np.zeros((n_batch, nx, ngrad), dtype=np.float64, order="C"))

            if return_sens_lam:
                nlam = self.__acados_lib.ocp_nlp_dims_get_from_attr(
                    solver0.nlp_config, solver0.nlp_dims, solver0.nlp_out, s, "lam".encode('utf-8'))
                sens_lam.append(np.zeros((n_batch, nlam, ngrad), dtype=np.float64, order="C"))

            if return_sens_sl:
                ns = self.__acados_lib.ocp_nlp_dims_get_from_attr(
                    solver0.nlp_config, solver0.nlp_dims, solver0.nlp_out, s, "s".encode('utf-8'))
                sens_sl.append(np.zeros((n_batch, ns, ngrad), dtype=np.float64, order="C"))

            if return_sens_su:
                ns = self.__acados_lib.ocp_nlp_dims_get_from_attr(
                    solver0.nlp_config, solver0.nlp_dims, solver0.nlp_out, s, "s".encode('utf-8'))
                sens_su.append(np.zeros((n_batch, ns, ngrad), dtype=np.float64, order="C"))

            if s < N_horizon:
                if return_sens_u:
                    nu = self.__acados_lib.ocp_nlp_dims_get_from_attr(
                        solver0.nlp_config, solver0.nlp_dims, solver0.nlp_out, s, "u".encode('utf-8'))
                    sens_u.append(np.zeros((n_batch, nu, ngrad), dtype=np.float64, order="C"))

                if return_sens_pi:
                    npi = self.__acados_lib.ocp_nlp_dims_get_from_attr(
                        solver0.nlp_config, solver0.nlp_dims, solver0.nlp_out, s, "pi".encode('utf-8'))
                    sens_pi.append(np.zeros((n_batch, npi, ngrad), dtype=np.float64, order="C"))

        t0 = time.time()
        for k in range(ngrad):
            # Evaluate sensitivity for all batch instances.
            getattr(self.__shared_lib, f"{self.__name}_acados_batch_eval_param_sens")(
                self.__ocp_solvers_pointer, field, c_int(0), c_int(k), c_int(n_batch), c_int(self.__num_threads_in_batch_solve)
            )

            # Extract sensitivities for requested stages.
            for n, s in enumerate(stages_):
                if return_sens_x:
                    sens_x[n][:, :, k] = self.get(s, "sens_x", n_batch=n_batch)
                if return_sens_lam:
                    sens_lam[n][:, :, k] = self.get(s, "sens_lam", n_batch=n_batch)
                if return_sens_sl:
                    sens_sl[n][:, :, k] = self.get(s, "sens_sl", n_batch=n_batch)
                if return_sens_su:
                    sens_su[n][:, :, k] = self.get(s, "sens_su", n_batch=n_batch)

                if s < N_horizon:
                    if return_sens_u:
                        sens_u[n][:, :, k] = self.get(s, "sens_u", n_batch=n_batch)
                    if return_sens_pi:
                        sens_pi[n][:, :, k] = self.get(s, "sens_pi", n_batch=n_batch)

        self.time_solution_sens_solve = time.time() - t0

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


    def get(self, stage_: int, field_: str, n_batch: Optional[int] = None) -> np.ndarray:
        """
        Get field values at one stage for a batch of solvers.

        :param stage_: integer corresponding to shooting node
        :param field_: string in ['x', 'u', 'z', 'pi', 'lam', 'sl', 'su', 'p',
                       'sens_u', 'sens_pi', 'sens_x', 'sens_lam', 'sens_sl', 'sens_su']
        :param n_batch: number of batch entries to read. If None, uses n_batch_current.
        :returns: np.ndarray of shape (n_batch, field_dim)
        """
        out_fields = ['x', 'u', 'z', 'pi', 'lam', 'sl', 'su']
        in_fields = ['p']
        sens_fields = ['sens_u', 'sens_x', 'sens_pi', 'sens_lam', 'sens_sl', 'sens_su']
        all_fields = out_fields + in_fields + sens_fields

        if field_ not in all_fields:
            raise ValueError(f"AcadosOcpBatchSolver.get(stage={stage_}, field={field_}): '{field_}' is an invalid argument.\n"
                             f" Possible values are {all_fields}.")

        if not isinstance(stage_, int):
            raise TypeError(f"AcadosOcpBatchSolver.get(stage={stage_}, field={field_}): stage index must be an integer, got type {type(stage_)}.")

        N_horizon = self.__ocp_solvers[0].N
        if stage_ < 0 or stage_ > N_horizon:
            raise ValueError(f"AcadosOcpBatchSolver.get(stage={stage_}, field={field_}): stage index must be in [0, {N_horizon}], got: {stage_}.")

        if stage_ == N_horizon and field_ in ['u', 'pi', 'z', 'sens_u', 'sens_pi']:
            raise KeyError(f"AcadosOcpBatchSolver.get(stage={stage_}, field={field_}): field '{field_}' does not exist at final stage {stage_}.")

        if n_batch is None:
            n_batch = self.n_batch_current
        elif n_batch <= self.n_batch_current:
            self.__n_batch_current = n_batch
        else:
            raise ValueError("You are attempting to get more samples than problem instances initialized so far. "
                             "First initialize enough problem instances by using the setter methods.")

        field_dim = field_.replace('sens_', '') if field_ in sens_fields else field_
        dim = self.__ocp_solvers[0].dims_get(field_dim, stage_)

        out = np.zeros((n_batch, dim), dtype=np.float64, order="C")
        out_data = cast(out.ctypes.data, POINTER(c_double))

        getattr(self.__shared_lib, f"{self.__name}_acados_batch_get")(
            self.__ocp_solvers_pointer, field_.encode('utf-8'), c_int(stage_), out_data,
            c_int(n_batch * dim), c_int(n_batch), c_int(self.__num_threads_in_batch_solve)
        )

        return out


    def set_flat(self, field_: str, value_: np.ndarray) -> None:
        """
        Set concatenation solver initialization for the first `value.shape[0]` solvers.

        :param field_: string in ['x', 'u', 'z', 'pi', 'lam', 'sl', 'su', 'p']
        :param value_: np.array of shape (value.shape[0], n_field_total)
        """

        field = field_.encode('utf-8')
        if field_ not in ['x', 'u', 'z', 'pi', 'lam', 'sl', 'su', 'p']:
            raise ValueError(f'AcadosOcpSolver.get_flat(field={field_}): \'{field_}\' is an invalid argument.')

        dim = self.ocp_solvers[0].get_dim_flat(field_)
        n_batch = value_.shape[0]

        if value_.shape != (n_batch, dim):
            raise ValueError(f'AcadosOcpBatchSolver.set_flat(field={field_}, value): value has wrong shape, expected ({n_batch}, {dim}), got {value_.shape}.')

        if n_batch <= self.n_batch_current:
            self.__n_batch_current = n_batch
        else:
            self.__n_batch_current = n_batch
            self._create_missing_solvers(n_batch)

        value_ = value_.reshape((-1,), order='C')
        N_data = value_.shape[0]

        value_ = value_.astype(float)
        value_data = cast(value_.ctypes.data, POINTER(c_double))

        getattr(self.__shared_lib, f"{self.__name}_acados_batch_set_flat")(self.__ocp_solvers_pointer, field, value_data, N_data, n_batch, self.__num_threads_in_batch_solve)


    def get_flat(self, field_: str, n_batch: Optional[int] = None) -> np.ndarray:
        """
        Get concatenation of all stages of a field of the solver. If `n_batch` is None, the batch size is given by `n_batch_current`.

            :param field: string in ['x', 'u', 'z', 'pi', 'lam', 'sl', 'su', 'p']
            :returns: numpy array of shape (N_batch, n_field_total)
        """
        if field_ not in ['x', 'u', 'z', 'pi', 'lam', 'sl', 'su', 'p']:
            raise ValueError(f'AcadosOcpSolver.get_flat(field={field_}): \'{field_}\' is an invalid argument.')

        if n_batch is None:
            n_batch = self.n_batch_current
        elif n_batch <= self.n_batch_current:
            self.__n_batch_current = n_batch
        else:
            raise ValueError("You are attempting to get more samples than problem instances initialized so far. "
                             "First initialize enough problem instances by using the setter methods.")

        field = field_.encode('utf-8')

        dim = self.ocp_solvers[0].get_dim_flat(field_)

        out = np.ascontiguousarray(np.zeros((n_batch, dim,)), dtype=np.float64)
        out_data = cast(out.ctypes.data, POINTER(c_double))

        getattr(self.__shared_lib, f"{self.__name}_acados_batch_get_flat")(self.__ocp_solvers_pointer, field, out_data, n_batch*dim, n_batch, self.__num_threads_in_batch_solve)

        return out

    @deprecated(version="0.5.4", reason="store_iterate_to_flat_obj is deprecated, use get_flat_iterate instead.")
    def store_iterate_to_flat_obj(self, n_batch: Optional[int] = None) -> AcadosOcpFlattenedBatchIterate:
        return self.get_flat_iterate(n_batch)

    def get_flat_iterate(self, n_batch: Optional[int] = None) -> AcadosOcpFlattenedBatchIterate:
        """
        Returns the current iterate of the first `n_batch` OCP solvers as an AcadosOcpFlattenedBatchIterate.
        The batch size is given by `n_batch_current` if `n_batch` is None.
        """
        if n_batch is None:
            n_batch = self.n_batch_current
        elif n_batch <= self.n_batch_current:
            self.__n_batch_current = n_batch
        else:
            raise ValueError("You are attempting to get more iterates than problem instances initialized so far. "
                             "First initialize enough problem instances by using the setter methods.")

        return AcadosOcpFlattenedBatchIterate(x = self.get_flat("x", n_batch),
                                              u = self.get_flat("u", n_batch),
                                              z = self.get_flat("z", n_batch),
                                              sl = self.get_flat("sl", n_batch),
                                              su = self.get_flat("su", n_batch),
                                              pi = self.get_flat("pi", n_batch),
                                              lam = self.get_flat("lam", n_batch),
                                              N_batch=n_batch)

    @deprecated(version="0.5.4", reason="load_iterate_from_flat_obj is deprecated, use set_iterate instead.")
    def load_iterate_from_flat_obj(self, iterate: AcadosOcpFlattenedBatchIterate) -> None:
        self.set_iterate(iterate)

    def set_iterate(self, iterate: AcadosOcpFlattenedBatchIterate) -> None:
        """
        Loads the provided iterate into the first `n_batch` OCP solvers.
        n_batch is determined by the iterate object.

        Note: The iterate object does not contain the parameters.
        """
        n_batch = iterate.N_batch
        if n_batch <= self.n_batch_current:
            self.__n_batch_current = n_batch
        else:
            self.__n_batch_current = n_batch
            self._create_missing_solvers(n_batch)

        self.set_flat("x", iterate.x)
        self.set_flat("u", iterate.u)
        self.set_flat("z", iterate.z)
        self.set_flat("sl", iterate.sl)
        self.set_flat("su", iterate.su)
        self.set_flat("pi", iterate.pi)
        self.set_flat("lam", iterate.lam)

    def _create_missing_solvers(self, n_batch: int):

        n_batch_max_old = len(self.ocp_solvers)
        n_missing = n_batch - n_batch_max_old

        if n_missing > 0:
            template_solver = self.ocp_solvers[0]
            self.__ocp_solvers.extend([AcadosOcpSolver(template_solver.acados_ocp,
                                                    build=False,
                                                    generate=False,
                                                    verbose=self.verbose if n==0 else False,
                                                    check_reuse_possible=False
                                                    )
                                        for n in range(n_missing)])
            self.__ocp_solvers_pointer = (c_void_p * n_batch)()
            for i in range(len(self.ocp_solvers)):
                self.__ocp_solvers_pointer[i] = self.ocp_solvers[i].capsule

            # Recreate status array
            status_old = self.__status
            self.__status = np.zeros((n_batch,), dtype=np.intc, order="C")
            self.__status_p = cast(self.__status.ctypes.data, POINTER(c_int))
            self.__status[:n_batch_max_old] = status_old

    def constraints_set(self, stage_: int, field_: str, value_: np.ndarray):
        """
        Set numerical data in the constraint module of the solvers.

        :param stage: integer corresponding to shooting node
        :param field: string in ['lbx', 'ubx', 'lbu', 'ubu', 'lg', 'ug', 'lh', 'uh', 'uphi', 'C', 'D']
        :param value: of shape (n_batch, value_dim) for vector-valued fields,
                      or (n_batch, rows, cols) for matrix-valued fields (e.g., 'C', 'D').

        Note: For matrix-valued fields, values are expected in column-major (Fortran) order per sample,
        corresponding to the 'new' API variant of :py:meth:`~acados_template.acados_ocp_solver.AcadosOcpSolver.constraints_set`.
        """
        n_batch = value_.shape[0]

        if value_.ndim < 2:
            raise ValueError(f"Expected batched input with at least two dimensions, got shape {value_.shape}.")

        if n_batch <= self.n_batch_current:
            self.__n_batch_current = n_batch
        else:
            self.__n_batch_current = n_batch
            self._create_missing_solvers(n_batch)

        # For 3D input (matrix-valued fields like C, D), flatten each sample in Fortran (column-major) order.
        # For 2D input (vector-valued fields), data is already flat per sample.
        if value_.ndim == 3:
            data = np.ascontiguousarray(value_.transpose(0, 2, 1).reshape(n_batch, -1), dtype=np.float64)
        else:
            data = np.ascontiguousarray(value_, dtype=np.float64)

        field = field_.encode('utf-8')
        N_data = data.size
        data_p = cast(data.ctypes.data, POINTER(c_double))

        getattr(self.__shared_lib, f"{self.__name}_acados_batch_constraints_set")(
            self.__ocp_solvers_pointer, field, c_int(stage_), data_p, c_int(N_data),
            c_int(n_batch), c_int(self.__num_threads_in_batch_solve)
        )

    def set_p_global_and_precompute_dependencies(self, data_: np.ndarray):
        """
        Sets values of p_global and precomputes all parts of the CasADi graphs of all other functions that only depend on p_global.

        :param data: the global parameters of shape (n_batch, p_global_dim)
        """
        n_batch = data_.shape[0]

        if data_.ndim < 2:
            raise ValueError(f"Expected batched input with at least two dimensions, got shape {data_.shape}.")

        if n_batch <= self.n_batch_current:
            self.__n_batch_current = n_batch
        else:
            self.__n_batch_current = n_batch
            self._create_missing_solvers(n_batch)

        for i, solver in enumerate(self.ocp_solvers[:n_batch]):
            solver.set_p_global_and_precompute_dependencies(data_[i])

    def reset(self, n_batch: Optional[int] = None, reset_qp_solver_mem: bool = True, reset_numerical_values: bool = False, reset_solver_options: bool = False, reset_x_to_x0_bar: bool = False):
        """
        Resets the first n_batch solvers.
        """
        if n_batch is None:
            n_batch = self.n_batch_current
        elif n_batch <= self.n_batch_current:
            self.__n_batch_current = n_batch
        else:
            self.__n_batch_current = n_batch
            self._create_missing_solvers(n_batch)

        for solver in self.ocp_solvers[:n_batch]:
            solver.reset(reset_qp_solver_mem=reset_qp_solver_mem, reset_numerical_values=reset_numerical_values, reset_solver_options=reset_solver_options, reset_x_to_x0_bar=reset_x_to_x0_bar)

    def set(self,  stage_: int, field_: str, data_: np.ndarray):
        """
        Set numerical data inside the solvers using a parallelized batch C call.

        :param stage: integer corresponding to shooting node
        :param field: string in ['x', 'u', 'pi', 'lam', 'z', 'sl', 'su', 'p', 'xdot_guess', 'z_guess', 'sens_x', 'sens_u']
        :param data: the data of shape (n_batch, field_dim).

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
        valid_fields = ['x', 'u', 'pi', 'lam', 'z', 'sl', 'su', 'p', 'xdot_guess', 'z_guess', 'sens_x', 'sens_u']
        if field_ not in valid_fields:
            raise ValueError(f"AcadosOcpBatchSolver.set(): '{field_}' is not a valid argument.\n"
                             f" Possible values are {valid_fields}.")

        n_batch = data_.shape[0]

        if data_.ndim < 2:
            raise ValueError(f"Expected batched input with at least two dimensions, got shape {data_.shape}.")

        if n_batch <= self.n_batch_current:
            self.__n_batch_current = n_batch
        else:
            self.__n_batch_current = n_batch
            self._create_missing_solvers(n_batch)

        field = field_.encode('utf-8')
        data_ = np.ascontiguousarray(data_, dtype=np.float64)
        N_data = data_.size
        data_p = cast(data_.ctypes.data, POINTER(c_double))

        getattr(self.__shared_lib, f"{self.__name}_acados_batch_set")(
            self.__ocp_solvers_pointer, field, c_int(stage_), data_p, c_int(N_data), c_int(n_batch), c_int(self.__num_threads_in_batch_solve)
        )
