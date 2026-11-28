#
# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# Licensed under the 2-Clause BSD License.
#
"""
fatrop backend for acados OCPs.

`AcadosFatropOcpSolver` is an `AcadosOcpSolver` whose `solve` runs fatrop's interior point method
instead of acados' NLP solver. The OCP is set up, evaluated and accessed through acados as usual:
fatrop's NLP interface is implemented on top of acados' NLP modules
(c_templates_tera/fatrop/acados_fatrop_adapter.cpp), so dynamics, cost and constraints are
evaluated by the same code as in acados' SQP.
"""

import ctypes
import os
import shutil
import subprocess
from typing import Optional

from .acados_ocp import AcadosOcp
from .acados_ocp_solver import AcadosOcpSolver
from .utils import fatrop_return_flag_to_status, get_shared_lib, get_shared_lib_ext, get_shared_lib_prefix

_SOURCE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'c_templates_tera', 'fatrop')
# acados option -> fatrop option
_OPTIONS = {'nlp_solver_max_iter': 'max_iter', 'nlp_solver_tol_stat': 'tolerance',
            'nlp_solver_tol_eq': 'constr_viol_tol', 'nlp_solver_tol_ineq': 'constr_viol_tol',
            'nlp_solver_tol_comp': 'compl_inf_tol', 'print_level': 'print_level'}
_STATS = ('nlp_iter', 'time_tot', 'status', 'fatrop_status')


class AcadosFatropOcpSolver(AcadosOcpSolver):
    """
    Solve an `AcadosOcp` with fatrop.

    Everything except `solve` is inherited from `AcadosOcpSolver`: set the initial guess, parameters,
    cost and constraint data with `set`, `cost_set`, `constraints_set`, read the solution with `get`
    and `get_cost`. `solve` starts from the iterate in the solver and stores the result in it.

    Constraint rows (bounds, g and h) that acados marks as equalities (idxe, e.g. the initial state
    via idxbxe_0) are equality constraints for fatrop. Other rows without finite bound when the
    solver is created are dropped, all others are inequalities. This structure is fixed: `solve`
    raises if a row changed between bounded and unbounded.

    Not supported: soft constraints (slacks), algebraic variables z, BGP constraints, multi-phase
    OCPs, fixed_hess, nlp_solver_type other than SQP / SQP_RTI.
    Ignored (fatrop uses its own algorithm): the QP solver, globalization, regularization and
    Levenberg-Marquardt options of acados.

    Options: nlp_solver_max_iter -> max_iter - 1 (fatrop stops after max_iter + 1 iterations, so at
    most nlp_solver_max_iter iterations are done), nlp_solver_tol_stat -> tolerance (fatrop's overall
    scaled NLP error), min(nlp_solver_tol_eq, nlp_solver_tol_ineq) -> constr_viol_tol,
    nlp_solver_tol_comp -> compl_inf_tol, print_level > 0 -> fatrop's iteration output.

    Status: fatrop's return flag is mapped by :py:func:`acados_template.utils.fatrop_return_flag_to_status`;
    the raw flag is available as get_stats('fatrop_status').

    fatrop has to be built against the BLASFEO of acados, see the documentation.

    :param ocp: the OCP
    :param fatrop_dir: directory containing fatropConfig.cmake (default: env FATROP_DIR, else
        CMake's default search path)
    :param fatrop_options: dict of fatrop options, applied after the acados options
    :param build: build the acados solver and the fatrop backend; other arguments as in `AcadosOcpSolver`
    """

    def __init__(self, ocp: AcadosOcp, fatrop_dir: Optional[str] = None,
                 fatrop_options: Optional[dict] = None, build: bool = True, generate: bool = True,
                 verbose: bool = True):
        if not isinstance(ocp, AcadosOcp):
            raise TypeError('ocp should be of type AcadosOcp.')
        self._fatrop_handle = None
        super().__init__(ocp, build=build, generate=generate, verbose=verbose)
        self._check_supported()

        code_dir = os.path.join(os.path.abspath(ocp.code_gen_options.code_export_directory), 'fatrop')
        lib_path = os.path.join(code_dir, f'{get_shared_lib_prefix()}acados_fatrop{get_shared_lib_ext()}')
        if build or not os.path.exists(lib_path):
            self._build_fatrop(code_dir, fatrop_dir)
        self._load_fatrop(lib_path)

        opts = ocp.solver_options
        for name in ('nlp_solver_max_iter', 'nlp_solver_tol_stat', 'nlp_solver_tol_comp', 'print_level'):
            self.options_set(name, getattr(opts, name))
        self.options_set('nlp_solver_tol_eq', min(opts.nlp_solver_tol_eq, opts.nlp_solver_tol_ineq))
        for name, value in (fatrop_options or {}).items():
            self.options_set(name, value)

    # ------------------------------------------------------------------ setup
    def _check_supported(self):
        ocp = self.ocp
        dims, opts, c = ocp.dims, ocp.solver_options, ocp.constraints
        if opts.N_horizon == 0:
            raise NotImplementedError('AcadosFatropOcpSolver: N_horizon = 0 is not supported.')
        if opts.nlp_solver_type not in ('SQP', 'SQP_RTI'):
            raise NotImplementedError('AcadosFatropOcpSolver: nlp_solver_type must be SQP or SQP_RTI.')
        if opts.fixed_hess:
            raise NotImplementedError('AcadosFatropOcpSolver: fixed_hess is not supported.')
        if dims.nz > 0:
            raise NotImplementedError('AcadosFatropOcpSolver: algebraic variables z are not supported.')
        if dims.ns_0 + dims.ns + dims.ns_e > 0:
            raise NotImplementedError('AcadosFatropOcpSolver: soft constraints (slacks) are not supported.')
        for attr in ('constr_type_0', 'constr_type', 'constr_type_e'):
            if getattr(c, attr) != 'BGH':
                raise NotImplementedError(f'AcadosFatropOcpSolver: {attr} = {getattr(c, attr)} not supported, use BGH.')

    def _build_fatrop(self, code_dir: str, fatrop_dir: Optional[str]):
        """CMake build; the number of parallel jobs follows CMAKE_BUILD_PARALLEL_LEVEL."""
        os.makedirs(code_dir, exist_ok=True)
        for f in ('acados_fatrop_adapter.cpp', 'CMakeLists.txt'):
            shutil.copy(os.path.join(_SOURCE_DIR, f), code_dir)
        cgo = self.ocp.code_gen_options
        build_dir = os.path.join(code_dir, 'build')
        cfg = ['cmake', '-S', code_dir, '-B', build_dir, '-DCMAKE_BUILD_TYPE=Release',
               f'-DACADOS_INCLUDE_PATH={cgo.acados_include_path}', f'-DACADOS_LIB_PATH={cgo.acados_lib_path}']
        fatrop_dir = fatrop_dir or os.environ.get('FATROP_DIR')
        if fatrop_dir:
            cfg.append(f'-Dfatrop_DIR={fatrop_dir}')
        for cmd in (cfg, ['cmake', '--build', build_dir, '--config', 'Release']):
            res = subprocess.run(cmd, capture_output=True, text=True)
            if res.returncode != 0:
                raise RuntimeError(f'AcadosFatropOcpSolver: build failed:\n{res.stdout}\n{res.stderr}')

    def _load_fatrop(self, lib_path: str):
        lib = get_shared_lib(lib_path)
        vp, cp, dp = ctypes.c_void_p, ctypes.c_char_p, ctypes.POINTER(ctypes.c_double)
        lib.acados_fatrop_create.argtypes = [vp, vp, vp]
        lib.acados_fatrop_create.restype = vp
        lib.acados_fatrop_free.argtypes = [vp]
        lib.acados_fatrop_option_set.argtypes = [vp, cp, ctypes.c_double, cp]
        lib.acados_fatrop_solve.argtypes = [vp]
        lib.acados_fatrop_get_stat.argtypes = [vp, cp, dp]
        self._fatrop_lib = lib
        self._fatrop_handle = lib.acados_fatrop_create(self.nlp_solver, self.nlp_in, self.nlp_out)
        if not self._fatrop_handle:
            raise RuntimeError('AcadosFatropOcpSolver: creating the fatrop solver failed.')

    def __del__(self):
        # the fatrop solver refers to the acados solver memory, free it first
        if getattr(self, '_fatrop_handle', None):
            self._fatrop_lib.acados_fatrop_free(self._fatrop_handle)
            self._fatrop_handle = None
        super().__del__()

    # ------------------------------------------------------------------ public API
    def options_set(self, field: str, value):
        """
        Set a fatrop option: the acados options nlp_solver_max_iter, nlp_solver_tol_stat,
        nlp_solver_tol_eq, nlp_solver_tol_ineq, nlp_solver_tol_comp and print_level are mapped to
        fatrop (see class docstring); any other name is passed to fatrop unchanged. Numbers are
        converted to the type of the fatrop option.
        """
        is_str = isinstance(value, str)
        if not is_str and not isinstance(value, (bool, int, float)) and not hasattr(value, '__float__'):
            raise TypeError(f'AcadosFatropOcpSolver: option {field} must be a number, bool or string, got {value!r}.')
        if field == 'print_level' and not is_str:
            value = 5 if value > 0 else 0
        if field == 'nlp_solver_max_iter' and not is_str:
            value = value - 1
        name = _OPTIONS.get(field, field)
        ret = self._fatrop_lib.acados_fatrop_option_set(self._fatrop_handle, name.encode(),
                                                        0. if is_str else float(value),
                                                        value.encode() if is_str else None)
        if ret == 2:
            raise ValueError(f'AcadosFatropOcpSolver: unknown option {field}.')
        if ret == 3:
            raise TypeError(f'AcadosFatropOcpSolver: wrong type for option {field}: {value!r}.')
        if ret:
            raise ValueError(f'AcadosFatropOcpSolver: could not set option {field} to {value!r}.')

    def solve(self) -> int:
        """Solve the OCP with fatrop; returns the acados status, see the class docstring."""
        flag = int(self._fatrop_lib.acados_fatrop_solve(self._fatrop_handle))
        self._status = fatrop_return_flag_to_status(flag)
        if flag == -2:
            raise ValueError('AcadosFatropOcpSolver: the constraint structure changed: inequality rows without '
                             'finite bound at creation must stay unbounded, all others must keep a finite bound.')
        return self._status

    def get_stats(self, field: str):
        """Fields: nlp_iter, time_tot [s], status (acados), fatrop_status (fatrop's return flag)."""
        if field not in _STATS:
            raise NotImplementedError(f'AcadosFatropOcpSolver.get_stats: field {field} not supported, use one of {_STATS}.')
        if field == 'status':
            return self._status
        out = ctypes.c_double()
        name = 'status' if field == 'fatrop_status' else field
        if self._fatrop_lib.acados_fatrop_get_stat(self._fatrop_handle, name.encode(), ctypes.byref(out)):
            raise RuntimeError(f'AcadosFatropOcpSolver: getting {field} failed.')
        return out.value if field == 'time_tot' else int(out.value)

    def print_statistics(self):
        """Print a summary of the last solve."""
        print(f"fatrop: status {self._status} (fatrop flag {self.get_stats('fatrop_status')}), "
              f"{self.get_stats('nlp_iter')} iterations, {1e3 * self.get_stats('time_tot'):.3f} ms, "
              f"cost {self.get_cost():.6e}")
