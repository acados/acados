#
# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# Licensed under the 2-Clause BSD License.

import json
import subprocess
import sys
from pathlib import Path

import pytest


TEST_DIR = Path(__file__).parent
PARAMETER_NAMES = (
    "cost_module_terminal",
    "cost_module",
    "hessian_approximation",
    "integrator_type",
    "qp_solver",
    "regularization",
    "solver_type",
)


def _load_combinations(filename: str):
    with (TEST_DIR / "test_data" / filename).open(encoding="utf-8") as file:
        return json.load(file)


CASES = (
    _load_combinations("test_combinations_pendulum_GAUSS_NEWTON.json")
    + _load_combinations("test_combinations_pendulum_EXACT.json")
)


@pytest.mark.parametrize(
    PARAMETER_NAMES,
    [pytest.param(*case, id="-".join(case)) for case in CASES],
)
def test_pendulum_ocp_formulation(
    cost_module_terminal,
    cost_module,
    hessian_approximation,
    integrator_type,
    qp_solver,
    regularization,
    solver_type,
):
    subprocess.run(
        [
            sys.executable,
            str(TEST_DIR / "ocp_setting.py"),
            "--COST_MODULE_N",
            cost_module_terminal,
            "--COST_MODULE",
            cost_module,
            "--HESS_APPROX",
            hessian_approximation,
            "--INTEGRATOR_TYPE",
            integrator_type,
            "--QP_SOLVER",
            qp_solver,
            "--REGULARIZATION",
            regularization,
            "--SOLVER_TYPE",
            solver_type,
        ],
        cwd=TEST_DIR,
        check=True,
    )
