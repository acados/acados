#
# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# Licensed under the 2-Clause BSD License.

import sys
from pathlib import Path


def add_example_paths():
    examples_dir = Path(__file__).parents[1]
    for directory in (
        examples_dir / "pendulum_on_cart" / "common",
        examples_dir / "chain_mass",
        examples_dir / "getting_started",
        examples_dir / "solution_sensitivities_convex_example",
    ):
        path = str(directory)
        if path not in sys.path:
            sys.path.insert(0, path)
