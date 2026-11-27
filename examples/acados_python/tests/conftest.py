#
# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# Licensed under the 2-Clause BSD License.

import sys
from pathlib import Path


EXAMPLES_DIR = Path(__file__).parents[1]
for directory in (
    EXAMPLES_DIR / "pendulum_on_cart" / "common",
    EXAMPLES_DIR / "chain_mass",
    EXAMPLES_DIR / "getting_started",
    EXAMPLES_DIR / "solution_sensitivities_convex_example",
):
    sys.path.insert(0, str(directory))
