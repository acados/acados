#
# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# Licensed under the 2-Clause BSD License.

from pathlib import Path


if __name__ == "__main__":
    import pytest

    raise SystemExit(
        pytest.main([str(Path(__file__).with_name("test_ocp_formulations.py"))])
    )

