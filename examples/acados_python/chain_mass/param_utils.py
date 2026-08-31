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

import numpy as np
import casadi as ca
from typing import Tuple, Union


class ParamLayout:
    """Describes how named (and possibly repeated) parameter blocks are laid out inside
    one flat vector.

    Each entry is a dict with:
        name:   str
        shape:  int (vector of that length) or tuple(int, int) (matrix, column-major)
        repeat: int, how many times this named block is repeated (default 1)
    """

    def __init__(self, entries: list[dict]):
        self._offsets: dict[str, list[Tuple[int, int, Union[int, Tuple[int, int]]]]] = {}
        offset = 0
        for e in entries:
            name = e["name"]
            shape = e["shape"]
            repeat = e.get("repeat", 1)
            size = shape[0] * shape[1] if isinstance(shape, tuple) else shape
            blocks = []
            for _ in range(repeat):
                blocks.append((offset, size, shape))
                offset += size
            self._offsets[name] = blocks
        self.size = offset

    def block(self, name: str, i: int = 0) -> Tuple[int, int, Union[int, Tuple[int, int]]]:
        return self._offsets[name][i]

    def indices(self, name: str, i: int = 0) -> list[int]:
        start, size, _ = self.block(name, i)
        return list(range(start, start + size))

    def label_indices(self, label: str) -> list[int]:
        """Parse a label such as 'Q', 'm_3', or 'C_2_0' (name[_repeat_index[_component]])
        and return the corresponding flat indices. This replaces the old
        `find_idx_for_labels(struct.cat, label)` helper, which parsed the printed
        string representation of a casadi struct -- here we just read the layout.
        """
        parts = label.split("_")
        name = parts[0]
        if name not in self._offsets:
            raise KeyError(f"Unknown parameter name '{name}' in label '{label}'")

        if len(parts) == 1:
            # whole entry (across all repeats)
            idxs: list[int] = []
            for start, size, _ in self._offsets[name]:
                idxs.extend(range(start, start + size))
            return idxs
        elif len(parts) == 2:
            i = int(parts[1])
            return self.indices(name, i)
        elif len(parts) == 3:
            i, j = int(parts[1]), int(parts[2])
            start, _, _ = self.block(name, i)
            return [start + j]
        else:
            raise ValueError(f"Cannot parse parameter label '{label}'")


class ParamVector:
    """Wraps a flat vector (SX for symbolic use, numpy array for numeric use) together
    with a `ParamLayout`, giving dict/struct-like named access:

        pv["D", i]          # block i of entry "D" (vector or reshaped matrix)
        pv["D", i] = value   # set block i of entry "D" (numeric backing only)
        pv["Q"]              # whole entry "Q" (all repeats concatenated)
        pv.cat                # underlying flat vector
    """

    def __init__(self, layout: ParamLayout, vec):
        self.layout = layout
        self.cat = vec  # SX column vector (size, 1) or numpy array (size, 1)

    def _split_key(self, key):
        return key if isinstance(key, tuple) else (key, None)

    def __getitem__(self, key):
        name, i = self._split_key(key)
        blocks = self.layout._offsets[name]

        if i is None:
            start = blocks[0][0]
            end = blocks[-1][0] + blocks[-1][1]
            return self.cat[start:end, :]

        start, size, shape = blocks[i]
        s = self.cat[start:start + size, :]
        if isinstance(shape, tuple):
            if isinstance(self.cat, np.ndarray):
                s = np.reshape(s, shape, order="F")
            else:
                s = ca.reshape(s, shape[0], shape[1])
        return s

    def __setitem__(self, key, value):
        if not isinstance(self.cat, np.ndarray):
            raise TypeError("Symbolic ParamVector entries cannot be assigned numeric values.")

        name, i = self._split_key(key)
        blocks = self.layout._offsets[name]

        if i is None:
            start = blocks[0][0]
            end = blocks[-1][0] + blocks[-1][1]
        else:
            start, size, _ = blocks[i]
            end = start + size

        arr = np.array(value, dtype=float).flatten(order="F").reshape(-1, 1)
        self.cat[start:end, :] = arr