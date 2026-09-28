#!/usr/bin/env python3
# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# Licensed under the 2-Clause BSD License.

"""Replace full acados license headers with compact notices."""

from pathlib import Path
import subprocess


REDISTRIBUTION_CLAUSE = (
    "Redistributions in binary form must reproduce" + " the above copyright notice"
)
COPYRIGHT = "Copyright (c) The acados authors."
LICENSE_NOTICE = "Licensed under the 2-Clause BSD License."


def compact_header(comment, newline):
    if comment == "block":
        lines = [
            "/*",
            f" * {COPYRIGHT}",
            " *",
            " * This file is part of acados.",
            " *",
            f" * {LICENSE_NOTICE}",
            " */",
        ]
    elif comment == "plain":
        lines = [
            COPYRIGHT,
            "",
            "This file is part of acados.",
            "",
            LICENSE_NOTICE,
        ]
    else:
        lines = [
            f"{comment} {COPYRIGHT}",
            comment,
            f"{comment} This file is part of acados.",
            comment,
            f"{comment} {LICENSE_NOTICE}",
        ]
    return [line + newline for line in lines]


def shorten_header(path):
    with path.open(encoding="utf-8", newline="") as source:
        lines = source.read().splitlines(keepends=True)
    text = "".join(lines)
    if REDISTRIBUTION_CLAUSE not in text:
        return False

    clause_index = next(
        index for index, line in enumerate(lines) if REDISTRIBUTION_CLAUSE in line
    )
    footer_index = next(
        (
            index
            for index in range(clause_index, len(lines))
            if "POSSIBILITY OF SUCH DAMAGE" in lines[index]
        ),
        None,
    )
    copyright_index = next(
        (
            index
            for index in range(clause_index)
            if "Copyright" in lines[index]
        ),
        None,
    )
    if footer_index is None or copyright_index is None:
        raise ValueError(f"Could not identify license header in {path}")

    copyright_line = lines[copyright_index]
    stripped = copyright_line.lstrip()
    if stripped.startswith("*"):
        comment = "block"
        start_index = copyright_index
        if copyright_index and lines[copyright_index - 1].strip() == "/*":
            start_index -= 1
    elif stripped.startswith("#"):
        comment = "#"
        start_index = copyright_index
    elif stripped.startswith("%"):
        comment = "%"
        start_index = copyright_index
    else:
        comment = "plain"
        start_index = copyright_index

    end_index = footer_index + 1
    if end_index < len(lines):
        following = lines[end_index].strip()
        if (comment == "block" and following == "*/") or following == comment:
            end_index += 1

    newline = next(
        (line[len(line.rstrip("\r\n")) :] for line in lines if line.endswith(("\n", "\r"))),
        "\n",
    )
    replacement = compact_header(comment, newline)
    with path.open("w", encoding="utf-8", newline="") as destination:
        destination.write("".join(lines[:start_index] + replacement + lines[end_index:]))
    return True


def main():
    root = Path(__file__).resolve().parents[2]
    result = subprocess.run(
        ["git", "grep", "-l", "-F", REDISTRIBUTION_CLAUSE],
        cwd=root,
        check=False,
        capture_output=True,
        text=True,
    )
    if result.returncode not in (0, 1):
        raise subprocess.CalledProcessError(result.returncode, result.args)

    changed = 0
    for name in result.stdout.splitlines():
        if name == "LICENSE":
            continue
        changed += shorten_header(root / name)
    print(f"Shortened license headers in {changed} files.")


if __name__ == "__main__":
    main()
