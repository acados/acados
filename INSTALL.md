# Installation Guide

## Quick Start with pip

The simplest way to install acados is using pip:

```bash
pip install acados
```

This will:
- Download the acados source code
- Build the C libraries (requires cmake and a C compiler)
- Install the Python interface

### Setting Library Path

After installation, you need to set the library path so the system can find the acados shared libraries:

**Linux:**
```bash
export LD_LIBRARY_PATH=$(python -c "from acados_template import get_lib_path; print(get_lib_path())"):$LD_LIBRARY_PATH
```

**macOS:**
```bash
export DYLD_LIBRARY_PATH=$(python -c "from acados_template import get_lib_path; print(get_lib_path())"):$DYLD_LIBRARY_PATH
```

**Add to your shell profile** (e.g., `~/.bashrc` or `~/.zshrc`) to make it permanent:
```bash
# Add to ~/.bashrc or ~/.zshrc
export LD_LIBRARY_PATH=$(python -c "from acados_template import get_lib_path; print(get_lib_path())"):$LD_LIBRARY_PATH
```

### Alternative: Set in Python

You can also set the library path programmatically in your Python scripts:

```python
import os
from acados_template import get_lib_path

# Set library path
lib_path = get_lib_path()
os.environ['LD_LIBRARY_PATH'] = lib_path + os.pathsep + os.environ.get('LD_LIBRARY_PATH', '')
```

Note: Setting `LD_LIBRARY_PATH` in Python may not work if the libraries are already loaded. It's better to set it before running Python.

## Prerequisites

- Python >= 3.8
- CMake >= 3.7
- C compiler (gcc, clang, or MSVC)
- Build tools (make on Linux/macOS)

## Development Installation

For development, you can install in editable mode:

```bash
git clone https://github.com/acados/acados.git
cd acados
git submodule update --init --recursive
pip install -e .
```

This will build the C libraries and install the Python package in development mode.

## Traditional Installation

If you prefer the traditional installation method, see the detailed instructions in [docs/installation/index.md](docs/installation/index.md).

## Verifying Installation

Test your installation:

```bash
python -c "from acados_template import AcadosOcp; print('acados imported successfully')"
```

## Troubleshooting

### "libacados.so: cannot open shared object file"

This means the library path is not set. Follow the "Setting Library Path" instructions above.

### CMake not found

Install CMake:
- Ubuntu/Debian: `sudo apt-get install cmake`
- macOS: `brew install cmake`
- Windows: Download from https://cmake.org/download/

### Build fails

Make sure you have:
- A C compiler installed
- Git submodules initialized: `git submodule update --init --recursive`
- Sufficient disk space (the build requires ~500MB)

For more help, visit our forum at https://discourse.acados.org/
