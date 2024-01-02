# Python Interface

<!-- ``` eval_rst
.. automodule:: acados_template.
    :members:
    :private-members:
    :undoc-members:
``` -->

`acados_template` is a Python package that can be used to specify optimal control problems from Python and to generate self-contained C code to solve them using `acados`.
The `pip` package is based on templated code (C files, Header files and Makefiles), which are rendered from Python using the templating engine `Tera`.
The genereated C code can be compiled into a self-contained C library that can be deployed on an embedded system.

One can interact with the generated solver using the Python wrapper.
There is a `ctypes` based wrapper which is the default and a `cython` based wrapper which allows for faster interaction with the C code, to allow deployment of the `acados` solver in a Python framework with less overhead.

## Optimal Control Problem description
The Python interface relies on the same problem formulation as the MATLAB interface [see here](https://github.com/acados/acados/blob/master/docs/problem_formulation/problem_formulation_ocp_mex.pdf).
Currently, Python >= 3.8 is required.

## Installation
1. Compile and install `acados` by following the [`CMake` installation instructions](../installation/index.md).

2. Optional: Recommended.
    Create a Python virtual environment using `virtualenv`.
    ```
    virtualenv env --python=/usr/bin/python3
    # Source the environment
    source env/bin/activate
    ```
    Note: There are known path issues with more high level virtual environment managers. Such as `conda`, `miniconda`, `Pycharm`.
    It is not recommended to use them.
    However, if you need to do so and have issues, please have a look in the [`acados` forum](https://discourse.acados.org/).

2. Install `acados_template` Python package:
    ```
    pip install -e <acados_root>/interfaces/acados_template
    ```
    Note: The option `-e` makes the installation editable, so you can seamlessly switch to a later `acados` version and make changes in the Python interface yourself.

3. Add the path to the compiled shared libraries `libacados.so, libblasfeo.so, libhpipm.so` to `LD_LIBRARY_PATH` (default path is `<acados_root/lib>`) by running:
    ```bash
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:"<acados_root>/lib"
    export ACADOS_SOURCE_DIR="<acados_root>"
    ```
    NOTE: On MacOS `DYLD_LIBRARY_PATH` should be used instead of `LD_LIBRARY_PATH`.
    Hint: you can add these lines to your `.bashrc`/`.zshrc`.

4. Run a Python example to check that everything works.
    We suggest to get started with the example
    `<acados_root>/examples/acados_python/getting_started/minimal_example_ocp.py`.
    The example has to be run from the folder in which it is located.

5. Optional: Can be done automatically through the interface:
    In order to be able to successfully render C code templates, you need to download the `t_renderer` binaries for your platform from <https://github.com/acados/tera_renderer/releases/> and place them in `<acados_root>/bin` (please strip the version and platform from the binaries (e.g.`t_renderer-v0.0.34 -> t_renderer`).
    Notice that you might need to make `t_renderer` executable.
    Run `export ACADOS_SOURCE_DIR=<acados_root>` such that the location of acados will be known to the Python package at run time.

6. Optional: Set `acados_lib_path`, `acados_include_path`.
    If you want the generated Makefile to refer to a specific path (e.g. when cross-compiling or compiling from a location different from the one where you generate the C code), you will have to set these paths accordingly in the generating Python code.

### Windows 10+ (WSL) and VSCode

#### Virtual Environment

Follow the installation instructions above, using a virtual environment.

#### Connect to VS-Code

- Install VS-Code on Windows

- Install the [Remote Development Extension](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.vscode-remote-extensionpack).

- Inside the WSL Ubuntu shell run

```bash
code .
```

to start VS-Code in the current folder.

- Select the created virtual environment as Python interpreter.

### Windows

The Python interface of acados can run natively under Windows using the CMake process.
Take a look at the [CMake installation instructions for Workflow with Microsoft Visual C Compiler (MSVC)](../installation/index.html#workflow-with-microsoft-visual-c-compiler-msvc).
We suggest to get started with the example
`<acados_root>/examples/acados_python/getting_started/minimal_example_ocp_cmake.py`. It uses CMake instead of Make by providing an instance of the [CMakeBuilder class](#acados_template.builders.CMakeBuilder) to the constructor of [`AcadosOcpSolver`](#acados_template.acados_ocp_solver.AcadosOcpSolver) or [`AcadosSimSolver`](#acados_template.acados_sim_solver.AcadosSimSolver). Depending on your Visual Studio Version the CMake generator might have to be adapted (see [`CMakeBuilder.generator`](#acados_template.builders.CMakeBuilder.generator)).

## Overview
The following image shows an overview of the available classes in the `acados` Python interface and their dependencies.
``` eval_rst
.. graphviz:: py_acados_classes.dot
    :name: sphinx.ext.graphviz.pyclasses
    :caption: Python API classes overview
    :alt: Overview of acados Python classes
    :align: center
```

## acados OCP solver interface -- `AcadosOcp` and `AcadosOcpSolver`
The class [`AcadosOcp`](#acados_template.acados_ocp.AcadosOcp) can be used to formulate optimal control problems (OCP), for which an acados solver ([`AcadosOcpSolver`](#acados_template.acados_ocp_solver.AcadosOcpSolver)) can be created.

The interface to interact with the acados OCP solver in C is based on [ctypes](https://cython.org/).

Additionally, there is the option to use a Cython based wrapper around the C part.
This offers faster interaction with the solver, because getter and setter calls, which include shape checking are done in compiled C code.
The cython based wrapper is called [`AcadosOcpSolverCython`](#acados_template.acados_ocp_solver.AcadosOcpSolverCython) and was added in [PR #791](https://github.com/acados/acados/pull/791).
`AcadosOcpSolverCython` and `AcadosOcpSolver` offer the same functionality and equivalent behavior is ensured in CI tests.

``` eval_rst
.. automodule:: acados_template.acados_ocp_solver
    :members:
    :private-members:
    :exclude-members: make_ocp_dims_consistent, get_ocp_nlp_layout
```


<!-- ## acados OCP -->
``` eval_rst
.. automodule:: acados_template.acados_ocp
    :members:
    :private-members:
    :exclude-members:
```

``` eval_rst
.. automodule:: acados_template.acados_ocp_cost
    :members:
    :private-members:
    :exclude-members:
```

``` eval_rst
.. automodule:: acados_template.acados_ocp_constraints
    :members:
    :private-members:
    :exclude-members:
```

``` eval_rst
.. automodule:: acados_template.acados_ocp_options
    :members:
    :private-members:
    :exclude-members:
```

``` eval_rst
.. automodule:: acados_template.acados_dims
    :members:
    :private-members:
    :exclude-members:
```



<!-- ## acados model -->
``` eval_rst
.. automodule:: acados_template.acados_model
    :members:
    :private-members:
    :exclude-members:
```



## acados integrator interface -- `AcadosSim` and `AcadosSimSolver`

The class `AcadosSim` can be used to formulate a simulation problem, for which an acados integrator (`AcadosSimSolver`) can be created.
The interface to interact with the acados integrator in C is based on ctypes.

``` eval_rst
.. automodule:: acados_template.acados_sim
    :members:
    :private-members:
    :exclude-members:
```

<!-- ## acados sim solver -->
``` eval_rst
.. automodule:: acados_template.acados_sim_solver
    :members:
    :private-members:
    :exclude-members: make_ocp_dims_consistent, get_ocp_nlp_layout
```

## Builders
The default builder for `acados` is `make` using a `Makefile` generated by `acados`.
If cross-platform compatibility is required `CMake` can be used to build the binaries required by the solvers.
``` eval_rst
.. automodule:: acados_template.builders
    :members:
    :exclude-members:
```


## acados multi-phase OCP
Advanced feature interface to formulate multi-phase OCPs.

Added in
[#1004](https://github.com/acados/acados/pull/1004) and
[#1007](https://github.com/acados/acados/pull/1007).

``` eval_rst
.. automodule:: acados_template.acados_multiphase_ocp
    :members:
    :private-members:
    :exclude-members:
```
