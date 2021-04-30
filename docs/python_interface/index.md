# Python Interface

<!-- ``` eval_rst
.. automodule:: acados_template.
    :members:
    :private-members:
    :undoc-members:
``` -->

`acados_template` is a Python package that can be used to specify optimal control problems from Python and to generate self-contained C code to solve them using `acados`.
The pip package is based on templated code (C files, Header files and Makefiles), which are rendered from Python using the templating engine `Tera`.
The genereated C code can be compiled into a self-contained C library that can be deployed on an embedded system.

## Optimal Control Problem description
The Python interface relies on the same problem formulation as the MATLAB interface [see here](https://github.com/acados/acados/blob/master/docs/problem_formulation/problem_formulation_ocp_mex.pdf).

## Installation
1. Compile and install `acados` by following the [`CMake` installation instructions](../installation/index.md).

2. Optional: Recommended.
    Create a Python virtual environment using `virtualenv`.
    ```
    virtualenv env --python=/usr/bin/python3.7
    # Source the environment
    source env/bin/activate
    ```
    Note: There are known path issues with more high level virtual environment managers. Such as `conda`, `miniconda`, `Pycharm`.
    It is not recommended to use them.
    However, if you need to do so and have issues, please have a look in the [`acados` forum](discourse.acados.org/).

2. Install `acados_template` Python package:
    ```
    pip3 install -e <acados_root>/interfaces/acados_template
    ```
    Note: If you are working with a virtual Python environment, use the `pip` corresponding to this Python environment instead of `pip3`.
    Note: The option `-e` makes the installation editable, so you can seeminglessly switch to a later `acados` version and make changes in the Python interface yourself.

3. Add the path to the compiled shared libraries `libacados.so, libblasfeo.so, libhpipm.so` to `LD_LIBRARY_PATH` (default path is `<acados_root/lib>`) by running:
    ```bash
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:"<acados_root>/lib"
    export ACADOS_SOURCE_DIR="<acados_root>"
    ```
    Hint: you can add these lines to your `.bashrc`/`.zshrc`.

4. Run a Python example to check that everything works.
    We suggest to get started with the example
    `<acados_root>/examples/acados_python/getting_started/minimal_example_ocp.py`.

5. Optional: Can be done automatically through the interface:
    In order to be able to successfully render C code templates, you need to download the `t_renderer` binaries for your platform from <https://github.com/acados/tera_renderer/releases/> and place them in `<acados_root>/bin` (please strip the version and platform from the binaries (e.g.`t_renderer-v0.0.34 -> t_renderer`).
    Notice that you might need to make `t_renderer` executable.
    Run `export ACADOS_SOURCE_DIR=<acados_root>` such that the location of acados will be known to the Python package at run time.

6. Optional: Set `acados_lib_path`, `acados_include_path`.
    If you want the generated Makefile to refer to a specific path (e.g. when cross-compiling or compiling from a location different from the one where you generate the C code), you will have to set these paths accordingly in the generating Python code.


## Overview
The following image shows an overview of the available classes in the `acados` Python interface and their dependencies.
``` eval_rst
.. graphviz:: py_acados_classes.dot
    :name: sphinx.ext.graphviz.pyclasses
    :caption: Python API classes overview
    :alt: Overview of acados Python classes
    :align: center
```

## acados OCP solver
The class `AcadosOcp` can be used to formulate optimal control problems, for which an acados solver (`AcadosOcpSolver`) can be created.

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


<!-- ## acados model -->
``` eval_rst
.. automodule:: acados_template.acados_model
    :members:
    :private-members:
    :exclude-members:
```



## acados integrator interface
The class `AcadosSim` can be used to formulate a simulation problem, for which an acados integrator (`AcadosSimSolver`) can be created.
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
