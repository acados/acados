# Interfaces


## C API
The C API of `acados` is an efficient interface to the core functionalities of `acados`.
It provides setters and getters that can be used to interact with the core of `acados` with
negligible computation overhead. In order to learn about the `acados` C API, you
can look at the examples in
[`acados/examples/c/`](https://github.com/acados/acados/tree/master/examples/c).

## MATLAB/Octave
In order to use `acados` from Octave or Matlab, you need to create the `acados` shared libraries  using either the `CMake` or `Make` build system, as described in [here](../installation/index.md).

Additionally, the examples require an installation of `CasADi` to generate the model functions.
Note, that the `getting_started` example offers the option to attempt to automatically download the correct version in the right folder.
Thus, the next subsection can be skipped, but is recommended if there are problems finding `CasADi`.

### Download CasADi:
To create external function for your problem, we suggest to use CasADi from the folder `<acados_root_folder>/external`.
Depending on the environment you want to use to generate CasADi functions from, proceed with the corresponding paragraph (MATLAB, Octave).
For Python, CasADi is automatically downloaded when installing the interface using `pip`.

#### **Matlab**
Put CasADi binaries into `<acados_root_folder>/external/casadi-matlab` :
```
cd external
wget -q -nc --show-progress https://github.com/casadi/casadi/releases/download/3.4.0/casadi-linux-matlabR2014b-v3.4.0.tar.gz
mkdir -p casadi-matlab
tar -xf casadi-linux-matlabR2014b-v3.4.0.tar.gz -C casadi-matlab
cd ..
```

#### **Octave version 4.4 or later**
Put CasADi binaries into `<acados_root_folder>/external/casadi-octave` :
```
cd external
wget -q -nc --show-progress https://github.com/casadi/casadi/releases/download/3.4.5/casadi-linux-octave-4.4.1-v3.4.5.tar.gz
mkdir -p casadi-octave
tar -xf casadi-linux-octave-4.4.1-v3.4.5.tar.gz -C casadi-octave
```

#### **Octave version 4.2 or earliear**
Put CasADi binaries into `<acados_root_folder>/external/casadi-octave` :
```
cd external
wget -q -nc --show-progress https://github.com/casadi/casadi/releases/download/3.4.0/casadi-linux-octave-v3.4.0.tar.gz
mkdir -p casadi-octave
tar -xf casadi-linux-octave-v3.4.0.tar.gz -C casadi-octave
cd ..
```


### Native MEX interface (Rapid Prototyping)
This interface makes a broad set of `acados` functionalities available from Matlab and Octave.
<!-- As of now, this closely tracks the latest developments in the core of acados, e.g.
exact Hessians, adjoint corrections, regularization, etc. -->

To get started with this interface we recommend the examples in `<acados_root>/examples/acados_matlab_octave/getting_started`.

The problem formulation is stated in [this PDF](https://github.com/acados/acados/tree/master/docs/problem_formulation/problem_formulation_ocp_mex.pdf).

A table and explanation of the various options of the interface can be found in a spreadsheet [at this link](https://docs.google.com/spreadsheets/d/1rVRycLnCyaWJLwnV47u30Vokp7vRu68og3OhlDbSjDU/edit?usp=sharing) (thanks to [@EnricaSo](https://github.com/EnricaSo)).

To run the examples, navigate into the selected folder, and there run the command
```
export ACADOS_INSTALL_DIR="<acados_root>"
source env.sh # Which can be found in the folder of one of the examples
```

If `ACADOS_INSTALL_DIR` is not specified, it will be assumed that the examples are run from the sub-folders of the current folder (i.e. `acados` main folder is 2 folders up from the current folder).

Afterwards, launch `Matlab` or Octave from the same shell.

If you want to run the examples in a different folder, please close the current shell and open a new one to repeat the procedure: this ensures the correct setting of the environment variables.

### Templates
There is the option to generate embeddable `C` code from Matlab.
The workflow uses the same templates as the Python interface (see below) and the `Tera` renderer.
After creating an acados solver `ocp`, you can use the routine `ocp.generate_c_code()` to generate `C` code which can be used for embedded applications.
These templates can be found in [`<acados_root>/interfaces/acados_template/acados_template/c_templates_tera`](https://github.com/acados/acados/tree/master/interfaces/acados_template/acados_template/c_templates_tera)

Note: This part of the MATLAB/Octave interface does not yet support all features of the one mentioned before.

### Simulink
The templates mentioned above also contain templated Sfunctions and corresponding make functions for Matlab for both the OCP solver and the acados integrator.
An example can be found in [`<acados_root>/examples/acados_python/getting_started/simulink_example.m`](https://github.com/acados/acados/blob/master/examples/acados_matlab_octave/getting_started/simulink_example.m)
Note that the simulink blocks do not offer the option to change all the numerical values of the OCP description on the fly, as it can be done from Matlab, Octave and Python,
since each kind value that needs to be updated also needs a separate input for the Simulink block.
If you want a more advanced interaction with the `acados` solver via Simulink, feel free to edit the corresponding templates in [`<acados_root>/interfaces/acados_template/acados_template/c_templates_tera`](https://github.com/acados/acados/tree/master/interfaces/acados_template/acados_template/c_templates_tera) to add more inputs or outputs.

## Python (templates)
`acados_template` is a Python package that can be used to specify optimal control problems from Python and to generate self-contained C code that uses the acados solvers to solve them.
The pip package is based on templated code (C files, Header files and Makefiles), which are rendered from Python using the templating engine `Tera`.
The genereated C code can be compiled into a self-contained C library that can be deployed on an embedded system.

### Optimal Control Problem description
The Python interface relies on the same problem formulation as the MATLAB interface [see here](https://github.com/acados/acados/blob/master/docs/problem_formulation/problem_formulation_ocp_mex.pdf).

### Installation
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
    Tipp: you can add these lines to your `.bashrc`/`.zshrc`.

4. Run a Python example to check that everything works.
    We suggest to get started with the example
    `<acados_root>/examples/acados_python/getting_started/minimal_example_ocp.py`.

5. Optional: Can be done automatically through the interface:
    In order to be able to successfully render C code templates, you need to download the `t_renderer` binaries for your platform from <https://github.com/acados/tera_renderer/releases/> and place them in `<acados_root>/bin` (please strip the version and platform from the binaries (e.g.`t_renderer-v0.0.34 -> t_renderer`).
    Notice that you might need to make `t_renderer` executable.
    Run `export ACADOS_SOURCE_DIR=<acados_root>` such that the location of acados will be known to the Python package at run time.

6. Optional: Set `acados_lib_path`, `acados_include_path`.
    If you want the generated Makefile to refer to a specific path (e.g. when cross-compiling or compiling from a location different from the one where you generate the C code), you will have to set these paths accordingly in the generating Python code.

For more information contact `@zanellia`.

### Python API
``` eval_rst
.. automodule:: acados_template.
    :members:
    :private-members:
    :undoc-members:
```

``` eval_rst
.. automodule:: acados_template.acados_ocp
    :members:
    :private-members:
    :exclude-members:
```

``` eval_rst
.. automodule:: acados_template.acados_model
    :members:
    :private-members:
    :exclude-members:
```
<!-- OCP -->


``` eval_rst
.. automodule:: acados_template.acados_ocp_solver
    :members:
    :private-members:
    :exclude-members: make_ocp_dims_consistent, get_ocp_nlp_layout
```

<!-- SIM -->
``` eval_rst
.. automodule:: acados_template.acados_sim
    :members:
    :private-members:
    :exclude-members:
```

``` eval_rst
.. automodule:: acados_template.acados_sim_solver
    :members:
    :private-members:
    :exclude-members: make_ocp_dims_consistent, get_ocp_nlp_layout
```