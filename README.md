# acados
<!-- [![Travis Status](https://secure.travis-ci.org/acados/acados.png?branch=master)](http://travis-ci.org/acados/acados) -->
[![Appveyor status](https://ci.appveyor.com/api/projects/status/q0b2nohk476u5clg?svg=true)](https://ci.appveyor.com/project/roversch/acados)
![Github actions full build workflow](https://github.com/acados/acados/actions/workflows/full_build.yml/badge.svg)
<!-- [![codecov](https://codecov.io/gh/acados/acados/branch/master/graph/badge.svg)](https://codecov.io/gh/acados/acados) -->

Fast and embedded solvers for nonlinear optimal control.

## General
- `acados` offers interfaces to the programming languages `C`, `Python`, `MATLAB` and `Octave`
- Documentation can be found on [docs.acados.org](https://docs.acados.org/)
- Forum: If you have any acados-related question, feel free to post on our forum [discourse.acados.org](https://discourse.acados.org/).
- Citing acados: references can be found [here](https://docs.acados.org/citing).

## Installation

1. Initialize all submodules
    ```
    git submodule update --recursive --init
    ```

1. Build and install `acados`.
Both a CMake and a Makefile based build system is supported at the moment.
Please choose one and proceed with the corresponding paragraph.

    ### **CMake**
    Set the `BLASFEO_TARGET` in `<acados_root>/CMakeLists.txt`.
    For a list of supported targets, we refer to https://github.com/giaf/blasfeo/blob/master/README.md .
    Install acados as follows
    ```
    mkdir -p build
    cd build
    cmake .. # with optional arguments e.g. -DACADOS_WITH_OSQP=OFF/ON -DACADOS_INSTALL_DIR=<path_to_acados_installation_folder>
    make install
    ```

    ### **Make**
    Set the `BLASFEO_TARGET` in `<acados_root>/Makefile.rule`.
    For a list of supported targets, we refer to https://github.com/giaf/blasfeo/blob/master/README.md .
    Install acados as follows
    ```
    make shared_library
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:<path_to_acados_folder>/lib
    make examples_c
    make run_examples_c
    ```

## `acados` interfaces
`acados` written in `C` and offers interfaces to the programming languages `C`, `Python`, `MATLAB` and `Octave`.

### `C` interface
- This includes all functionality in `<acados_root>/interfaces/acados_c`.
- Documentation can be found here: [C API](https://docs.acados.org/c_interface/index.html)

### `MATLAB` and `Octave` interface
- The interface is written in `MEX`
- It deals with the problem formulation described in [this PDF](https://github.com/acados/acados/tree/master/docs/problem_formulation/problem_formulation_ocp_mex.pdf)
- To get started
    - install `acados` as described above
    - Octave users: provide a [`CasADi version`](https://web.casadi.org/get/) in `<acados_root>/external/casadi_octave/`.
      We recommend version 3.4.5. For Matlab, a `CasADi` version should be provided in `<acados_root>/external/casadi_matlab/`, however, this is done automatically when using the examples in the `getting_started` folder.
      On Linux machines with Octave 4.4.1. or later this can be done as follows:
      ```
        cd external
        wget -q -nc --show-progress https://github.com/casadi/casadi/releases/download/3.4.5/casadi-linux-octave-4.4.1-v3.4.5.tar.gz
        mkdir -p casadi-octave
        tar -xf casadi-linux-octave-4.4.1-v3.4.5.tar.gz -C casadi-octave
        ```
    - Linux and MacOS users:
        - in a terminal navigate to `<acados_root>/examples/acados_matlab_octave/getting_started`
        - run `source env.sh` and start Matlab/Octave
    - Windows users:
        - start Matlab
        - run `acados_examples_env.m` in `<acados_root>/examples/acados_matlab_octave`
    - enjoy the examples in `<acados_root>/examples/acados_matlab_octave/getting_started`
- More documentation can be found on [docs.acados.org/interfaces/](https://docs.acados.org/interfaces/)


### `Python` interface
- The interface is based on templated `C` code, header files and a Makefile, which are rendered with the templating engine [`Tera`](https://docs.rs/tera/0.10.0/tera/struct.Tera.html).
- The generated `C` code can be easily deployed on embedded platforms.
- The same problem formulation as for the MATLAB interface is used (see [here](https://github.com/acados/acados/blob/master/docs/problem_formulation/problem_formulation_ocp_mex.pdf)).
- The `ctypes` package is used to interact with the rendered Code.
- `CasADi` is required and is downloaded automatically, [when installing the Python interface using `pip`](https://docs.acados.org/interfaces/index.html#python-templates)
- More documentation can be found on [docs.acados.org/interfaces/](https://docs.acados.org/interfaces/)
