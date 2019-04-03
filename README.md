# acados
[![Travis Status](https://secure.travis-ci.org/acados/acados.png?branch=master)](http://travis-ci.org/acados/acados)
[![Appveyor status](https://ci.appveyor.com/api/projects/status/q0b2nohk476u5clg?svg=true)](https://ci.appveyor.com/project/roversch/acados)
[![codecov](https://codecov.io/gh/acados/acados/branch/master/graph/badge.svg)](https://codecov.io/gh/acados/acados)

Fast and embedded solvers for nonlinear optimal control.

### Optional requirements
Some functionalities in acados require CasADi (version 3.4.0) to be installed on your system.
To install CasADi, you can follow the installation instructions [here](https://github.com/casadi/casadi/wiki/InstallationInstructions)

### Installation
Both a CMake and a Makefile based build system are supported at the moment.

1. Initialize all submodules
    ```
    git submodule update --recursive --init
    ```

1. Download CasADi into the `<acados_root_folder>/external` folder:
    ```
    cd external
    ```
    and, depending on your preferred CasADi interface (Python, MATLAB, Octave):

    ```
    wget -q -nc --show-progress https://github.com/casadi/casadi/releases/download/3.4.0/casadi-linux-py35-v3.4.0-64bit.tar.gz
    mkdir -p casadi-py35-v3.4.0-64bit
    tar -xf casadi-linux-py35-v3.4.0-64bit.tar.gz -C casadi-py35-v3.4.0-64bit
    cd ..
    ```

    ```
    wget -q -nc --show-progress https://github.com/casadi/casadi/releases/download/3.4.0/casadi-linux-matlabR2014b-v3.4.0.tar.gz
    mkdir -p casadi-matlabR2014b-v3.4.0
    tar -xf casadi-linux-matlabR2014b-v3.4.0.tar.gz -C casadi-matlabR2014b-v3.4.0
    cd ..
    ```

    ```
    wget -q -nc --show-progress https://github.com/casadi/casadi/releases/download/3.4.0/casadi-linux-octave-v3.4.0.tar.gz
    mkdir -p casadi-octave-v3.4.0
    tar -xf casadi-linux-octave-v3.4.0.tar.gz -C casadi-octave-v3.4.0
    cd ..
    ```

1. Build and install `acados`.
    When using the CMake-based build sytem:
    ```
    mkdir -p build
    cd build
    cmake .. (with optional arguments e.g. -DACADOS_WITH_OSQP=OFF/ON -DACADOS_INSTALL_DIR=<path_to_acados_installation_folder>)
    make install
    ```

    When using the Makefile-based build sytem:
    ```
    make acados_shared
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:<path_to_acados_folder>/lib
    make examples_c
    make run_examples_c
    ```

* soon: binaries for all operating systems available for download (see Releases)
