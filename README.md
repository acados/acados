# acados
[![Travis Status](https://secure.travis-ci.org/acados/acados.png?branch=master)](http://travis-ci.org/acados/acados)
[![Appveyor status](https://ci.appveyor.com/api/projects/status/q0b2nohk476u5clg?svg=true)](https://ci.appveyor.com/project/roversch/acados)
[![codecov](https://codecov.io/gh/acados/acados/branch/master/graph/badge.svg)](https://codecov.io/gh/acados/acados)

Fast and embedded solvers for nonlinear optimal control.

### Installation
If you are on Ubuntu (tested with 16.04), you can run `./install.sh`. You can also
follow the [manual installation instructions](#manual-installation) below.

* soon: binaries for all operating systems available for download (see Releases)

### Getting started

#### MATLAB
First, add CasADi and acados to your MATLAB path. From a MATLAB command window
```
addpath <path_to_acados_root_folder>/external/casadi-matlabR2014b-v3.4.0
addpath <path_to_acados_installation_folder>
% To permanently add these paths:
savepath
```
Run an acados example, from `<path_to_acados_root_folder>/examples/matlab/`:
```
ocp_qp_example.m
```

#### Python

acados only supports `Python3`. Add CasADi and acados to the `PYTHONPATH` environment variable (add those lines to your `.bashrc` or `.zshrc` to set the paths permanently):
```
export PYTHONPATH=<path_to_acados_root_folder>/external/casadi-py35-v3.4.0-64bit:$PYTHONPATH
export PYTHONPATH=<path_to_acados_installation_folder>:$PYTHONPATH
```
To run a Python example from the acados root folder:
```
python examples/python/ocp_nlp.py
```

### Manual Installation

Otherwise, follow the steps below:

1. Install the dependencies:
    ```
    sudo apt-get install libgsl0-dev liblapack-dev libopenblas-dev liboctave-dev libeigen3-dev python3-tk
    sudo apt-get install byacc # for swig
    sudo apt-get install python3-scipy python3-numpy python3-matplotlib
    ```

1. Download CasADi into the `<acados_root_folder>/external` folder:
    ```
    cd external
    wget -q -nc http://files.casadi.org/download/3.4.0/casadi-linux-octave-v3.4.0.tar.gz
    mkdir -p casadi-octave-v3.4.0
    tar -xf casadi-linux-octave-v3.4.0.tar.gz -C casadi-octave-v3.4.0

    wget -q -nc http://files.casadi.org/download/3.4.0/casadi-linux-py35-v3.4.0-64bit.tar.gz
    mkdir -p casadi-py35-v3.4.0-64bit
    tar -xf casadi-linux-py35-v3.4.0-64bit.tar.gz -C casadi-py35-v3.4.0-64bit

    wget -q -nc http://files.casadi.org/download/3.4.0/casadi-linux-matlabR2014b-v3.4.0.tar.gz
    mkdir -p casadi-matlabR2014b-v3.4.0
    tar -xf casadi-linux-matlabR2014b-v3.4.0.tar.gz -C casadi-matlabR2014b-v3.4.0
    cd ..
    ```

1. Initialize all submodules
    ```
    git submodule update --recursive --init
    ```

1. Build and install `swig`. Make sure you don't have an older version installed (e.g. via the package system):
    ```
    cd external/swig
    ./autogen.sh
    ./configure --prefix=$(pwd)/swig_install --enable-silent-rules
    make
    make install > /dev/null    # quiet installation
    export PATH=$(pwd):$PATH    # add swig to PATH
    cd ../.. # back to acados root folder
    ```

1. Build and install `acados`. By default, `acados` is installed in `$HOME/local/lib`. If you want to install `acados` elsewhere, pass `-D ACADOS_INSTALL_DIR=<path_to_acados_installation_folder>` to `cmake` below.
    ```
    mkdir -p build
    cd build
    cmake -D SWIG_MATLAB=1 -D SWIG_PYTHON=1 ..   # set SWIG_MATLAB=0 if you don't have MATLAB installed
    make install
    ```
