# acados
[![Build Status](https://secure.travis-ci.org/acados/acados.png?branch=master)](http://travis-ci.org/acados/acados)
[![codecov](https://codecov.io/gh/acados/acados/branch/master/graph/badge.svg)](https://codecov.io/gh/acados/acados)

Fast optimal control problem solvers

### Installation
If you are on Ubuntu (tested with 16.04), you can use the `install.sh` script that does everything for you.

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
    wget -q -nc http://files.casadi.org/3.1.1/linux/casadi-octave-v3.1.1.tar.gz
    mkdir -p casadi-octave-v3.1.1
    tar -xf casadi-octave-v3.1.1.tar.gz -C casadi-octave-v3.1.1

    wget -q -nc http://files.casadi.org/3.1.1/linux/casadi-py35-np1.9.1-v3.1.1.tar.gz
    mkdir -p casadi-py35-np1.9.1-v3.1.1
    tar -xf casadi-py35-np1.9.1-v3.1.1.tar.gz -C casadi-py35-np1.9.1-v3.1.1

    wget -q -nc https://sourceforge.net/projects/casadi/files/CasADi/3.1.1/linux/casadi-matlabR2014b-v3.1.1.tar.gz
    mkdir -p casadi-matlabR2014b-v3.1.1
    tar -xf casadi-matlabR2014b-v3.1.1.tar.gz -C casadi-matlabR2014b-v3.1.1
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
    make install > /dev/null # quiet installation
    export PATH=$(pwd):$PATH # add to path environment variable to make it visible for the acados build process
    cd ../.. # back to acados root folder
    ```

1. Build `acados`
    ```
    mkdir -p build
    cd build
    cmake -D SWIG_MATLAB=1 -D SWIG_PYTHON=1 .. # set SWIG_MATLAB=0 if you don't have matlab installed
    make install # installs into ~/local/lib
    ```

### Getting started

To use acados, you have to add acados and Casadi paths to the respective environment variables (add those lines to your .bashrc / .zshrc to set the paths permanently):
```
export MATLABPATH=<path_to_acados_root_folder>/external/casadi-matlabR2014b-v3.1.1:$MATLABPATH
export PYTHONPATH=<path_to_acados_root_folder>/external/casadi-py35-np1.9.1-v3.1.1:$PYTHONPATH
export PYTHONPATH=<path_to_acados_installation_folder>:$PYTHONPATH # usually ~/local/lib
```

Now you can run the python example:
```
python3 examples/python/ocp_nlp.py
```
