# acados
[![Build Status](https://secure.travis-ci.org/acados/acados.png?branch=master)](http://travis-ci.org/acados/acados)
[![codecov](https://codecov.io/gh/acados/acados/branch/master/graph/badge.svg)](https://codecov.io/gh/acados/acados)

Fast optimal control problem solvers

### Getting started
Make sure `cmake` is installed.
From the `acados` root folder:

    source install_packages.sh  # for Ubuntu
    git submodule update --recursive --init
    mkdir build
    cd build
    cmake -D SWIG_MATLAB=1 -D SWIG_PYTHON=1 ..
    make install
    cd ..
    
    cd examples/python
    python3 ocp_nlp.py
