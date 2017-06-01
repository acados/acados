#!/bin/bash -xe

# Install dependencies
sudo apt-get install libgsl0-dev liblapack-dev libopenblas-dev liboctave-dev libeigen3-dev python3-tk
sudo apt-get install byacc # swig
sudo apt-get install python3-scipy python3-numpy python3-matplotlib

# Get CasADi for octave, python and matlab
pushd external
wget -q -nc http://files.casadi.org/3.1.1/linux/casadi-octave-v3.1.1.tar.gz
mkdir -p casadi-octave-v3.1.1
tar -xf casadi-octave-v3.1.1.tar.gz -C casadi-octave-v3.1.1

wget -q -nc http://files.casadi.org/3.1.1/linux/casadi-py35-np1.9.1-v3.1.1.tar.gz
mkdir -p casadi-py35-np1.9.1-v3.1.1
tar -xf casadi-py35-np1.9.1-v3.1.1.tar.gz -C casadi-py35-np1.9.1-v3.1.1
export PYTHONPATH=$(pwd)/casadi-py35-np1.9.1-v3.1.1:$PYTHONPATH

wget -q -nc https://sourceforge.net/projects/casadi/files/CasADi/3.1.1/linux/casadi-matlabR2014b-v3.1.1.tar.gz
mkdir -p casadi-matlabR2014b-v3.1.1
tar -xf casadi-matlabR2014b-v3.1.1.tar.gz -C casadi-matlabR2014b-v3.1.1
export MATLABPATH=$(pwd)/casadi-matlabR2014b-v3.1.1:$MATLABPATH


# Get all git submodules
git submodule update --recursive --init


# Install swig
pushd swig
./autogen.sh
./configure --prefix=$(pwd)/swig_install --enable-silent-rules
make
make install > /dev/null # quiet installation
export PATH=$(pwd):$PATH
popd # swig
popd # external

# Build acados
mkdir -p build
pushd build
cmake -D SWIG_MATLAB=1 -D SWIG_PYTHON=1 ..
make install
popd # build
