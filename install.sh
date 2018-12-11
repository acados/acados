#!/bin/bash -xe
INSTALLPATH=/usr/local

# Install dependencies
# sudo add-apt-repository ppa:george-edison55/cmake-3.x
# sudo apt-add-repository ppa:octave/stable
sudo apt-get update
sudo apt-get install octave liboctave-dev
sudo apt-get install libgsl0-dev liblapack-dev libopenblas-dev libeigen3-dev python3-tk automake libpcre3-dev git cmake python3-dev
sudo apt-get install byacc # swig
sudo apt-get install python3-scipy python3-numpy python3-matplotlib

# Get CasADi for octave, python and matlab
pushd external
wget -q -nc https://github.com/casadi/casadi/releases/download/3.4.0/casadi-linux-octave-v3.4.0.tar.gz
mkdir -p casadi-octave-v3.4.0
tar -xf casadi-linux-octave-v3.4.0.tar.gz -C casadi-octave-v3.4.0

wget -q -nc https://github.com/casadi/casadi/releases/download/3.4.0/casadi-linux-py35-v3.4.0-64bit.tar.gz
mkdir -p casadi-py35-v3.4.0-64bit
tar -xf casadi-linux-py35-v3.4.0-64bit.tar.gz -C casadi-py35-v3.4.0-64bit
export CASADIPATH=$(pwd)/casadi-py35-v3.4.0-64bit
export PYTHONPATH=$CASADIPATH:$PYTHONPATH
# will not work with custom install dir
export PYTHONPATH=$INSTALLPATH/lib:$PYTHONPATH

wget -q -nc https://github.com/casadi/casadi/releases/download/3.4.0/casadi-linux-matlabR2014b-v3.4.0.tar.gz
mkdir -p casadi-matlabR2014b-v3.4.0
tar -xf casadi-linux-matlabR2014b-v3.4.0.tar.gz -C casadi-matlabR2014b-v3.4.0
export MATLABPATH=$(pwd)/casadi-matlabR2014b-v3.4.0:$MATLABPATH


# Get all git submodules
git submodule update --recursive --init


# Install swig
pushd swig
./autogen.sh
./configure --prefix=$(pwd)/swig_install --enable-silent-rules
make -j4 -l4
make install > /dev/null # quiet installation
export PATH=$(pwd):$PATH
popd # swig
popd # external

# Build acados
mkdir -p build
pushd build
cmake -D SWIG_MATLAB=1 -D SWIG_PYTHON=1 -D ACADOS_INSTALL_DIR=$INSTALLPATH ..
make -j4 -l4
make install
popd # build
