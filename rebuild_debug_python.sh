#!/bin/bash -xe

# Install swig
pushd external

export CASADIPATH=$(pwd)/casadi-py35-v3.4.0-64bit
export PYTHONPATH=$CASADIPATH:$PYTHONPATH
# will not work with custom install dir
export PYTHONPATH=/usr/local/lib:$PYTHONPATH
export MATLABPATH=$(pwd)/casadi-matlabR2014b-v3.4.0:$MATLABPATH

pushd swig
./autogen.sh
./configure --prefix=$(pwd)/swig_install --enable-silent-rules
make
make install > /dev/null # quiet installation
export PATH=$(pwd):$PATH
popd # swig
popd # external

# Build acados
sudo rm -rf build
mkdir -p build
pushd build
cmake -D SWIG_PYTHON=1 -DCMAKE_BUILD_TYPE=DEBUG ..
sudo make install -j4 -l4
popd # build
