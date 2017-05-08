#!/usr/bin/env bash -xe

sudo apt-get install libgsl0-dev liblapack-dev libopenblas-dev cppcheck

sudo add-apt-repository -y ppa:octave/stable
sudo apt-get update -yq
sudo apt-get install octave liboctave-dev valgrind

pip install numpy scipy matplotlib

pushd external
wget http://bitbucket.org/eigen/eigen/get/3.2.10.tar.gz
mkdir -p eigen
tar -xvf 3.2.10.tar.gz --strip-components=1 -C eigen

wget http://files.casadi.org/3.1.1/linux/casadi-octave-v3.1.1.tar.gz
mkdir -p casadi-octave-v3.1.1
tar -xvf casadi-octave-v3.1.1.tar.gz -C casadi-octave-v3.1.1

wget http://files.casadi.org/3.1.1/linux/casadi-py35-np1.9.1-v3.1.1.tar.gz
mkdir -p casadi-py35-np1.9.1-v3.1.1
tar -xvf casadi-py35-np1.9.1-v3.1.1.tar.gz -C casadi-py35-np1.9.1-v3.1.1
export PYTHONPATH=$(pwd)/casadi-py35-np1.9.1-v3.1.1:$PYTHONPATH

wget https://sourceforge.net/projects/casadi/files/CasADi/3.1.1/linux/casadi-matlabR2014b-v3.1.1.tar.gz
mkdir -p casadi-matlabR2014b-v3.1.1
tar -xvf casadi-matlabR2014b-v3.1.1.tar.gz -C casadi-matlabR2014b-v3.1.1
export MATLABPATH=$(pwd)/casadi-matlabR2014b-v3.1.1:$MATLABPATH

pushd swig

./configure --prefix=

popd # swig

popd # external
