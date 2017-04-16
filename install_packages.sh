#!/bin/sh

pip install numpy
sudo add-apt-repository -y ppa:octave/stable
sudo apt-get update -yq
sudo apt-get install octave liboctave-dev valgrind

pushd external
wget http://bitbucket.org/eigen/eigen/get/3.2.10.tar.gz
mkdir eigen
tar -xvf 3.2.10.tar.gz --strip-components=1 -C eigen

wget http://files.casadi.org/3.1.1/linux/casadi-octave-v3.1.1.tar.gz
mkdir casadi-octave-v3.1.1
tar -xvf casadi-octave-v3.1.1.tar.gz -C casadi-octave-v3.1.1

wget http://files.casadi.org/3.1.1/linux/casadi-py35-np1.9.1-v3.1.1.tar.gz
mkdir casadi-py35-np1.9.1-v3.1.1
tar -xvf casadi-py35-np1.9.1-v3.1.1.tar.gz -C casadi-py35-np1.9.1-v3.1.1

export PYTHONPATH=$(pwd)/casadi-py35-np1.9.1-v3.1.1:$PYTHONPATH
popd
