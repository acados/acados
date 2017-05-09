#!/usr/bin/env bash -xe

sudo add-apt-repository -y ppa:octave/stable
sudo add-apt-repository -y ppa:ubuntu-toolchain-r/test
if [[ "$CC" == "clang-3.7" ]];
then
    echo 'deb http://apt.llvm.org/precise/ llvm-toolchain-precise-3.7 main' > /tmp/myppa.list
    sudo cp /tmp/myppa.list /etc/apt/sources.list.d/
    rm /tmp/myppa.list
fi
sudo apt-get update --force-yes -yq
sudo apt-get install $CXX $CC $GFORTRAN libgsl0-dev liblapack-dev libopenblas-dev liboctave-dev

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

./autogen.sh
./configure --prefix=$(pwd)/swig_install --enable-silent-rules
make
make install
export PATH=$(pwd):$PATH

popd # swig

popd # external
