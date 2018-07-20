#!/usr/bin/env bash -xe

sudo add-apt-repository -y ppa:octave/stable
sudo apt-get update -yqq
sudo apt-get --allow-unauthenticated install -yqq $CXX $CC $COVERAGE libgsl0-dev liblapack-dev libopenblas-dev liboctave-dev mingw-w64 bsdtar

pip install numpy scipy matplotlib

# Windows libs for openblas
pushd $HOME
wget -q https://sourceforge.net/projects/openblas/files/v0.2.19/OpenBLAS-v0.2.19-Win64-int32.zip
mkdir -p WindowsLibs
bsdtar xvf OpenBLAS-v0.2.19-Win64-int32.zip --strip-components=1 -C WindowsLibs
popd

pushd external
wget http://bitbucket.org/eigen/eigen/get/3.2.10.tar.gz
mkdir -p eigen
tar -xf 3.2.10.tar.gz --strip-components=1 -C eigen

wget https://github.com/casadi/casadi/releases/download/3.4.0/casadi-linux-octave-v3.4.0.tar.gz
mkdir -p casadi-octave-v3.4.0
tar -xf casadi-linux-octave-v3.4.0.tar.gz -C casadi-octave-v3.4.0

wget https://github.com/casadi/casadi/releases/download/3.4.0/casadi-linux-py35-v3.4.0-64bit.tar.gz
mkdir -p casadi-linux-py35-v3.4.0-64bit
tar -xf casadi-linux-py35-v3.4.0-64bit.tar.gz -C casadi-linux-py35-v3.4.0-64bit
export PYTHONPATH=$(pwd)/casadi-linux-py35-v3.4.0-64bit:$PYTHONPATH
export CASADIPATH=$(pwd)/casadi-linux-py35-v3.4.0-64bit

wget https://github.com/casadi/casadi/releases/download/3.4.0/casadi-linux-matlabR2014b-v3.4.0.tar.gz
mkdir -p casadi-linux-matlabR2014b-v3.4.0
tar -xf casadi-linux-matlabR2014b-v3.4.0.tar.gz -C casadi-linux-matlabR2014b-v3.4.0
export MATLABPATH=$(pwd)/casadi-linux-matlabR2014b-v3.4.0:$MATLABPATH

# wget -q http://icl.cs.utk.edu/lapack-for-windows/libraries/VisualStudio/3.7.0/Dynamic-MINGW/Win64/libblas.lib
# sudo mv libblas.lib /usr/x86_64-w64-mingw32/lib/libblas.a
# wget -q http://icl.cs.utk.edu/lapack-for-windows/libraries/VisualStudio/3.7.0/Dynamic-MINGW/Win64/liblapack.lib
# sudo mv liblapack.lib /usr/x86_64-w64-mingw32/lib/liblapack.a

pushd swig

./autogen.sh
./configure --prefix=$(pwd)/swig_install --enable-silent-rules
make
make install > /dev/null # quiet installation
export PATH=$(pwd):$PATH

popd # swig

popd # external
