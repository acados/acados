#!/usr/bin/env bash -xe

sudo add-apt-repository -y ppa:octave/stable
sudo apt-get update -yqq
sudo apt-get install -yqq $CXX $CC $COVERAGE libgsl0-dev liblapack-dev libopenblas-dev liboctave-dev mingw-w64 bsdtar

pip install numpy scipy matplotlib

# Windows libs for openblas
pushd $HOME
wget -q https://sourceforge.net/projects/openblas/files/v0.2.19/OpenBLAS-v0.2.19-Win64-int32.zip
mkdir -p WindowsLibs
bsdtar xvf OpenBLAS-v0.2.19-Win64-int32.zip --strip-components=1 -C WindowsLibs
popd

pushd external
wget -q http://bitbucket.org/eigen/eigen/get/3.2.10.tar.gz
mkdir -p eigen
tar -xf 3.2.10.tar.gz --strip-components=1 -C eigen

wget -q https://sourceforge.net/projects/casadi/files/CasADi/3.2.2/linux/casadi-octave-v3.2.2.tar.gz
mkdir -p casadi-octave-v3.2.2
tar -xf casadi-octave-v3.2.2.tar.gz -C casadi-octave-v3.2.2

wget -q https://sourceforge.net/projects/casadi/files/CasADi/3.1.1/linux/casadi-py35-np1.9.1-v3.1.1.tar.gz
mkdir -p casadi-py35-np1.9.1-v3.1.1
tar -xf casadi-py35-np1.9.1-v3.1.1.tar.gz -C casadi-py35-np1.9.1-v3.1.1
export PYTHONPATH=$(pwd)/casadi-py35-np1.9.1-v3.1.1:$PYTHONPATH

wget -q https://sourceforge.net/projects/casadi/files/CasADi/3.2.3/linux/casadi-matlabR2014b-v3.2.3.tar.gz
mkdir -p casadi-matlabR2014b-v3.2.3
tar -xf casadi-matlabR2014b-v3.2.3.tar.gz -C casadi-matlabR2014b-v3.2.3
export MATLABPATH=$(pwd)/casadi-matlabR2014b-v3.2.3:$MATLABPATH

wget -q http://icl.cs.utk.edu/lapack-for-windows/libraries/VisualStudio/3.7.0/Dynamic-MINGW/Win64/libblas.lib
sudo mv libblas.lib /usr/x86_64-w64-mingw32/lib/libblas.a
wget -q http://icl.cs.utk.edu/lapack-for-windows/libraries/VisualStudio/3.7.0/Dynamic-MINGW/Win64/liblapack.lib
sudo mv liblapack.lib /usr/x86_64-w64-mingw32/lib/liblapack.a

pushd swig

./autogen.sh
./configure --prefix=$(pwd)/swig_install --enable-silent-rules
make
make install > /dev/null # quiet installation
export PATH=$(pwd):$PATH

popd # swig

popd # external
