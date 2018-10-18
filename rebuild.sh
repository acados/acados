# Get all git submodules
# git submodule update --recursive --init


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
rm -rf build
mkdir -p build
pushd build
cmake -D SWIG_MATLAB=1 -D SWIG_PYTHON=1 ..
sudo make install
popd # build
