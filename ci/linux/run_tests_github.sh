# source virtualenv
source ../acados/acadosenv/bin/activate;
which python;

# run tests
ctest -C $BUILD_TYPE -V;
