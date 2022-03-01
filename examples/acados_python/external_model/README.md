# How to run this example

## Clone and build CasADi locally
```
# somewhere on your machine
git clone git@github.com:casadi/casadi.git
mkdir build
cd build
cmake ..
sudo make install -j4
```

## Build external library
```
cd <to folder of this README>/test_external_lib
mkdir build
cd build
cmake ..
make -j4
```

## Run example
```
python minimal_example_external_ode.py
```