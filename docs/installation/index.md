# Installation

## Optional requirements
Some functionalities in acados require CasADi (version 3.4.0) to be installed on your system.
To install CasADi, you can follow the installation instructions [here](https://github.com/casadi/casadi/wiki/InstallationInstructions)

## Installation
Both a CMake and a Makefile based build system are supported at the moment.

1. Clone acados
    ```
    git clone https://github.com/acados/acados.git
    ```

1. Initialize all submodules
    ```
    git submodule update --recursive --init
    ```

1. Download CasADi into the `<acados_root_folder>/external` folder:
    ```
    cd external
    ```
    and, depending on your preferred CasADi interface (Python, MATLAB, Octave):

    ```
    wget -q -nc http://files.casadi.org/download/3.4.0/casadi-linux-py35-v3.4.0-64bit.tar.gz
    mkdir -p casadi-py35-v3.4.0-64bit
    tar -xf casadi-linux-py35-v3.4.0-64bit.tar.gz -C casadi-py35-v3.4.0-64bit
    ```

    ```
    wget -q -nc http://files.casadi.org/download/3.4.0/casadi-linux-matlabR2014b-v3.4.0.tar.gz
    mkdir -p casadi-matlabR2014b-v3.4.0
    tar -xf casadi-linux-matlabR2014b-v3.4.0.tar.gz -C casadi-matlabR2014b-v3.4.0
    cd ..
    ```

    ```
    wget -q -nc http://files.casadi.org/download/3.4.0/casadi-linux-octave-v3.4.0.tar.gz
    mkdir -p casadi-octave-v3.4.0
    tar -xf casadi-linux-octave-v3.4.0.tar.gz -C casadi-octave-v3.4.0
    ```

1. To build and install `acados` library you can either use `Makefile`
   or `CMake`


### CMake

```bash
mkdir -p build
cd build
cmake ..
make install
```

Optional CMake arguments:
* `-DACADOS_WITH_OSQP=OFF/ON`
* `-DACADOS_INSTALL_DIR=<path_to_acados_installation_folder>`


### Makefile

```bash
make acados_shared
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:<path_to_acados_folder>/lib
make examples_c
make run_examples_c
```
### Interfaces installation
For the installation of Python/MATLAB interfaces, please refer to the [Interfaces](../interfaces/index.md) page.
