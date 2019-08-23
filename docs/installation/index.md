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

1. Download CasADi:
To create external function for your problem, we suggest to use CasADi and use it from `<acados_root_folder>/external`.
Depending on the environment you want to use to generate CasADi functions from, proceed with the corresponding paragraph (Python, MATLAB, Octave):

    ### **Python**

    ```
    cd external
    wget -q -nc --show-progress https://github.com/casadi/casadi/releases/download/3.4.0/casadi-linux-py35-v3.4.0-64bit.tar.gz
    mkdir -p casadi-py35-v3.4.0-64bit
    tar -xf casadi-linux-py35-v3.4.0-64bit.tar.gz -C casadi-py35-v3.4.0-64bit
    cd ..
    ```

    ### **Matlab**
    Put CasADi binaries into `<acados_root_folder>/external/casadi-matlab` :
    ```
    cd external
    wget -q -nc --show-progress https://github.com/casadi/casadi/releases/download/3.4.0/casadi-linux-matlabR2014b-v3.4.0.tar.gz
    mkdir -p casadi-matlab
    tar -xf casadi-linux-matlabR2014b-v3.4.0.tar.gz -C casadi-matlab
    cd ..
    ```

    ### **Octave version 4.4 or later**
    Put CasADi binaries into `<acados_root_folder>/external/casadi-octave` :
    ```
    cd external
    wget -q -nc --show-progress https://github.com/casadi/casadi/releases/download/3.4.5/casadi-linux-octave-4.4.1-v3.4.5.tar.gz
    mkdir -p casadi-octave
    tar -xf casadi-linux-octave-4.4.1-v3.4.5.tar.gz -C casadi-octave
    ```

    ### **Octave version 4.2 or earliear**
    Put CasADi binaries into `<acados_root_folder>/external/casadi-octave` :

    ```
    cd external
    wget -q -nc --show-progress https://github.com/casadi/casadi/releases/download/3.4.0/casadi-linux-octave-v3.4.0.tar.gz
    mkdir -p casadi-octave
    tar -xf casadi-linux-octave-v3.4.0.tar.gz -C casadi-octave
    cd ..
    ```

1. Build and install `acados`.
Both a CMake and a Makefile based build system is supported at the moment.
Please choose one and proceed with the corresponding paragraph.

    ### **CMake**
    Set the `BLASFEO_TARGET` in `<acados_root_folder>/CMakeLists.txt`.
    For a list of supported targets, we refer to https://github.com/giaf/blasfeo/blob/master/README.md .
    Install acados as follows
    ```
    mkdir -p build
    cd build
    cmake .. (with optional arguments e.g. -DACADOS_WITH_OSQP=OFF/ON -DACADOS_INSTALL_DIR=<path_to_acados_installation_folder>)
    make install
    ```

    ### **Make**
    Set the `BLASFEO_TARGET` in `<acados_root_folder>/Makefile.rule`.
    For a list of supported targets, we refer to https://github.com/giaf/blasfeo/blob/master/README.md .
    Install acados as follows
    ```
    make shared_library
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:<path_to_acados_folder>/lib
    make examples_c
    make run_examples_c
    ```

### Interfaces installation
For the installation of Python/MATLAB interfaces, please refer to the [Interfaces](../interfaces/index.md) page.
