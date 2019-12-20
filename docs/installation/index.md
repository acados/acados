# Installation

## Linux/Mac

### Prerequisites
We assume you have: git, make, cmake installed on your system.

### Clone acados
Clone acados and its submodules by running:
```
git clone https://github.com/acados/acados.git
git submodule update --recursive --init
```

### Build and install `acados`
Both a CMake and a Makefile based build system is supported at the moment.
Please choose one and proceed with the corresponding paragraph.

#### **CMake**
Set the `BLASFEO_TARGET` in `<acados_root_folder>/CMakeLists.txt`.
For a list of supported targets, we refer to https://github.com/giaf/blasfeo/blob/master/README.md .
Install acados as follows
```
mkdir -p build
cd build
cmake .. (with optional arguments e.g. -DACADOS_WITH_OSQP=OFF/ON -DACADOS_INSTALL_DIR=<path_to_acados_installation_folder>)
make install
```

#### **Make**
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
For the installation of Python/MATLAB/Octave interfaces, please refer to the [Interfaces](../interfaces/index.md) page.


### Download CasADi:
To create external function for your problem, we suggest to use CasADi from the folder `<acados_root_folder>/external`.
Depending on the environment you want to use to generate CasADi functions from, proceed with the corresponding paragraph (Python, MATLAB, Octave):

#### **Python**

```
cd external
wget -q -nc --show-progress https://github.com/casadi/casadi/releases/download/3.4.0/casadi-linux-py35-v3.4.0-64bit.tar.gz
mkdir -p casadi-py35-v3.4.0-64bit
tar -xf casadi-linux-py35-v3.4.0-64bit.tar.gz -C casadi-py35-v3.4.0-64bit
cd ..
```

#### **Matlab**
Put CasADi binaries into `<acados_root_folder>/external/casadi-matlab` :
```
cd external
wget -q -nc --show-progress https://github.com/casadi/casadi/releases/download/3.4.0/casadi-linux-matlabR2014b-v3.4.0.tar.gz
mkdir -p casadi-matlab
tar -xf casadi-linux-matlabR2014b-v3.4.0.tar.gz -C casadi-matlab
cd ..
```

#### **Octave version 4.4 or later**
Put CasADi binaries into `<acados_root_folder>/external/casadi-octave` :
```
cd external
wget -q -nc --show-progress https://github.com/casadi/casadi/releases/download/3.4.5/casadi-linux-octave-4.4.1-v3.4.5.tar.gz
mkdir -p casadi-octave
tar -xf casadi-linux-octave-4.4.1-v3.4.5.tar.gz -C casadi-octave
```

#### **Octave version 4.2 or earliear**
Put CasADi binaries into `<acados_root_folder>/external/casadi-octave` :
```
cd external
wget -q -nc --show-progress https://github.com/casadi/casadi/releases/download/3.4.0/casadi-linux-octave-v3.4.0.tar.gz
mkdir -p casadi-octave
tar -xf casadi-linux-octave-v3.4.0.tar.gz -C casadi-octave
cd ..
```



## Windows (for use with Matlab)

### Prerequisites
You should have the following software installed on your machine.
- Recent Matlab version, with 
- CMake for Windows
- [Windows Git Client](https://git-scm.com/download/win)

### Prepare acados build
- Locate the `cmake.exe` file. The default location is `C:\Program Files\CMake311\bin`.
- Add this path to your environment variable PATH, using the GUI.

- Install mingw from MATLAB add-ons manager.
- Locate this mingw installation. The default location is `C:\ProgramData\MATLAB\SupportPackages\R2018a\3P.instrset\mingw_w64.instrset`.
- Add the subfolders `bin` and `x86_64-w64-mingw32\bin` of the above mentioned mingw installation to your environment variable PATH.

### Clone acados
Clone acados and its submodules by running the following from your Git shell:
```
git clone https://github.com/acados/acados.git
git submodule update --recursive --init
```

### Build acados
Run the following from your terminal in the `<acados_root_folder>`:
```
mkdir -p build
cd build
cmake.exe -G "MinGW Makefiles" -D BLASFEO_TARGET=GENERIC -D HPIPM_TARGET=GENERIC -D ACADOS_INSTALL_DIR=.. -DBUILD_SHARED_LIBS=OFF -DACADOS_EXAMPLES=OFF -DACADOS_UNIT_TESTS=OFF ..
mingw32-make.exe -j4
mingw32-make.exe install
```

### Try a Matlab example
- Open Matlab
- selet MinGW compiler using `mbuild`/ `mex`
- go to `<acados_root_folder>/examples/acados_matlab_octave`
- run `acados_env_variables_windows`
- go to the `getting_started` subfolder
- run `minimal_example_ocp`