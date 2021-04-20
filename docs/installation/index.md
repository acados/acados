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
Install `acados` as follows:
```
mkdir -p build
cd build
cmake -DACADOS_WITH_QPOASES=ON ..
# add more optional arguments e.g. -DACADOS_WITH_OSQP=OFF/ON -DACADOS_INSTALL_DIR=<path_to_acados_installation_folder> above
make install
```
NOTE: you can set the `BLASFEO_TARGET` in `<acados_root_folder>/CMakeLists.txt`.
For a list of supported targets, we refer to https://github.com/giaf/blasfeo/blob/master/README.md .
The default is `X64_AUTOMATIC`, which attempts to determine the best available target for your machine.

#### **Make**
Set the `BLASFEO_TARGET` in `<acados_root_folder>/Makefile.rule`.
Since some of the `C` examples use `qpOASES`, also set `ACADOS_WITH_QPOASES = 1` in  `<acados_root_folder>/Makefile.rule`.
For a list of supported targets, we refer to https://github.com/giaf/blasfeo/blob/master/README.md .
Install `acados` as follows:
```
make shared_library
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:<path_to_acados_folder>/lib
make examples_c
make run_examples_c
```

### Interfaces installation
For the installation of Python/MATLAB/Octave interfaces, please refer to the [Interfaces](../interfaces/index.md) page.


## Windows (for use with Matlab)

### Prerequisites
You should have the following software installed on your machine.
- Recent Matlab version, with 
- CMake for Windows
- [Windows Git Client](https://git-scm.com/download/win)

### Prepare acados build
- Locate the `cmake.exe` file. The default location is `C:\Program Files\CMake311\bin`.
- Add this path to your environment variable PATH, using the Windows GUI. To open the GUI press Windows key and type "env".
- Install mingw from MATLAB add-ons manager.
- Locate this mingw installation. The default location is `C:\ProgramData\MATLAB\SupportPackages\R2018a\3P.instrset\mingw_w64.instrset`.
- Add the subfolders `bin` and `x86_64-w64-mingw32\bin` of the above mentioned mingw installation to your environment variable PATH.

### Clone acados
Clone `acados` and its submodules by running the following from your Git shell:
```
git clone https://github.com/acados/acados.git
cd acados
git submodule update --recursive --init
```

### Build acados
Run the following from your terminal in the `<acados_root_folder>`:
```
mkdir -p build
cd build
```

Configure the `cmake` command if you want to use other external QP solvers or change the `HPIPM` and `BLASFEO` targets.
```
cmake.exe -G "MinGW Makefiles" -DACADOS_INSTALL_DIR=.. -DBUILD_SHARED_LIBS=OFF -DACADOS_WITH_QPOASES=ON -DACADOS_WITH_OSQP=ON ..
# useful options to add above:
# -DACADOS_WITH_QPOASES=ON/OFF -DACADOS_WITH_OSQP=ON/OFF -DACADOS_WITH_QPDUNES=ON/OFF ..
# -DBLASFEO_TARGET=GENERIC -DHPIPM_TARGET=GENERIC
```

In a powershell, navigate to the folder `<acados_root_folder>/build` and execute
```
mingw32-make.exe -j4
mingw32-make.exe install
```

### Try a Matlab example
- Open Matlab
- select MinGW compiler using `mbuild`/ `mex`
- go to `<acados_root_folder>/examples/acados_matlab_octave`
- run `acados_env_variables_windows`
- go to the `getting_started` subfolder
- run `minimal_example_ocp`