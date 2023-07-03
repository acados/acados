# Installation

## Linux/Mac

### Prerequisites
We assume you have: git, make, cmake installed on your system.

### Clone acados
Clone acados and its submodules by running:
```
git clone https://github.com/acados/acados.git
cd acados
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
make install -j4
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
NOTE: On MacOS `DYLD_LIBRARY_PATH` should be used instead of `LD_LIBRARY_PATH`.

### Interfaces installation
For the installation of Python/MATLAB/Octave interfaces, please refer to the [Interfaces](../interfaces/index.md) page.

## Windows 10+ (WSL)

### Prerequisites

- Install [Ubuntu on WSL](https://ubuntu.com/wsl) using the [Microsoft Store](https://apps.microsoft.com/store/detail/ubuntu/9PDXGNCFSCZV).

- Start Ubuntu, the shell should pop up.

- Update apt

```bash
apt-get update
```

- Install cmake, build-essentials, pip, virtualenv

```bash
apt-get install cmake build-essential python3-pip python3-virtualenv
```

### Clone, Build and Install acados

- Navigate to the directory where you would like to install acados. For example

```
cd /mnt/c/Users/Documents/
```

- Follow the [Linux/Mac workflow above and install acados using cmake](#linux-mac).

### Interfaces installation

For the installation of Python/MATLAB/Octave interfaces, please refer to the [Interfaces](../interfaces/index.md) page.

## Windows (for use with Matlab)

Disclaimer: The high-level interfaces on Windows are not tested on Github Actions.

### Prerequisites
You should have the following software installed on your machine.
- Recent Matlab version
- CMake for Windows
- [Windows Git Client](https://git-scm.com/download/win)

### Clone acados
Clone `acados` and its submodules by running the following from your Git shell:
```
git clone https://github.com/acados/acados.git
cd acados
git submodule update --recursive --init
```

### Prepare acados build (minGW)
- Install mingw from MATLAB add-ons manager.
- Add the following paths to your environment variable `PATH`, using the Windows GUI. To open the GUI press Windows key and type "env".
1. The path to the `cmake.exe` file. The default location is `C:\Program Files\CMake\bin`.
2. The path to the MATLAB installation of `mingw32-make.exe` file. The default location is `C:\ProgramData\MATLAB\SupportPackages\R2018a\3P.instrset\mingw_w64.instrset\bin`.
- You can check whether the modification to `PATH` variable is in effect by executing in cmd `echo %PATH%` or in PowerShell `echo $env:path`.

### Automatic build of acados (minGW)
Run the following in Matlab in the folder `<acados_root_folder>/interfaces/acados_matlab_octave`:
```
acados_install_windows()
```

This will build acados with the standard options and install the external dependencies; casadi and terra renderer

The `acados_install_windows(CMakeConfigString)` script can take an optional argument:
- CMake configuration string: Configuration options for CMake. The default is `-DBUILD_SHARED_LIBS=OFF -DACADOS_WITH_OSQP=OFF`

### Build acados manually (minGW)
If the automated install procedure does not work acados can be built manually using these steps.

Run the following from a powershell in the `<acados_root_folder>`:
```
$ACADOS_INSTALL_DIR=$(pwd)
mkdir -p build
cd build
```

Configure the `cmake` command if you want to use other external QP solvers or change the `HPIPM` and `BLASFEO` targets.
```
cmake.exe -G "MinGW Makefiles" -DACADOS_INSTALL_DIR="$ACADOS_INSTALL_DIR" -DBUILD_SHARED_LIBS=OFF -DACADOS_WITH_OSQP=ON ..
# useful options to add above:
# -DACADOS_WITH_QPOASES=ON/OFF -DACADOS_WITH_OSQP=ON/OFF -DACADOS_WITH_QPDUNES=ON/OFF ..
# -DBLASFEO_TARGET=GENERIC -DHPIPM_TARGET=GENERIC
# NOTE: check the output of cmake: -- Installation directory: should be <acados_root_folder>,
#     if this is not the case, set -DACADOS_INSTALL_DIR=<acados_root_folder> explicitly above.
```

In a powershell, navigate to the folder `<acados_root_folder>/build` and execute
```
mingw32-make.exe -j4
mingw32-make.exe install
```

### Try a Matlab example
- Open Matlab
- select MinGW compiler using `mex`
- go to `<acados_root_folder>/examples/acados_matlab_octave`
- run `acados_env_variables_windows`
- go to the `getting_started` subfolder
- run `minimal_example_ocp`

### Workflow with Microsoft Visual C Compiler (MSVC)
Note: this workflow is preliminary and not thoroughly tested.
(Tested once with MSVC 2017 and MSVC 2019 on May 2021)

- clone acados (see above)
- use the `Developer Command Prompt for VS`, navigate to `<acados_root_folder>/build` and run
```
cmake -G "Visual Studio 15 2017 Win64" -DBLASFEO_TARGET=GENERIC -DACADOS_INSTALL_DIR=.. -DBUILD_SHARED_LIBS=OFF ..
# respectively for MVSC 2019
# cmake -G "Visual Studio 16 2019" -DBLASFEO_TARGET=GENERIC -DACADOS_INSTALL_DIR=.. -DBUILD_SHARED_LIBS=OFF ..
cmake --build . -j10 --target INSTALL --config Release
```
- In Matlab, run `mex -setup C` and select the same `MSVC` version.
- Try a Matlab example (see above).
