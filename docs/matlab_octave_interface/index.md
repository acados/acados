# Matlab + Simulink and Octave Interface

In order to use `acados` from Octave or Matlab, you need to create the `acados` shared libraries  using either the `CMake` or `Make` build system, as described [on the installation page](../installation/index.md).

## Getting started
Check out the examples in [`<acados_root>/examples/acados_matlab_octave/getting_started`](https://github.com/acados/acados/tree/master/examples/acados_matlab_octave/getting_started) to get started with this interface.

The examples require an installation of `CasADi` to generate the model functions.
The `getting_started` example offers the option to attempt to automatically download the correct version in the recommended folder.
Detailed instructions for a manual installation can be found in the last section of this page [Setup CasADi](#setup-casadi).

The problem formulation is stated in [this PDF](https://github.com/acados/acados/tree/master/docs/problem_formulation/problem_formulation_ocp_mex.pdf).

### Linux / macOS
To run the examples, navigate into the folder of the example you want to run and execute the following command:
```
source env.sh # Which can be found in the folder of one of the examples
```

If you want to run an `acados` example from another folder, you need to export the environment variable `ACADOS_INSTALL_DIR` properly.
In the `env.sh` file it is assumed that `ACADOS_INSTALL_DIR` is two folders above the directory, in which the example is located.

Afterwards, launch `Matlab` or `Octave` from the same shell.

If you want to run the examples in a different folder, please close the current shell and open a new one to repeat the procedure: this ensures the correct setting of the environment variables.

### Windows
1. Open `Matlab` and navigate into [`<acados_root>/examples/acados_matlab_octave`](https://github.com/acados/acados/blob/master/examples/acados_matlab_octave).
2. Run [`acados_env_variables_windows`](https://github.com/acados/acados/blob/master/examples/acados_matlab_octave/acados_env_variables_windows.m) to export the environment variable `ACADOS_INSTALL_DIR`.
3. Navigate into [`<acados_root>/examples/acados_matlab_octave/getting_started`](https://github.com/acados/acados/tree/master/examples/acados_matlab_octave/getting_started) and run one of the examples.


## Interface structure
The interface consists of two parts:
1. the [native MEX interface for rapid prototyping](#native-mex-rapid-prototyping)
2. the [template based interface for code generation, use in Simulink and embedded deployment](#templates)

The structure is visualized below:
```eval_rst
.. image:: ./mex_interface_overview.png
.. This is a comment.
```

### Native MEX (rapid prototyping)
This interface makes a broad set of `acados` functionalities available from `Matlab` and `Octave`.

### Templates
There is the option to generate embeddable `C` code from Matlab.
The workflow uses the same templates as the Python interface (see [`Python interface`](../python_interface/index.md)) and the `Tera` renderer.
After creating an acados solver `ocp`, you can use the routine `ocp.generate_c_code()` to generate `C` code which can be used for embedded applications.
These templates can be found in [`<acados_root>/interfaces/acados_template/acados_template/c_templates_tera`](https://github.com/acados/acados/tree/master/interfaces/acados_template/acados_template/c_templates_tera).

Note: This part of the Matlab/Octave interface does not yet support all features of the one mentioned before.

## Options documentation
A table and explanation of various options of the native `MEX` interface can be found in [this spreadsheet](https://docs.google.com/spreadsheets/d/1rVRycLnCyaWJLwnV47u30Vokp7vRu68og3OhlDbSjDU/edit?usp=sharing) (thanks to [@EnricaSo](https://github.com/EnricaSo)).

For the template based part of the `Matlab` interface, we refer to [the docstring based documentation of the Python interface](../python_interface/index.md).

## Simulink
The templates mentioned [above](#templates) also contain templated S-functions and corresponding make functions for Matlab for both the OCP solver and the acados integrator.

A basic Simulink example can be found in [`<acados_root>/examples/acados_python/getting_started/simulink_example.m`](https://github.com/acados/acados/blob/master/examples/acados_matlab_octave/getting_started/simulink_example.m)

A more advanced Simulink example which showcases how to customize the inputs and outputs of the Simulink block corrsponding to the solver can be found in [`<acados_root>/examples/acados_python/getting_started/simulink_example.m`](https://github.com/acados/acados/blob/master/examples/acados_matlab_octave/getting_started/simulink_example_advanced.m)

If you want a more advanced interaction with the `acados` solver via Simulink, feel free to edit the corresponding templates in [`<acados_root>/interfaces/acados_template/acados_template/c_templates_tera`](https://github.com/acados/acados/tree/master/interfaces/acados_template/acados_template/c_templates_tera) to add more inputs or outputs.


## Setup CasADi
To create external function for your problem, we suggest to use `CasADi` from the folder `<acados_root_folder>/external`.
Depending on the environment you want to use to generate `CasADi` functions from, proceed with the corresponding paragraph (Matlab, Octave).

### **Matlab**
Download and extract the `CasADi` binaries into `<acados_root_folder>/external/casadi-matlab`:
```
cd external
wget -q -nc --show-progress https://github.com/casadi/casadi/releases/download/3.4.0/casadi-linux-matlabR2014b-v3.4.0.tar.gz
mkdir -p casadi-matlab
tar -xf casadi-linux-matlabR2014b-v3.4.0.tar.gz -C casadi-matlab
cd ..
```

### **Octave version 4.4 or later**
Download and extract the `CasADi` binaries into `<acados_root_folder>/external/casadi-octave`:
```
cd external
wget -q -nc --show-progress https://github.com/casadi/casadi/releases/download/3.4.5/casadi-linux-octave-4.4.1-v3.4.5.tar.gz
mkdir -p casadi-octave
tar -xf casadi-linux-octave-4.4.1-v3.4.5.tar.gz -C casadi-octave
```

### **Octave version 4.2 or earliear**
Download and extract the `CasADi` binaries into `<acados_root_folder>/external/casadi-octave`:
```
cd external
wget -q -nc --show-progress https://github.com/casadi/casadi/releases/download/3.4.0/casadi-linux-octave-v3.4.0.tar.gz
mkdir -p casadi-octave
tar -xf casadi-linux-octave-v3.4.0.tar.gz -C casadi-octave
cd ..
```