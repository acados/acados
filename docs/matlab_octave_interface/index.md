# MATLAB + Simulink and Octave Interface

In order to use `acados` from Octave or MATLAB, you need to create the `acados` shared libraries using either the `CMake` or `Make` build system, as described [on the installation page](../installation/index.md).

## Getting started
Check out the examples [`minimal_example_ocp.m`](https://github.com/acados/acados/tree/main/examples/acados_matlab_octave/getting_started/minimal_example_ocp.m) and [`minimal_example_sim.m`](https://github.com/acados/acados/tree/main/examples/acados_matlab_octave/getting_started/minimal_example_sim.m) to get started with the MATLAB interface of `acados`.
Note that `acados` currently supports both an old MATLAB interface (< v0.4.0) as well as the new one (>= v0.4.0).
Unfortunately, not all MATLAB examples have been ported to the new interface yet.
If you are new to `acados` please start with [those examples](https://github.com/acados/acados/issues/1196#issuecomment-2311822122) that use the new interface already.


The examples require an installation of `CasADi` to generate the model functions.
The `getting_started` example offers the option to attempt to automatically download the correct version in the recommended folder.
Detailed instructions for a manual installation can be found in the last section of this page [Setup CasADi](#setup-casadi).

The problem formulation is stated in [this PDF](https://github.com/acados/acados/tree/main/docs/problem_formulation/problem_formulation_ocp_mex.pdf).



## Export environment variables
In order to run the examples, some environment variables need to be exported.
Instead of running the scripts below, you can modify an `rc` file, like `.bashrc` when launching MATLAB from bash,
[`.matlab7rc.sh`](https://discourse.acados.org/t/matlab-mex-more-elegant-way-to-setup-env-sh/62/4) or `startup.m` to always have those environment variables defined when starting `MATLAB`.

### Linux / macOS
Navigate into the folder of the example you want to run and execute the following command:
```
source env.sh # Which can be found in the folder of one of the examples
```

If you want to run an `acados` example from another folder, you need to export the environment variable `ACADOS_INSTALL_DIR` properly.
In the `env.sh` file it is assumed that `ACADOS_INSTALL_DIR` is two folders above the directory, in which the example is located.

Afterwards, launch `MATLAB` or `Octave` from the same shell.

If you want to run the examples in a different folder, please close the current shell and open a new one to repeat the procedure: this ensures the correct setting of the environment variables.

### Windows
1. Open `MATLAB` and navigate into [`<acados_root>/examples/acados_matlab_octave`](https://github.com/acados/acados/blob/main/examples/acados_matlab_octave).
2. Run [`acados_env_variables_windows`](https://github.com/acados/acados/blob/main/examples/acados_matlab_octave/acados_env_variables_windows.m) to export the environment variable `ACADOS_INSTALL_DIR`.
3. Navigate into [`<acados_root>/examples/acados_matlab_octave/getting_started`](https://github.com/acados/acados/tree/main/examples/acados_matlab_octave/getting_started) and run one of the examples.


## Interface structure
The interface allows one to conveniently and compactly formulate an OCP (or IVP) and specify solver options.
The nonlinear problem functions can be formulated using CasADi symbolics which are generated as C code with the required derivatives using automatic differentiation.
The whole problem description is written to a json-file which is then used to render different templates, via the `Tera` renderer.
These are the same templates as in the Python interface (see [`Python interface`](../python_interface/index.md)).
In addition to a `MEX` wrapper it contains all the `C` code that is needed for embedded deployment.
These templates can be found in [`<acados_root>/interfaces/acados_template/acados_template/c_templates_tera`](https://github.com/acados/acados/tree/main/interfaces/acados_template/acados_template/c_templates_tera).

## Options documentation
For the template based part of the `MATLAB` interface, we refer to [the docstring based documentation of the Python interface](../python_interface/index.md).

## Simulink
The templates mentioned [above](#templates) also contain templated S-functions and corresponding make functions for both the OCP solver and the acados integrator.

A basic Simulink example can be found in [`<acados_root>/examples/acados_python/getting_started/simulink_example.m`](https://github.com/acados/acados/blob/main/examples/acados_matlab_octave/getting_started/simulink_example.m)

A more advanced Simulink example which showcases how to customize the inputs and outputs of the Simulink block corresponding to the solver can be found in [`<acados_root>/examples/acados_python/getting_started/simulink_example_advanced.m`](https://github.com/acados/acados/blob/main/examples/acados_matlab_octave/getting_started/simulink_example_advanced.m)

### List of possible inputs
This is a list of possible inputs to the Simulink block of an OCP solver which can be activated by setting the corresponding values in the acados Simulink options.

| Port Name         | Description                                                                                       | Shape                 | MOCP Support |
|-------------------|-----------------------------------------------------------------------------------------------|-----------------------|--------------|
| `lbx_0`           | Lower bound on x for stage 0                                                                  | `nbx_0`               | Yes          |
| `ubx_0`           | Upper bound on x for stage 0                                                                  | `nbx_0`               | Yes          |
| `parameter_traj`  | Parameters - concatenated for all stages 0 to N                                               | `sum(np_i), i=0,..., N`            | Yes          |
| `p_global`        | Global parameters - first value indicates if update should be performed (0 = no update), followed by new numerical values of `p_global` | `1 + np_global`        | Yes          |
| `y_ref_0`         | Reference `y_ref` at stage 0                                                                  | `ny_0`                | Yes          |
| `y_ref`           | `y_ref` concatenated for stages 1 to N-1                                                      | `(N-1) * ny`          | No           |
| `y_ref_e`         | Reference `y_ref` at stage N                                                                  | `ny_e`                | Yes          |
| `lbx`             | Lower bound x values concatenated for stages 1 to N-1                                         | `sum(nbx_i), for i = 1,..., N-1`         | Yes           |
| `ubx`             | Upper bound x values concatenated for stages 1 to N-1                                         | `sum(nbx_i), for i = 1,..., N-1`         | Yes           |
| `lbx_e`           | Lower bound x at shooting node N                                                              | `nbx_e`               | Yes          |
| `ubx_e`           | Upper bound x at shooting node N                                                              | `nbx_e`               | Yes          |
| `lbu`             | Lower bound u values concatenated for stages 0 to N-1                                         | `sum(nbu_i), for i = 0,..., N-1` | Yes          |
| `ubu`             | Upper bound u values concatenated for stages 0 to N-1                                         | `sum(nbu_i), for i = 0,..., N-1` | Yes          |
| `lg`              | Lower bound g values concatenated for stages 0 to N-1                                         | `sum(ng_i), for i = 0,..., N-1`  | No           |
| `ug`              | Upper bound g values concatenated for stages 0 to N-1                                         | `sum(ng_i), for i = 0,..., N-1`  | No           |
| `lh`              | Lower bound h values concatenated for stages 1 to N-1                                         | `sum(nh_i), for i = 1,..., N-1`  | Yes           |
| `uh`              | Upper bound h values concatenated for stages 1 to N-1                                         | `sum(nh_i), for i = 1,..., N-1`  | Yes           |
| `lh_0`            | Lower bound h at stage 0                                                                      | `nh_0`                | Yes          |
| `uh_0`            | Upper bound h at stage 0                                                                      | `nh_0`                | Yes          |
| `lh_e`            | Lower bound h at stage N                                                                      | `nh_e`                | Yes          |
| `uh_e`            | Upper bound h at stage N                                                                      | `nh_e`                | Yes          |
| `cost_W_0`        | Cost matrix `W_0` in column-major format                                                      | `ny_0 * ny_0`         | No           |
| `cost_W`          | Cost matrix `W` in column-major format                                                        | `ny * ny`             | No           |
| `cost_W_e`        | Cost matrix `W_e` in column-major format                                                      | `ny_e * ny_e`         | No           |
| `cost_zl`         | Cost `zl` for all nodes 0 to N                                                                | `sum(ns_i) for i = 0,..., N`     | Yes          |
| `cost_zu`         | Cost `zu` for all nodes 0 to N                                                                | `sum(ns_i) for i = 0,..., N`     | Yes          |
| `cost_Zl`         | Cost `Zl` for all nodes 0 to N                                                                | `sum(ns_i) for i = 0,..., N`     | Yes          |
| `cost_Zu`         | Cost `Zu` for all nodes 0 to N                                                                | `sum(ns_i) for i = 0,..., N`     | Yes          |
| `reset_solver`    | Determines if the solver's iterate is set to all zeros before other initializations            | `1`                   | Yes          |
| `ignore_inits`    | Determines if initialization (`x_init`, `u_init`, `pi_init`, `slacks_init`) is set (0) or ignored (1) | `1`                   | Yes          |
| `x_init`          | Initialization of x for all stages                                                            | `sum(nx_i), i=0,..., N`            | Yes          |
| `u_init`          | Initialization of u for stages 0 to N-1                                                       | `sum(nu_i), i=0,..., N-1`            | Yes          |
| `pi_init`         | Initialization of pi for stages 0 to N-1                                                      | `sum(npi_i), i=0,..., N-1`           | Yes          |
| `slacks_init`     | Initialization of slack values for all stages (0 to N)                                        | `2 * ns_total`        | Yes          |
| `rti_phase`       | Real-time iteration phase                                                                     | `1`                   | Yes          |


### List of possible outputs
This is a list of possible outputs of the Simulink block of an OCP solver which can be activated by setting the corresponding values in the acados Simulink options.

| Port Name        | Description                                                                                 | Shape                           | MOCP Support |
|------------------|---------------------------------------------------------------------------------------------|----------------------------------|--------------|
| `u0`             | Control input at node 0                                                                     | `nu_0`                          | Yes          |
| `utraj`          | Control input concatenated for nodes 0 to N-1                                                | `sum(nu_i), i=0,..., N-1`        | Yes          |
| `xtraj`          | State concatenated for nodes 0 to N                                                         | `sum(nx_i), i=0,..., N`          | Yes          |
| `ztraj`          | Algebraic states concatenated for nodes 0 to N-1                                             | `sum(nz_i), i=0,..., N-1`        | Yes          |
| `pi_all`         | Equality Lagrange multipliers concatenated for nodes 0 to N-1                                | `sum(npi_i), i=0,..., N-1`       | Yes          |
| `slack_values`   | Slack values concatenated in order [sl_0, su_0, ..., sl_N, su_N]                             | `2 * ns_total`                   | Yes          |
| `solver_status`  | Acados solver status (0 = SUCCESS)                                                          | `1`                              | Yes          |
| `cost_value`     | Cost function value                                                                         | `1`                              | Yes          |
| `KKT_residual`   | KKT residual                                                                                | `1`                              | Yes          |
| `KKT_residuals`  | KKT residuals, size [4] (stat, eq, ineq, comp)                                               | `4`                              | Yes          |
| `x1`             | State at node 1                                                                             | `nx_1`                           | Yes          |
| `CPU_time`       | CPU time                                                                                    | `1`                              | Yes          |
| `CPU_time_sim`   | CPU time for integrator                                                                     | `1`                              | Yes          |
| `CPU_time_qp`    | CPU time for QP solution                                                                    | `1`                              | Yes          |
| `CPU_time_lin`   | CPU time for linearization (including integrator)                                            | `1`                              | Yes          |
| `sqp_iter`       | NLP solver iterations                                                                       | `1`                              | Yes          |
| `parameter_traj` | Parameter trajectory                                                                        | `sum(np_i), i=0,..., N`          | Yes          |



### Developing extensions
If you want a more advanced interaction with the `acados` solver via Simulink, feel free to edit the corresponding templates in [`<acados_root>/interfaces/acados_template/acados_template/c_templates_tera/matlab_templates`](https://github.com/acados/acados/tree/main/interfaces/acados_template/acados_template/c_templates_tera/matlab_templates) to add more inputs or outputs.

### S-function Mask
The Simulink S-functions interface can be improved using the following mask commands, which shows the names of input and output ports, facilitating debugging operations.
These commands use global variables that are automatically generated by the `make_sfun` command.
1. OCP S-function mask
```
global sfun_input_names sfun_output_names
for i = 1:length(sfun_input_names)
	port_label('input', i, sfun_input_names{i})
end
for i = 1:length(sfun_output_names)
	port_label('output', i, sfun_output_names{i})
end
```
2. SIM S-function mask
```
global sfun_sim_input_names sfun_sim_output_names
for i = 1:length(sfun_sim_input_names)
	port_label('input', i, sfun_sim_input_names{i})
end
for i = 1:length(sfun_sim_output_names)
	port_label('output', i, sfun_sim_output_names{i})
end
```
To use the mask command just copy-paste it in the "icon drawing commands" field, accessible by right clicking on the S-function block - Mask - Edit Mask.

## Setup CasADi
To create external function for your problem, we suggest to use `CasADi` from the folder `<acados_root_folder>/external`.
Depending on the environment you want to use to generate `CasADi` functions from, proceed with the corresponding paragraph (MATLAB, Octave).

Any CasADi version between 3.4.0 and 3.6.7 should work.
If you don't have CasADi yet, you can install it as described below.

### **MATLAB**
Download and extract the `CasADi` binaries into `<acados_root_folder>/external/casadi-matlab`:
```
cd external
wget -q -nc --show-progress https://github.com/casadi/casadi/releases/download/3.4.0/casadi-linux-matlabR2014b-v3.4.0.tar.gz
mkdir -p casadi-matlab
tar -xf casadi-linux-matlabR2014b-v3.4.0.tar.gz -C casadi-matlab
cd ..
```

### **Octave version 6.2 or later**
Download and extract the `CasADi` binaries into `<acados_root_folder>/external/casadi-octave`:
```
cd external
wget -O casadi-linux-octave.zip https://github.com/casadi/casadi/releases/download/3.6.7/casadi-3.6.7-linux64-octave7.3.0.zip
mkdir -p casadi-octave
unzip casadi-linux-octave.zip -d ./casadi-octave;
```


<!-- ### **Octave version 4.4 or later**
Download and extract the `CasADi` binaries into `<acados_root_folder>/external/casadi-octave`:
```
cd external
wget -q -nc --show-progress https://github.com/casadi/casadi/releases/download/3.4.5/casadi-linux-octave-4.4.1-v3.4.5.tar.gz
mkdir -p casadi-octave
tar -xf casadi-linux-octave-4.4.1-v3.4.5.tar.gz -C casadi-octave
```

### **Octave version 4.2 or earlier**
Download and extract the `CasADi` binaries into `<acados_root_folder>/external/casadi-octave`:
```
cd external
wget -q -nc --show-progress https://github.com/casadi/casadi/releases/download/3.4.0/casadi-linux-octave-v3.4.0.tar.gz
mkdir -p casadi-octave
tar -xf casadi-linux-octave-v3.4.0.tar.gz -C casadi-octave
cd ..
``` -->
