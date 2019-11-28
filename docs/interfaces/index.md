# Interfaces


## C API
The C API of `acados` is an efficient interface to the core functionalities of `acados`.
It provides setters and getters that can be used to interact with the core of `acados` with
negligible computation overhead. In order to learn about the `acados` C API, you
can look at the examples in
[`acados/examples/c/`](https://github.com/acados/acados/tree/master/examples/c).


## acados MATLAB (rapid prototyping)

This interface makes a broad set of `acados` functionalities available from Matlab or Octave.
As of now, this closely tracks the latest developments in the core of acados, e.g.
exact Hessians, adjoint corrections, regularization, etc.

To get started with this interface we recommend the examples in `<acados_root>/examples/acados_matlab_octave/getting_started`.

The problem formulation is stated in [this PDF](https://github.com/acados/acados/tree/master/docs/problem_formulation/problem_formulation_ocp_mex.pdf)

This interface uses the shared libraries created using the make command from the main `acados` folder

```bash
make shared_library
```

To run the examples, navigate into the selected folder, and there run the command
```
export ACADOS_INSTALL_DIR="<acados_root>"
source env.sh # Which can be found in the folder of one of the examples
```

If `ACADOS_INSTALL_DIR` is not specified, it will be assumed that the examples are run from the sub-folders of the current folder (i.e. `acados` main folder is 2 folders up from the current folder).

Afterwards, launch `Matlab` or Octave from the same shell.

If you want to run the examples in a different folder, please close the current shell and open a new one to repeat the procedure: this ensures the correct setting of the environment variables.

For a detailed explanation about the options available from this
interface see: [matlab_interface_spreadsheet](https://docs.google.com/spreadsheets/d/1rVRycLnCyaWJLwnV47u30Vokp7vRu68og3OhlDbSjDU/edit?usp=sharing) (thanks to [@EnricaSo](https://github.com/EnricaSo)).

### acados MATLAB (going embedded)
There is the option to generate embeddable `C` code from Matlab.
The workflow uses the same templates as the Python interface (see below) and the `Tera` renderer.
After creating an acados solver `ocp`, you can use the routine `ocp.generate_c_code` to generate `C` code which can be used for embedded applications.

Note: This part of the MATLAB/Octave interface does not yet support all features of the one mentioned before.


## acados embedded - Python

`acados_template` is a Python package that can be used to specify optimal control problems from Python and to generate self-contained C code that uses the acados solvers to solve them.
In comparison with the MATLAB interface for rapid prototyping (see above), it supports less features, but it allows the user to generate a self-contained C library
that can be easily deployed on an embedded system.

The framework is based on templated C files which are rendered from Python using the templating engine `Tera`.

The currently supported formulation reads as

```math
\begin{equation}
\begin{aligned}
&\underset{\begin{subarray}{c}
    x(\cdot),\,u(\cdot), \, z(\cdot)
\end{subarray}}{\min}	    &&\int_0^T l(x(\tau), u(\tau), z(\tau), p)\mathrm{d}\tau + m(x(T), z(T), p)\\ 
&\,\,\,\quad \text{s.t.}    &&x(0) - \bar{x}_0 = 0, &&\\
&&&F(x(t), \dot{x}(t), u(t), z(t), p) = 0, &&\quad t \in [0,\,T),\\
&&&\underline{h} \leq h(x(t), u(t), p) \leq \bar{h}, &&\quad t \in [0,\,T),\\
&&&\underline{x} \leq \Pi_{x}x(t) \leq \bar{x}, &&\quad t \in [0,\,T),\\
&&&\underline{u} \leq \Pi_{u}u(t) \leq \bar{u}, &&\quad t \in [0,\,T),\\
&&&\underline{c} \leq Cx(t) + Du(t)\leq \bar{c}, &&\quad t \in [0,\,T), \\
&&& && \\[-1em]
&&&\underline{h}^e \leq h^e(x(T), p) \leq \bar{h}^e, &&\\
&&&\underline{x}^e \leq \Pi_{x}^e x(T) \leq \bar{x}^{e}, &&\\
&&&\underline{c}^e \leq C^e x(T)\leq \bar{c}^e, &&\\
\end{aligned}
\end{equation}
```
```eval_rst
Where:

* :math:`l: \mathbb{R}^{n_x}\times\mathbb{R}^{n_u}\times\mathbb{R}^{n_z} \rightarrow \mathbb{R}` is the Lagrange objective term.
* :math:`m: \mathbb{R}^{n_x}\times\mathbb{R}^{n_z} \rightarrow \mathbb{R}` is the Mayer objective term.

* :math:`F: \mathbb{R}^{n_x}\times\mathbb{R}^{n_x}\times\mathbb{R}^{n_u}\times\mathbb{R}^{n_z}\times\mathbb{R}^{n_p} \rightarrow \mathbb{R}^{n_x+n_z}` describes the (potentially) fully implicit dynamics.

* :math:`h: \mathbb{R}^{n_x}\times\mathbb{R}^{n_u}\times\mathbb{R}^{n_z}\times\mathbb{R}^{n_p} \rightarrow \mathbb{R}^{n_h}` and :math:`h^e: \mathbb{R}^{n_x}\times\mathbb{R}^{n_z}\times\mathbb{R}^{n_p} \rightarrow \mathbb{R}^{n_{h_e}}` are general nonlinear functions.

* :math:`C,\,D,\,C^e,\,\Pi_x,\,\Pi_u,\,\Pi_x^e` are matrices of appropriate dimensions defining the polytopic and box constraints.

The `acados_template` interface makes some limiting assumptions on the problem formulation above, namely:
* :math:`l` must be in linear least-squares form :math:`l = \frac{1}{2}\| V_x x(t) + V_u u(t) + V_z z(t) - y_{\text{ref}}\|_W^2`
* :math:`m` must be in linear least-squares form :math:`m = \frac{1}{2}\| V^e_x x(t) - y_{\text{ref}}^e\|_{W^e}^2`
* Constraints cannot depend on algebraic variables (yet)
```


### Installation
1. Compile and install `acados` by running:
```bash
cd <acados_root>/build
cmake -DACADOS_WITH_QPOASES=ON ..
make install -j4
```

2. Install acados_template Python package by running
```
pip3 install <acados_root>/interfaces/acados_template
```

(Notice that, of course, you might need to use `pip` instead, if you run, for example, from within a Python virtual
environment) You should now be able to import it as a Python module and use it as shown in the examples in `<acados_root>/examples/acados_template/python/<example_name>/generate_c_code.py`.

In order to be able to successfully render C code templates,
you need to download the `t_renderer` binaries for your platform
from <https://github.com/acados/tera_renderer/releases/> and
place them in `<acados_root>/bin` (please strip the version and platform from the binaries (e.g.
`t_renderer-v0.0.20 -> t_renderer`). Notice that you might need to make `t_renderer` executable. Run
`export ACADOS_SOURCE_DIR=<acados_root>` such that the location of acados will be known to the Python
package at run time. Additionally, you will have to make sure that the environment variable `LD_LIBRARY_PATH` contains the path to `libacados.so` (default path is `<acados_root/lib>`). Notice that, if you want to run the examples from a location that differs from '<acados_root>/interfaces/acados_template' or you want the generated Makefile to refer to a specific path (e.g. when cross-compiling or compiling from a location different from the one where you generate the C code), you will have to adapt 'ocp.acados_include_path' and 'ocp.acados_lib_path' accordingly in the generating Python code.

For more information contact `@zanellia`.

### optimal control description
``` eval_rst
.. automodule:: acados_template.casadi_functions
    :members:
    :private-members:
    :undoc-members:
```
``` eval_rst
.. automodule:: acados_template.acados_ocp_nlp
    :members:
    :private-members:
    :exclude-members: acados_ocp2json_layout, cast_ocp_nlp, dict2json_layout, dict2json_layout_rec, check_ra, json2dict_rec

```
