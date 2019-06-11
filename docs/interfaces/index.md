# Interfaces


## C API
The C API of acados is an efficient interface to the core functionalities of acados. 
It provides setters and getters that can be used to interact with the core of acados with 
negligible computation overhead. In order to learn about the acados C API, you 
can look at the examples in
[`acados/examples/c/`](https://github.com/acados/acados/tree/master/examples/c). 


## acados MATLAB (rapid prototyping)

This interface makes a broad set of `acados` functionalities available from Matlab or Octave 
for prototyping purpose. As of now, this closely tracks the latest developments in the core of acados, e.g.
exact Hessians, adjoint corrections, regularization, etc. However, for the time being, it will not be possible to 
generate a self-contained C library that can be deployed on an embedded system. For this purpose 
see the `acados emebedded` high-level interface below. 

Some examples for the use of this interface can be found in `<acados_dir>/examples/matlab_mex`

This interface uses the shared libraries created using the make command from the main acados folder

```bash
make acados_c_shared
```

To run the examples, navigate into the selected folder, and there run the command
```
export ACADOS_INSTALL_DIR="<your acados repo dir>"
source env.sh # Which can be found in the folder of one of the examples
```

If `ACADOS_INSTALL_DIR` is not specified, it will be assumed that the examples are run from the sub-folders of the current folder (i.e. acados main folder is 2 folders up from the current folder).

Afterwards, launch Matlab or Octave from the same shell.

If you want to run the examples in a different folder, please close the current shell and open a new one to repeat the procedure: this ensures the correct setting of the environment variables.



## acados embedded - Python


`acados_template` is an Python package that can be used to specify optimal control problems from Python and to generate self-contained C code that uses the acados solvers to solve them.
In comparison with the MATLAB interface for rapid prototyping (see above), it supports less features, but it allows the user to generate a self-contained C library  
that can be easily deployed on an embedded system.

The framework is based on templated C files which are rendered from Python using the templating engine `Jinja2`.

The currently supported formulations reads as

```math
\begin{equation}
\begin{aligned}
&\underset{\begin{subarray}{c}
    x(\cdot),\,u(\cdot), \, z(\cdot)
\end{subarray}}{\min}	    &&\int_0^T l(x(\tau), u(\tau), z(\tau), p)\mathrm{d}\tau + m(x(T), z(T), p)\\ 
                            &\,\,\,\quad \text{s.t.}    &&x(0) - \bar{x}_0 = 0, &&\\
                            & 						    &&F(x(t), \dot{x}(t), u(t), z(t), p) = 0, &&\quad t \in [0,\,T),\\
                            & 						    &&\underline{h} \leq h(x(t), u(t), p) \leq \bar{h}, &&\quad t \in [0,\,T),\\
                            & 						    &&\underline{x} \leq \Pi_{x}x(t) \leq \bar{x}, &&\quad t \in [0,\,T),\\
                            & 						    &&\underline{u} \leq \Pi_{u}u(t) \leq \bar{u}, &&\quad t \in [0,\,T),\\
                            & 						    &&\underline{c} \leq Cx(t) + Du(t)\leq \bar{c}, &&\quad t \in [0,\,T), \\
                            &                           &&                                                   && \\[-1em]
                            & 						    &&\underline{h}^e \leq h^e(x(T), p) \leq \bar{h}^e, &&\\
                            & 						    &&\underline{x}^e \leq \Pi_{x}^e x(T) \leq \bar{x}^{e}, &&\\
                            & 						    &&\underline{c}^e \leq C^e x(T)\leq \bar{c}^e, &&\\
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

Currently not yet implemented features:

* :math:`l` must be in linear least-squares form :math:`l = \frac{1}{2}\| V_x x(t) + V_u u(t) + V_z z(t) - y_{\text{ref}}\|_W^2`
* :math:`m` must be in linear least-squares form :math:`m = \frac{1}{2}\| V^e_x x(t) - y_{\text{ref}}^e\|_{W^e}^2`
* Constraints cannot depend on algebraic variables (yet)
```

You can check out the examples folder to learn about  how to use acados_template.
First of all, you need to compile and install acados without the qpDUNES, HPMPC and OSQP QP solvers running:
```bash
cd <cmake_build_dir>
cmake -DACADOS_WITH_QPDUNES=OFF -DACADOS_WITH_HPMPC=OFF -DACADOS_WITH_OSQP=OFF ..
make install
```

Then, you will need to install acados_template Python package by running
```
pip install <acados_repo_root>/interfaces/acados_template
```

You should now be able to import it as a Python module and specify the problem formulation as in `examples/<example_name>/generate_c_code.py`

For more information contact `@zanellia`

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
## acados embedded - MATLAB
