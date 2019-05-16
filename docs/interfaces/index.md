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
.. automodule:: acados_template.ocp_nlp_render_arguments
    :members:
    :private-members:
    :exclude-members: acados_ocp2json_layout, cast_ocp_nlp, dict2json_layout, dict2json_layout_rec, check_ra, json2dict_rec

```
## acados embedded - MATLAB
