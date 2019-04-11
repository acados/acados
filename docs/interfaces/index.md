# Interfaces


## C API
The C API of acados is an efficient interface to the core functionalities of acados. 
It provides setters and getters that can be used to interact with the core of acados with 
negligible computation overhead. In order to learn about the acados C API, you 
can look at the examples in
[`acados/examples/c/`](https://github.com/acados/acados/tree/master/examples/c). 


## acados MEX

This interface makes a subset of `acados` functionalities available from Matlab or Octave.

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



## acados Python


`acados_template` is an (experimental) Python package that can be used to specify optimal control problems from Python and to generate self-contained C code that uses the acados solvers to solve them.

The framework is based on templated C files which are rendered from Python using the templating engine `Jinja2`. Notice that, at the moment, some of the features are not yet implemented.


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



