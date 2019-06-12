`acados_template` is a Python package that can be used to specify optimal control problems from Python and to generate self-contained C code that uses the acados solvers to solve them. The framework is based on templated C files which are rendered from Python using the templating engine Jinja2. Notice that, at the moment, some of the features are not yet implemented (see below). 
## problem formulation 

![alt text](https://github.com/zanellia/acados/blob/master/interfaces/acados_template/docs/acados_template_docs-crop.png)

## usage
### Linux/macOs 
You can check out the examples folder to learn about  how to use acados_template. First of all, you need to compile and install acados without the qpDUNES, HPMPC and OSQP QP solvers running 
~~~
cmake -DACADOS_WITH_QPDUNES=OFF -DACADOS_WITH_HPMPC=OFF -DACADOS_WITH_OSQP=OFF .. & make install
~~~
Then, you will need to install acados_template Python package by running 'pip install .' from the root folder. You should now be able to import it as a Python module and specify the problem formulation as in examples/<example_name>/generate_c_code.py
### Windows
You can in principle install the acados_template package within your native Python shell, but we highly recommend 
using Windows Subsystems for Linux (https://docs.microsoft.com/en-us/windows/wsl/about) and to follow the 
Linux/macOS installation instruction.
