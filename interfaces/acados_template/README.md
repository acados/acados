`acados_template` is a Python package that can be used to specify optimal control problems from Python and to generate self-contained C code that uses the acados solvers to solve them.

## usage
### Linux/macOs 
You can check out the examples folder to learn about  how to use acados_template. First of all, you need to compile and install acados with the qpOASES solver running 
~~~
cmake -DACADOS_WITH_QPOASES=ON .. & make install
~~~
Then, you will need to install acados_template Python package by running 'pip3 install .' from the Python package root folder '<acados_root/interfaces/acados_template>'. You should now be able to import it as a Python module and specify the problem formulation as in examples/<example_name>/generate_c_code.py. Addiotionally, you will have to make sure that the environment variable `LD_LIBRARY_PATH` contains the path to `libacados.so` (default path is `<acados_root/lib>`).
### Windows
You can in principle install the acados_template package within your native Python shell, but we highly recommend 
using Windows Subsystems for Linux (https://docs.microsoft.com/en-us/windows/wsl/about) and to follow the 
Linux/macOS installation instruction.

For more information visit
https://docs.acados.org/interfaces/index.html#acados-embedded-python
