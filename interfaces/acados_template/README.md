acados_template is an (experimental) Python package that can be used to specify optimal control problems from Python and to generate self-contained C code that uses the acados solvers to solve them. The framework is based on templated C files which are rendered from Python using the templating engine Jinja2. Notice that, at the moment, some of the features are not yet implemented (see below). 
## problem formulation 

![alt text](https://github.com/zanellia/acados/blob/master/interfaces/acados_template/docs/acados_template_docs-crop.png)

## usage
You can check out the examples folder to learn about  how to use acados_template. First of all, you will need to install it by running 'pip install .' from the root folder. The you should be able to import it as a Python module and specify the problem formulation as in examples/<example_name>/generate_c_code.py
