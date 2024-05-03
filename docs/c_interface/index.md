# C Interface

The `C` Interface of `acados` is an efficient interface to the core functionalities of `acados`.
It provides setters and getters that can be used to interact with the core of `acados` with negligible computational overhead.
Setters and getters for matrices use column-major format.

If you want to use `acados` directly from `C` or `C++`, it is recommended to only use the functions in in [`interfaces/acados_c/`](https://github.com/acados/acados/tree/master/interfaces/acados_c), which encapsulate the `acados` core, i.e. the functionality implemented in [`acados/`](https://github.com/acados/acados/tree/master/acados).

Disclaimer: the `C` interface is NOT thoroughly documented using docstrings.
It is recommended to instead look at the header files in [`interfaces/acados_c/`](https://github.com/acados/acados/tree/master/interfaces/acados_c).


## Examples
Note that the `Matlab` and `Python` interfaces can be used to code generate `C` examples which are cleaner than some of the plain C examples in [`examples/c`](https://github.com/acados/acados/tree/master/examples/c).
A recommended workflow is thus to prototype an NMPC controller from one of the high-level interfaces and deploy the generated code with minor modifications in a `C`, `C++` or `ROS` framework.

A very important resource are the templated `C` files used by the template based interfaces (`Python` and code generation for `Matlab` and `Octave`), which show how to use the `C` interface properly.
These templates are actively maintained and tested using CI.
The templates can be found in [`interfaces/acados_template/acados_template/c_templates_tera`](https://github.com/acados/acados/tree/master/interfaces/acados_template/acados_template/c_templates_tera).

## docstring based documentation
```eval_rst
.. toctree::
   :maxdepth: 2

   dense_qp_interface
   ocp_qp_interface
   sim_interface
   ocp_nlp_interface
   condensing_interface
   external_function_interface
```

## Indices and tables

```eval_rst
* :ref:`genindex`
```
