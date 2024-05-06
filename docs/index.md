# acados


```eval_rst

|github-workflow-full-build|
|github-workflow-c_test_blasfeo_reference|
|appveyor-build|

.. |github-workflow-full-build| image:: https://github.com/acados/acados/actions/workflows/full_build.yml/badge.svg
    :target: https://github.com/acados/acados/actions/workflows/full_build.yml
    :alt: Github workflow status

.. |github-workflow-c_test_blasfeo_reference| image:: https://github.com/acados/acados/actions/workflows/c_test_blasfeo_reference.yml/badge.svg
    :target: https://github.com/acados/acados/actions/workflows/c_test_blasfeo_reference.yml
    :alt: Github workflow status


.. |appveyor-build| image:: https://ci.appveyor.com/api/projects/status/q0b2nohk476u5clg?svg=true
    :target: https://ci.appveyor.com/project/roversch/acados
    :alt: Appveyor workflow status

```

Fast and embedded solvers for nonlinear optimal control.

- `acados` source code is hosted on [Github](https://github.com/acados/acados).
Contributions via pull requests are welcome!
- `acados` has a discourse based [forum](https://discourse.acados.org/).
- `acados` is mainly developed by the [syscop group around Prof. Moritz Diehl, the Systems Control and Optimization Laboratory, at the University of Freiburg](https://www.syscop.de/).


## About `acados`

`acados` is a software package providing fast and embedded solvers for nonlinear optimal control.
Problems can be conveniently formulated using the [`CasADi`](https://web.casadi.org/) symbolic framework and the high-level `acados` interfaces.

`acados` provides a collection of computationally efficient building blocks tailored to optimal control structured problems, most prominently optimal control problems (OCP) and moving horizon estimation (MHE) problems.
Among others, `acados` implements:
- modules for the integration of ordinary differential equations (ODE) and differential-algebraic equations (DAE),
- interfaces to state-of-the-art QP solvers like [`HPIPM`](https://github.com/giaf/hpipm), `qpOASES`, [`DAQP`](https://github.com/darnstrom/daqp) and [`OSQP`](https://github.com/oxfordcontrol/osqp)
- (partial) condensing routines, provided by `HPIPM`
- nonlinear programming solvers for optimal control structured problems
- real-time algorithms, such as the real-time iteration (RTI) and advanced-step real-time iteration (AS-RTI) algorithms

The back-end of acados uses the high-performance linear algebra package [`BLASFEO`](https://github.com/giaf/blasfeo), in order to boost computational efficiency for small to medium scale matrices typical of embedded optimization applications.
`MATLAB`, `Octave` and `Python` interfaces can be used to conveniently describe optimal control problems and generate self-contained C code that can be readily deployed on embedded platforms.


## Problem Formulation

Since `acados` mainly aims on providing SQP type methods for optimal control, it naturally needs optimal control structured nonlinear programming formulations (OCP-NLP) and quadratic programming (QP) formulations to tackle the subproblems within SQP.

- __Optimal control structured NLP (OCP-NLP)__: The problem formulation targeted by `acados` OCP solver is stated [here](https://github.com/acados/acados/blob/master/docs/problem_formulation/problem_formulation_ocp_mex.pdf).


- __QP formulations (dense and OCP structured)__: `acados` relies on `HPIPM` for reformulating QP problems via (partial) condensing and expansion routines.
We thus use the flexible QP formulations from `HPIPM` for optimal control structured quadratic programming formulation (OCP-QP) and the dense QP formulation.
Both problem formulations are documented in the [`HPIPM guide`](https://github.com/giaf/hpipm/blob/master/doc/guide.pdf).


# Documentation page overview

```eval_rst
Documentation latest build: |today|
```

```eval_rst
.. toctree::
    :maxdepth: 2

    Home<self>
    citing/index
    installation/index
    list_of_projects/index
    developer_guide/index
```


```eval_rst
.. toctree::
    :maxdepth: 2
    :caption: Interfaces

    interfaces/index
    python_interface/index
    matlab_octave_interface/index
    c_interface/index
    embedded_workflow/index
```
