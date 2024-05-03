# acados

<!-- ![](https://secure.travis-ci.org/acados/acados.png?branch=master) -->

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


<!-- ![Github actions full build workflow](https://github.com/acados/acados/actions/workflows/full_build.yml/badge.svg?branch=master) -->
<!-- ![](https://ci.appveyor.com/api/projects/status/q0b2nohk476u5clg?svg=true) -->

Fast and embedded solvers for nonlinear optimal control.

- `acados` __source code__ is hosted on [Github](https://github.com/acados/acados).
Contributions via Pull requests are welcome!

- `acados` has a discourse based [__forum__](https://discourse.acados.org/).

- `acados` is mainly developed by the group around Prof. Moritz Diehl, the Systems Control and Optimization Laboratory (__syscop__), at the University of Freiburg. [More infos on the syscop web page](https://www.syscop.de/).


# About `acados`

`acados` is a software package for the efficient solution of optimal control and estimation problems.
<!-- It is the successor of the [`ACADO`](https://acado.github.io/) software package developed at KU Leuven and University of Freiburg by the team of Prof. Moritz Diehl. -->
It provides a collection of computationally efficient building blocks tailored to optimal control and estimation problems.
Among others, it implements:
- modules for the integration of ordinary differential equations (ODE) and differential-algebraic equations (DAE),
- interfaces to state-of-the-art QP solvers like [`HPIPM`](https://github.com/giaf/hpipm), `qpOASES`, [`DAQP`](https://github.com/darnstrom/daqp) and [`OSQP`](https://github.com/oxfordcontrol/osqp),
- (partial) condensing routines
- nonlinear programming solvers for optimal control structured problems
- real-time algorithms, such as the real-time iteration (RTI) and Advanced-Step Real-time iteration (AS-RTI) algorithms
The back-end of acados uses the high-performance linear algebra package [`BLASFEO`](https://github.com/giaf/blasfeo), in order to boost computational efficiency for small to medium scale matrices typical of embedded optimization applications.
`MATLAB`, `Octave` and `Python` interfaces can be used to conveniently describe optimal control problems and generate self-contained C code that can be readily deployed on embedded platforms.

## Further Reading
- or other papers stated under [Citing](https://docs.acados.org/citing.html) and


# Citing

## First journal publication on the `acados` software framework:
```latex
@Article{Verschueren2021,
  Title                    = {acados -- a modular open-source framework for fast embedded optimal control},
  Author                   = {Robin Verschueren and Gianluca Frison and Dimitris Kouzoupis and Jonathan Frey and Niels van Duijkeren and Andrea Zanelli and Branimir Novoselnik and Thivaharan Albin and Rien Quirynen and Moritz Diehl},
  Journal                  = {Mathematical Programming Computation},
  Year                     = {2021},
  Doi                      = {10.1007/s12532-021-00208-8},
  ISSN                     = {1867-2957},
  url = {https://doi.org/10.1007/s12532-021-00208-8},
}
```

## Publications on advanced `acados` features:


### Efficient Zero-Order Robust Optimization (zoRO) for Real-Time Model Predictive Control with acados
```latex
@InProceedings{Frey2024,
  Title                    = {Efficient Zero-Order Robust Optimization for Real-Time Model Predictive Control with acados},
  Author                   = {Jonathan Frey and Yunfan Gao and Florian Messerer and Amon Lahr and Melanie N Zeilinger and Moritz Diehl},
  Booktitle                = ECC,
  Year                     = {2024},
}
```


### Fast integrators with sensitivity propagation for use
This paper demonstrates the efficiency of `acados` integrators, compared with integrators available in CasADi.
Additionally it describes the [`casados-integrator`](https://github.com/FreyJo/casados-integrators) a wrapper which allows using the `acados` integrators within a `CasADi` NLP solver, like IPOPT.

```latex
@InProceedings{Frey2023,
  Title                    = {Fast integrators with sensitivity propagation for use in {C}as{AD}i},
  Author                   = {Frey, Jonathan and De Schutter, Jochem and Diehl, Moritz},
  Booktitle                = ECC,
  Year                     = {2023},
}
```

### Gauss-Newton Runge-Kutta (GNRK) integrators for efficient discretization of OCPs with long horizons and least-squares costs

Can be used by setting: the option `cost_discretization = 'INTEGRATOR'`.

```latex
@InProceedings{Frey2023b,
  Title                    = {{G}auss-{N}ewton {R}unge-{K}utta Integration for Efficient Discretization of Optimal Control Problems with Long Horizons and Least-Squares Costs},
  Author                   = {Jonathan Frey and Katrin Baumgärtner and Moritz Diehl},
  Booktitle                = {accepted for ECC 2024},
  Url                      = {https://arxiv.org/abs/2310.00618}
}
```

### Structure exploiting implicit Runge-Kutta method: GNSF
Can be used by setting: the option `integrator_tpye = 'GNSF'`.
```latex
@InProceedings{Frey2019,
  Title                    = {Detecting and Exploiting {G}eneralized {N}onlinear {S}tatic {F}eedback
Structures in {DAE} Systems for {MPC}},
  Author                   = {Jonathan Frey and Rien Quirynen and Dimitris Kouzoupis and Gianluca Frison and
Jens Geisler and Axel Schild and Moritz Diehl},
  Booktitle                = ECC,
  Year                     = {2019},
}
```



<!-- ## First acados publication in Proceedings of the IFAC Conference
```latex
@inproceedings{Verschueren2018,
    booktitle = {Proceedings of the IFAC Conference on Nonlinear Model Predictive Control (NMPC)},
    title = {Towards a modular software package for embedded optimization},
    year = {2018},
    author = {Robin Verschueren and Gianluca Frison and Dimitris Kouzoupis and Niels van Duijkeren and Andrea Zanelli and Rien Quirynen and Moritz Diehl},
}
``` -->

# Documentation page overview

```eval_rst
Documentation latest build: |today|
```

<!--
```eval_rst
.. toctree::
    Home<self>
```
 -->

```eval_rst
.. toctree::
   :maxdepth: 2
   :caption: Contents

   problem_formulation/index
   algorithm_overview/index
   installation/index
   interfaces/index
   python_interface/index
   matlab_octave_interface/index
   c_interface/index
   examples/index
   faq/index
   embedded_workflow/index
   list_of_projects/index
   developer_guide/index
```

