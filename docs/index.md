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

`acados` is a modular and efficient software package for solving nonlinear programs (NLP) with an optimal control problem (OCP) structure.
Such problems have to be solved repeatedly in **model predictive control (MPC)** and **moving horizon estimation (MHE)**.
The computational efficiency and modularity make `acados` an ideal choice for real-time applications.
It is designed for high-performance applications, embedded computations, and has been successfully used in [a wide range of applications](#fields-of-applications).

`acados` is written in `C`, but control problems can be conveniently formulated using the [`CasADi`](https://web.casadi.org/) symbolic framework via the high-level `acados` interfaces to the programming languages `Python`, `MATLAB` and `Octave`.

Some key features of `acados` are summarized in the following.
The [software design](#design-paradigms) allows to implement many algorithms beyond this list.
- **Nonlinear and economic model predictive control (NMPC)**: Solve challenging control problems with nonlinear dynamics and cost functions.
- **Moving horizon estimation (MHE)**: Estimate states and parameters of dynamic systems in real-time.
- **Support for differential algebraic equations (DAE)**: Efficiently handle systems with algebraic constraints.
- **Multiple shooting method**: Leverage the multiple shooting approach for time discretization, enabling fast and robust solutions.
- **Efficient integration methods**: Include advanced integrators for solving ODEs and DAEs, with support for first- and second-order sensitivities.
- **Real-time performance**: Optimized for high-frequency control loops, enabling reliable solutions for time-critical applications.
- **High-performance solvers**: Implement fast SQP-type solvers tailored for optimal control problems.
- **Modular design**: Easily extend and combine components for simulation, estimation, and control to fit diverse applications.
- **Solution sensitivity computation and combination with reinforcement learning (RL)**: The combination of MPC and RL is a hot research topic in control. Many learning algorithms can profit from the availability of solution sensitivities or in particular policy gradients.
`acados` offers the possibility to embed an NLP solver as a differentiable layer in an ML architecture as is demonstrated in the [`leap-c` project](https://github.com/leap-c/leap-c).

The back-end of acados uses the high-performance linear algebra package [`BLASFEO`](https://github.com/giaf/blasfeo), in order to boost computational efficiency for small to medium scale matrices typical of embedded optimization applications.
`MATLAB`, `Octave` and `Python` interfaces can be used to conveniently describe optimal control problems and generate self-contained C code that can be readily deployed on embedded platforms.


### Design paradigms:
The main design paradigms of `acados` are
- **efficiency**: realized by rigorously exploiting the OCP structure via tailored quadratic programming (QP) solvers, such as `HPIPM`, and (partial) condensing methods to transform QPs, enabling their efficient treatment.
Moreover, the common structure of slack variables, which for example occur when formulating soft constraints, can be exploited.
Additionally, a structure exploiting Runge-Kutta method is implemented, allowing to utilize linear dependencies within dynamical system models.
- **modularity**:
`acados` offers an extremely flexible problem formulation, allowing to not only formulate problems which occur in MPC and MHE.
More precisely, all problem functions and dimensions can vary between all stages.
Such problems are often called *multi-stage* or *multi-phase* problems.
Different NLP solvers, QP solvers, integration methods, regularization methods and globalization methods can be combined freely.
Moreover, cost and constraint functions can be declared by explicitly providing general *convex-over-nonlinear* structures, which can be exploited in the solvers.
- **usability**: The interfaces to Python, MATLAB, Simulink and Octave allow users to conveniently specify their problem in different domains and to specify their nonlinear expressions via the popular [`CasADi`](https://web.casadi.org/) symbolic software framework.
The interfaces allow to conveniently specify commonly used problem formulations via the `AcadosOcp` class and additionally expose the full flexibility of the internal `acados` problem formulation, via multi-phase formulations and `AcadosMultiphaseOcp`.

## Fields of applications
A non-exhaustive list of projects featuring `acados` is available at [docs.acados.org/list_of_projects](https://docs.acados.org/list_of_projects/index.html).
Contributions to this list are very welcome and allow to increase visibility of your work among other `acados` users.
- Robotics: Real-time NMPC for quadrotors, legged locomotion, and agile robotic platforms.
- Autonomous Vehicles: Used in projects like openpilot in driving assistance systems.
- Energy Systems: Optimization-based control for microgrids and wind turbines.
- Biomechanics: Optimal control in biomechanics through libraries like bioptim.
- Aerospace: Applications in trajectory optimization and control for drones and morphing-wing aircraft.


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

.. toctree::
    :maxdepth: 2
    :caption: Interfaces

    interfaces/index
    python_interface/index
    matlab_octave_interface/index
    embedded_workflow/index

.. toctree::
    :maxdepth: 2
    :caption: User Guide

    problem_formulation/index
    troubleshooting/index
    features/index
```

<!-- c_interface/index -->