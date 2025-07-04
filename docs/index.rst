acados
======

.. meta::
   :description: acados is an open-source software package for fast and embedded nonlinear model predictive control (MPC) and moving horizon estimation (MHE). It is written in C and has interfaces to Python, MATLAB, Octave, and Simulink. The software is designed for real-time applications and is used in various fields such as robotics, autonomous vehicles, energy systems, biomechanics, and aerospace.
   :keywords: optimal control, nonlinear programming, nonlinear model predictive control, embedded optimization, real-time optimization, moving horizon estimation, open-source software, software library, C library, Python interface, MATLAB interface, Octave interface, Simulink
   :google-site-verification: otz_joxGlxsY8wbpkkqOty47dLxqqWPa9JCC5HeyCg8
.. |github-workflow-full-build| image:: https://github.com/acados/acados/actions/workflows/full_build.yml/badge.svg
   :target: https://github.com/acados/acados/actions/workflows/full_build.yml
   :alt: Github workflow status

.. |github-workflow-c_test_blasfeo_reference| image:: https://github.com/acados/acados/actions/workflows/c_test_blasfeo_reference.yml/badge.svg
   :target: https://github.com/acados/acados/actions/workflows/c_test_blasfeo_reference.yml
   :alt: Github workflow status

.. |appveyor-build| image:: https://ci.appveyor.com/api/projects/status/q0b2nohk476u5clg?svg=true
   :target: https://ci.appveyor.com/project/roversch/acados
   :alt: Appveyor workflow status

|github-workflow-full-build|
|github-workflow-c_test_blasfeo_reference|
|appveyor-build|

Fast and embedded solvers for real-world applications of nonlinear optimal control.


Important links
----------------
|:cinema:| Get inspired by `real-world applications using acados <https://docs.acados.org/real_world_examples/index.html>`_ |:rocket:|

|:star:| The ``acados`` **source code** is hosted on `Github <https://github.com/acados/acados>`_.
Contributions via pull requests are welcome!

|:handshake:| ``acados`` has a discourse-based `forum <https://discourse.acados.org/>`_.

|:homes:| ``acados`` is mainly developed by the `syscop group around Prof. Moritz Diehl, at the University of Freiburg <https://www.syscop.de/>`_.

About ``acados``
----------------

``acados`` is a modular and efficient software package for solving nonlinear programs (NLP) with an optimal control problem (OCP) structure.
Such problems have to be solved repeatedly in **model predictive control (MPC)** and **moving horizon estimation (MHE)**.
The computational efficiency and modularity make ``acados`` an ideal choice for real-time applications.
It is designed for high-performance applications, embedded computations, and has been successfully used in `a wide range of applications <https://docs.acados.org/list_of_projects/index.html>`_.

``acados`` is written in ``C``, but control problems can be conveniently formulated using the `CasADi <https://web.casadi.org/>`_ symbolic framework via the high-level ``acados`` interfaces to the programming languages ``Python``, ``MATLAB``, and ``Octave``.

Some key features of ``acados`` are summarized in the following.
The `software design <#design-paradigms>`_ allows implementing many algorithms beyond this list:

- **Nonlinear and economic model predictive control (NMPC)**: Solve challenging control problems with nonlinear dynamics and cost functions.
- **Moving horizon estimation (MHE)**: Estimate states and parameters of dynamic systems in real-time.
- **Support for differential algebraic equations (DAE)**: Efficiently handle systems with algebraic constraints.
- **Multiple shooting method**: Leverage the multiple shooting approach for time discretization, enabling fast and robust solutions.
- **Efficient integration methods**: Include advanced integrators for solving ODEs and DAEs, with support for first- and second-order sensitivities.
- **Real-time performance**: Optimized for high-frequency control loops, enabling reliable solutions for time-critical applications.
- **High-performance solvers**: Implement fast SQP-type solvers tailored for optimal control problems.
- **Modular design**: Easily extend and combine components for simulation, estimation, and control to fit diverse applications.
- **Solution sensitivity computation and combination with reinforcement learning (RL)**: The combination of MPC and RL is a hot research topic in control. Many learning algorithms can profit from the availability of solution sensitivities or, in particular, policy gradients.
  ``acados`` offers the possibility to embed an NLP solver as a differentiable layer in an ML architecture, as demonstrated in the `leap-c project <https://github.com/leap-c/leap-c>`_.

The back-end of ``acados`` uses the high-performance linear algebra package `BLASFEO <https://github.com/giaf/blasfeo>`_, in order to boost computational efficiency for small to medium-scale matrices typical of embedded optimization applications.
``MATLAB``, ``Octave``, and ``Python`` interfaces can be used to conveniently describe optimal control problems and generate self-contained C code that can be readily deployed on embedded platforms.

Design paradigms
----------------

The main design paradigms of ``acados`` are:

- **Efficiency**: Realized by rigorously exploiting the OCP structure via tailored quadratic programming (QP) solvers, such as ``HPIPM``, and (partial) condensing methods to transform QPs, enabling their efficient treatment.
  Moreover, the common structure of slack variables, which, for example, occur when formulating soft constraints, can be exploited.
  Additionally, a structure-exploiting Runge-Kutta method is implemented, allowing the utilization of linear dependencies within dynamical system models.
- **Modularity**: ``acados`` offers an extremely flexible problem formulation, allowing not only the formulation of problems that occur in MPC and MHE.
  More precisely, all problem functions and dimensions can vary between all stages.
  Such problems are often called *multi-stage* or *multi-phase* problems.
  Different NLP solvers, QP solvers, integration methods, regularization methods, and globalization methods can be combined freely.
  Moreover, cost and constraint functions can be declared by explicitly providing general *convex-over-nonlinear* structures, which can be exploited in the solvers.
- **Usability**: The interfaces to Python, MATLAB, Simulink, and Octave allow users to conveniently specify their problem in different domains and to specify their nonlinear expressions via the popular `CasADi <https://web.casadi.org/>`_ symbolic software framework.
  The interfaces allow users to conveniently specify commonly used problem formulations via the ``AcadosOcp`` class and additionally expose the full flexibility of the internal ``acados`` problem formulation, via multi-phase formulations and ``AcadosMultiphaseOcp``.

Fields of applications
----------------------

``acados`` is used in a wide range of applications. Some examples include:
- **Robotics**: Real-time NMPC for quadrotors, legged locomotion, and agile robotic platforms.
- **Autonomous Vehicles**: Used in projects like openpilot in driving assistance systems.
- **Energy Systems**: Optimization-based control for microgrids and wind turbines.
- **Biomechanics**: Optimal control in biomechanics through libraries like bioptim.
- **Aerospace**: Applications in trajectory optimization and control for drones and morphing-wing aircraft.

Please also check this non-exhaustive `list of projects <https://docs.acados.org/list_of_projects/index.html>`_ featuring ``acados``.
Contributions to this list are very welcome and allow increasing the visibility of your work among other ``acados`` users.

Documentation page overview
---------------------------

Documentation latest build: |today|

.. toctree::
   :maxdepth: 2

   Home <self>
   real_world_examples/index
   citing/index
   installation/index
   list_of_projects/index

.. toctree::
   :maxdepth: 2
   :caption: Interfaces

   interfaces/index
   python_interface/index
   matlab_octave_interface/index

.. toctree::
   :maxdepth: 2
   :caption: User Guide

   problem_formulation/index
   troubleshooting/index
   features/index

.. toctree::
   :maxdepth: 2
   :caption: Advanced

   developer_guide/index
   embedded_workflow/index
