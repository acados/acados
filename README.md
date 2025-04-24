<!-- # acados -->

![](docs/_static/acados_logo.png)
<!-- [![Travis Status](https://secure.travis-ci.org/acados/acados.png?branch=master)](http://travis-ci.org/acados/acados) -->
[![Appveyor status](https://ci.appveyor.com/api/projects/status/q0b2nohk476u5clg?svg=true)](https://ci.appveyor.com/project/roversch/acados)
![Github actions full build workflow](https://github.com/acados/acados/actions/workflows/full_build.yml/badge.svg)
<!-- [![codecov](https://codecov.io/gh/acados/acados/branch/master/graph/badge.svg)](https://codecov.io/gh/acados/acados) -->

`acados` provides fast and embedded solvers for nonlinear optimal control, specifically designed for real-time applications and embedded systems.
It is written in `C` and offers interfaces to the programming languages `Python`, `MATLAB` and `Octave`.

## General
`acados` is a modular and efficient software package for solving nonlinear programs (NLP) with an optimal control problem (OCP) structure.
Such problems have to be solved repeatedly in **model predictive control (MPC)** and **moving horizon estimation (MHE)**.
The computational efficiency and modularity make `acados` an ideal choice for real-time applications.
It is designed for high-performance applications, embedded computations, and has been successfully used in [a wide range of applications](#fields-of-applications).

### Key Features:
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

## Documentation
- Documentation can be found on [docs.acados.org](https://docs.acados.org/)
- An overview of the interfaces can be found at [docs.acados.org/interfaces](https://docs.acados.org/interfaces)

## Forum
- If you have any `acados`-related questions, feel free to post on our forum at [discourse.acados.org](https://discourse.acados.org/)

## Citing
- References can be found at [docs.acados.org/citing](https://docs.acados.org/citing)

## Installation
- Instructions can be found at
[docs.acados.org/installation](https://docs.acados.org/installation)

### Design paradigms
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
