# Getting Started

`acados` is a software package for the efficient solution of
optimal control and estimation problems. It is the successor
of the [`ACADO`](https://acado.github.io/) software package
developed at KU Leuven and University of Freiburg by the team
of Prof. Moritz Diehl. It provides a collection of computationally
efficient building blocks tailored to optimal control and estimation
problems. Among others, it implements: modules for the integration
of ordinary differential equations and differential-algebraic equations,
interfaces to state-of-the-art QP solvers like qpOASES, HPIPM, qpDUNES
and `OSQP`, condensing routines and nonlinear programming solvers
based on the real-time iteration framework. The back-end of acados
uses the high-performance linear algebra package BLASFEO, in order
to boost computational efficiency for small to medium scale matrices
typical of embedded optimization applications. MATLAB/Octave and Python interfaces
can be used to conveniently describe optimal control problems and generate self-contained 
C code that can be readily deployed on embedded platforms.

## Problem formulation

The problem formulation targeted by acados is described [here](https://github.com/acados/acados/blob/master/docs/problem_formulation/problem_formulation_ocp_mex.pdf). 
