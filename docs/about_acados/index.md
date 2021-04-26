# About `acados`

`acados` is a software package for the efficient solution of optimal control and estimation problems.
It is the successor of the [`ACADO`](https://acado.github.io/) software package
developed at KU Leuven and University of Freiburg by the team of Prof. Moritz Diehl.
It provides a collection of computationally efficient building blocks tailored to optimal control and estimation problems.
Among others, it implements: modules for the integration
of ordinary differential equations (ODE) and differential-algebraic equations (DAE),
interfaces to state-of-the-art QP solvers like [`HPIPM`](https://github.com/giaf/hpipm), `qpOASES`, `qpDUNES`
and [`OSQP`](https://github.com/oxfordcontrol/osqp), condensing routines and nonlinear programming solvers
based on the real-time iteration framework.
The back-end of acados uses the high-performance linear algebra package [`BLASFEO`](https://blasfeo.syscop.de/), in order
to boost computational efficiency for small to medium scale matrices
typical of embedded optimization applications.
`MATLAB`, `Octave` and `Python` interfaces can be used to conveniently describe optimal control problems and generate self-contained 
C code that can be readily deployed on embedded platforms.

## Further Reading
For more information see
- the paper [acados – a modular open-source framework for fast embedded optimal control, Verschueren 2019](https://arxiv.org/abs/1910.13753)
- or other papers stated under [Citing](https://docs.acados.org/citing.html) and
- [Papers featuring `acados`](https://docs.acados.org/list_of_projects/index.html#papers-featuring-acados).
