# FAQ

A few frequently asked questions are answered here.

For further questions, we refer to the [`acados` forum](https://discourse.acados.org/)

## What is the difference between `acados` and `ACADO`?

`ACADO` was heavily based on code-generation, in `acados` code-generation
is used only for problem function derivatives.

In order to maximize performance `acados` is using [`BLASFEO`](https://blasfeo.syscop.de/), a basic linear algebra implementation with hand-optimized kernels for different CPU architectures.


## What is the difference between `acados` and `CasADi`?
Completely different.

`CasADi` is typically used by `acados` as a front-end to state nonlinear problem functions needed formulate optimal control problems (OCP) and Moving Horizon Estimation problems (MHE).
