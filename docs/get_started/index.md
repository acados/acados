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
typical of embedded optimization applications.

## Problem formulation

```math
\begin{equation}
\begin{aligned}
&\underset{\begin{subarray}{c}
    x(\cdot),\,u(\cdot), \, z(\cdot)
\end{subarray}}{\min}	    &&\int_0^T l(x(\tau), u(\tau), z(\tau), p)\mathrm{d}\tau + m(x(T), z(T), p)\\ 
                            &\,\,\,\quad \text{s.t.}    &&x(0) - \bar{x}_0 = 0, &&\\
                            & 						    &&\underline{x}_0 \leq \Pi_{x_0}x(0) \leq \bar{x}_0, && \\
                            & 						    &&\underline{u}_0 \leq \Pi_{u_0}u(0) \leq \bar{u}_0, && \\
                            & 						    &&\underline{z}_0 \leq \Pi_{z_0}z(0) \leq \bar{z}_0, && \\
                            & 						    &&\underline{c}_0 \leq C_0x(0) + D_0u(0) + E_0z(0)\leq \bar{c}_0, && \\
                            &                           &&                                                   && \\[-1em]
                            & 						    &&F(x(t), \dot{x}(t), u(t), z(t), p) = 0, &&\quad t \in [0,\,T),\\
                            & 						    &&\underline{h} \leq g(h(x(t), u(t), z(t), p)) \leq \bar{h}, &&\quad t \in [0,\,T),\\
                            & 						    &&\underline{x} \leq \Pi_{x}x(t) \leq \bar{x}, &&\quad t \in (0,\,T),\\
                            & 						    &&\underline{u} \leq \Pi_{u}u(t) \leq \bar{u}, &&\quad t \in (0,\,T),\\
                            & 						    &&\underline{z} \leq \Pi_{z}z(t) \leq \bar{z}, && \quad t \in (0,\,T),\\
                            & 						    &&\underline{c} \leq Cx(t) + Du(t) + Ez(t)\leq \bar{c}, &&\quad t \in (0,\,T), \\
                            &                           &&                                                   && \\[-1em]
                            & 						    &&F_T(x(T), z(T), p) = 0, &&\\
                            & 						    &&\underline{h}_T \leq g_T(h_T(x(T), z(T), p)) \leq \bar{h}_T, &&\\
                            & 						    &&\underline{x}_T \leq \Pi_{x_T}x(T) \leq \bar{u}_{T}, &&\\
                            & 						    &&\underline{z}_T \leq \Pi_{z_T}z(T) \leq \bar{z}_T, && \\
                            & 						    &&\underline{c}_T \leq C_Tx(T) + E_Tz(T)\leq \bar{c}_T, &&\\
\end{aligned}
\end{equation}
```

```eval_rst
Where:

* :math:`l: \mathbb{R}^{n_x}\times\mathbb{R}^{n_u}\times\mathbb{R}^{n_z} \rightarrow \mathbb{R}` is the Lagrange objective term.
* :math:`m: \mathbb{R}^{n_x}\times\mathbb{R}^{n_z} \rightarrow \mathbb{R}` is the Mayer objective term.

* :math:`F: \mathbb{R}^{n_x}\times\mathbb{R}^{n_x}\times\mathbb{R}^{n_u}\times\mathbb{R}^{n_z}\times\mathbb{R}^{n_p} \rightarrow \mathbb{R}^{n_x+n_z}` is the (potentially) fully implicit dynamics.
* :math:`F_T: \mathbb{R}^{n_x}\times\mathbb{R}^{n_z}\times\mathbb{R}^{n_p} \rightarrow \mathbb{R}^{n_x+n_z}` is the terminal algebraic constraint.

* :math:`h: \mathbb{R}^{n_x}\times\mathbb{R}^{n_u}\times\mathbb{R}^{n_z}\times\mathbb{R}^{n_p} \rightarrow \mathbb{R}^{n_h}` is the constraints general nonlinear function.
* :math:`h_T: \mathbb{R}^{n_x}\times\mathbb{R}^{n_z}\times\mathbb{R}^{n_p} \rightarrow \mathbb{R}^{n_{h_T}}` is the terminal constraints general nonlinear function.

* :math:`g: \mathbb{R}^{n_h} \rightarrow \mathbb{R}^{n_g}` is the constraints convex function.
* :math:`g_T: \mathbb{R}^{n_{h_T}} \rightarrow \mathbb{R}^{n_{g_T}}` is the terminal constraints convex function.

Currently not yet implemented features:

* :math:`l` must be in linear least-squares form :math:`l = \frac{1}{2}\| V_x x(t) + V_u u(t) + V_z z(t)\|_W^2`
* Support for soft constraints missing
* Constraints cannot depend on algebraic variables (yet)
```
