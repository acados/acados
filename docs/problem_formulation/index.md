# Problem Formulation

Since `acados` mainly aims on providing SQP type methods for optimal control, it naturally needs optimal control structured nonlinear programming formulations (OCP NLP) and quadratic programming (QP) formulations to tackle the subproblems within SQP.

## Optimal control structured NLP (OCP NLP)

The problem formulation targeted by `acados` OCP solver is stated here:
[https://github.com/acados/acados/blob/master/docs/problem_formulation/problem_formulation_ocp_mex.pdf](https://github.com/acados/acados/blob/master/docs/problem_formulation/problem_formulation_ocp_mex.pdf)


## QP formulations (dense and OCP structured)

`acados` relies on `HPIPM` for reformulating QP problems via (partial) condensing and expansion routines.
In order to use them efficiently, we use the flexible QP formulations from `HPIPM`.
Those are the optimal control structured quadratic programming formulation (OCP QP) and the
dense QP formulation.

Both problem formulations are documented in the [`HPIPM guide`](https://github.com/giaf/hpipm/blob/master/doc/guide.pdf).

