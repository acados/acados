## Problem Formulation

Since `acados` mainly aims on providing SQP type methods for optimal control, it naturally needs optimal control structured nonlinear programming formulations (OCP-NLP) and quadratic programming (QP) formulations to tackle the subproblems within SQP.

- __Optimal control structured NLP (OCP-NLP)__: The problem formulation targeted by `acados` OCP solver is stated [here](https://github.com/acados/acados/blob/main/docs/problem_formulation/problem_formulation_ocp_mex.pdf).

- __QP formulations (dense and OCP structured)__: `acados` relies on `HPIPM` for reformulating QP problems via (partial) condensing and expansion routines.
We thus use the flexible QP formulations from `HPIPM` for optimal control structured quadratic programming formulation (OCP-QP) and the dense QP formulation.
Both problem formulations are documented in the [`HPIPM guide`](https://github.com/giaf/hpipm/blob/master/doc/guide.pdf).

<!-- TODO: AcadosSim?! -->