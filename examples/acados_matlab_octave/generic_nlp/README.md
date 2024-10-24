### Solving a generic NLP with `acados`

This example solves a generic nonlinear program of the form:
```math
\min_x\quad f(x, p)
```
```math
\hspace{2.5cm}\text{s.t.}\quad g_{lb} \leq g(x, p) \leq g_{ub}
```
where $p$ is a vector of parameters. This can be done with the following:
- prediction horizon is equal to 1
- the initial state is unconstrained
- set $\dot{x} = 0$
- use $f(x,p)$ as a terminal cost term
- use a cheap integrator ERK with 1 stage and 1 step
- set $g(x,p)$ as a terminal constraint

Note that you might need a more careful initialization compared to, e.g., CasADi + IPOPT. On the other hand, there could be a significant speed up, depending on your problem.

The NLP solved in this example is:
```math
\hspace{2.3cm}\min_{x}\quad p_1(100(x_2-x_1^2)^2+(x_1-1)^2)
```
```math
\text{s.t.}\quad x_1^2+x_2^2-2p_2=0
```
More details and the original code can be found in [this forum post](https://discourse.acados.org/t/solving-simple-nlp-problem-exploit-blasfeo-performance/271).