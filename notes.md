## Examples

### For quick developement

- `acados/examples/acados_python/pendulum_on_cart/solution_sensitivities/forw_vs_adj_param_sens.py`
  - `generate_solvers = False` for quick
  - exactly what we need

### Ultimate test
- `examples/acados_python/pendulum_on_cart/solution_sensitivities/policy_gradient_example.py`
  - vary p_global and compute with 2 solvers
  - for quick testing: set `build = False, generate = False`, and maybe use less values of p

- One solver `acados/examples/acados_python/solution_sensitivities_convex_example/batch_adjoint_solution_sensitivity_example.py`

## TODOs
- HPIPM hottest start
  - giaf: implement HPIPM
  - python: set this to 3 at creation?
  - prepare

- setter for `tau_min` (`mu0`)
  - Python: set at runtime TODO!

- just eval_qp_matrices_and_factorize
  - skelleton is there;
  - giaf: HPIPM and C

## Part 2:
- check if HPIPM solve is needed in sens solver.


# Done
- getter for `tau_iter`
  - giaf: done;