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
- check if HPIPM solve is needed in sens solver.
  - last
- getter for `mu_target`
  - giaf
- HPIPM hottest start
  - giaf
- setter for `tau_min` (`mu0`)
  - setters are there;
  - OJ: TODO: set to same value;

### Optional
- just factorize OR for 1 adjoint, already use the correct rhs
  - OJ: prepare down to C.