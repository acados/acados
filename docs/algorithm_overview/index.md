# Algorithm overview

This page is work in progress..

<!-- ## SQP / SQP-RTI -->

## Regularization

The Hessian of the QP is computed in the `ocp_nlp_sqp`, `ocp_nlp_sqp_rti` module respectively.

The following steps are carried out:

- `ocp_nlp_approximate_qp_matrices()`
  - sets all Hessian blocks to 0.0
  - for `i in range(N+1)`
    - if `i<N:`
      - add to the diagonal of the Hessian of block `i` the term `in->Ts[i] * opts->levenberg_marquardt`
      - add the contribution of the dynamics module (can be turned off via `exact_hess_dyn`)
    - if `i==N:`
      -  add to the diagonal of the Hessian of block `N` the term `opts->levenberg_marquardt`.
  - add the cost contribution to the Hessian
    - Gauss-Newton Hessian (available in least-squares case)
    - or exact Hessian (always used with `EXTERNAL` cost module) if no "custom hessian" is set (see `cost_expr_ext_cost_custom_hess`, `cost_expr_ext_cost_custom_hess_0`, `cost_expr_ext_cost_custom_hess_e`)
  - add the inequality constraints contribution to the Hessian (can be turned off via `exact_hess_constr`)

- call the regularization module (`regularize`, see [`regularize_method`](https://docs.acados.org/python_interface/index.html?highlight=regularize#acados_template.acados_ocp_options.AcadosOcpOptions.regularize_method))

<!-- TODO: change this to have a seperate levenberg_marquardt term on the terminal stage (instead of 1 replacing Ts).
+ add the option to provide a vector that is added on diagonal, i.e. make levenberg_marquardt a vector of size nx+nu. -->