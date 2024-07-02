# Developer Guide
This page contains additional information for people who want to extend `acados`.

## Contributing
Contributions are handled via pull requests on Github
- Fork the project
- Push your changes to your fork
- Open a pull request https://github.com/acados/acados
- Describe what you changed and why
- Rather make minimal changes
- Rather more commits than less

## Memory management in `acados`
The following are guidelines on how memory should be assigned for an `acados` structure `astruct`.

There are two functions: `astruct_calculate_size()`, `astruct_assign()`.

### `astruct_calculate_size()`
Must return a multiple of 8 to keep the pointer aligned to 8 bytes when allocating substructures.
Thus, it should end with:
```
    make_int_multiple_of(8, &size);
```

### `astruct_assign()`
Should assign its members in the following order:

- Align to 8 bytes, i.e.:
```
    align_char_to(8, &c_ptr);
```

- Assign structure itself, i.e.:
```
    astruct *as_instance = (astruct *) c_ptr;
    c_ptr += sizeof(astruct);
```

- Assign pointers to substructures, if `astruct` contains an array of `bstruct`s, e.g.
```
    mem->bstructs = (void **) c_ptr;
    c_ptr += N*sizeof(void *);
```

- Align to 8 bytes, since `astruct` might contain `int`s and the pointers were assigned.

- Assign "substructures", i.e. structures that `astruct` has pointers to:
```
    // assume astruct contains an instance of bstruct
    as_instance->bstruct = bstruct_assign(c_ptr,...);
    c_ptr += bstruct_calculate_size(astruct);
```

- Note:
    - since calculate_size returns multiple of 8, `c_ptr` is still aligned to 8 bytes.
    - `blasfeo_dmat_structs`, `blasfeo_dvec_structs` can be seen as "substructures" here.


- Assign doubles (are 8 bytes anyway)
```
    assign_and_advance_double(n_doubles, &as_instance->doubles, &c_ptr);
```

- Assign pointers (4 bytes on 32 Bit)
```
    assign_and_advance_double_ptrs(n_pointers, &as_instance->double_pointer, &c_ptr);
```

- Align to 8 bytes, can be skipped if no pointers have been assigned

- Assign integers
```
    assign_and_advance_int(n_integers, &as_instance->ints, &c_ptr);
```

- Align to 64 bytes, i.e.:
```
    align_char_to(64, &c_ptr);
```

- Assign `blasfeo_dmat_mem` (are multiple of 64 Bytes)
```
    assign_and_advance_blasfeo_dmat_mem(nrow, ncol, &as_instance->blasfeo_mat, &c_ptr);
```

- Assign `blasfeo_vec_mem` (are multiple of 32 Bytes)
```
    assign_and_advance_blasfeo_dvec_mem(n, &as_instance->blasfeo_vec, &c_ptr);
```

- Align c_ptr to 8 byte here for nested assigns, see "substructures"
   - relevant if no `blasfeo_mems` are in `astruct`


## Regularization within SQP / SQP-RTI
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



## Dense QP solution: Populating `dense_qp_out`
After solving a dense QP, the solution should be stored in the `dense_qp_out` structure that is passed as an argument to the function `dense_qp_XXX` (where `XXX` is the name of the solver).
This structure has to be populated with three things:
- the primal solution,
- the dual solution,
- constraint slacks,

which corresponds to the fields `v`, `lam`, and `t`, respectively, in the `dense_qp_out` structure.
If there are, in addition, equality constraint their dual variables should be added to the field `pi`.

### Primal and soft slack variables
The field `v` should be populated with variables in the following order:
```
primal, lower soft slack, upper soft slack
```
with the corresponding dimensions `nv`,`ns`,`ns`.

### Dual variables
The field `lam` should be populated with dual variables corresponding to bounds in the following order:
```
lower bounds, lower general, upper bounds, upper general, lower soft slack, upper soft slack
```
with the corresponding dimensions `nb`,`ng`,`nb`,`ng`,`ns`,`ns`.

The sign convention used is that all dual variables are nonnegative, even the ones that correspond to lower bounds.

### Constraint slacks
If `v` has been set correctly, the constraint slack `t` can be computed using the auxiliary function `dense_qp_compute_t()`.

### Adding tests
In addition to specifying new tests using `add_test()` to the `interface/CMakeLists.txt`, please ensure that `set_tests_properties()` is also used to force serial execution for tests contained in the same directory.
