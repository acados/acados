
# Troubleshooting

As a first step, check the [solver status](https://docs.acados.org/python_interface/index.html#acados_template.acados_ocp_solver.AcadosOcpSolver.get_status) which is returned by `solve()` in Python and can be obtained with `solver.get_status()` and `solver.get('status')` in Python and MATLAB/Octave, respectively.

Next, check the NLP residuals by calling `solver.print_statistics()` in Python and `solver.print('stat')` in MATLAB/Octave.

In order to asses the QP solver status, you need to check the corresponding QP solver status definitions. For `HPIPM`, they are given [here](https://github.com/giaf/hpipm/blob/deb7808e49a3cc2b1bdb721cba23f13869c0a35c/include/hpipm_common.h#L57).

### Solver initialization
Use the `set` method of `AcadosOcpSolver` to provide different initializations.
Are all problem functions defined at the initial guess?
`NaN` errors can often be mitigated by solver initializations.

### QP diagnostics
- [Python](https://github.com/acados/acados/blob/main/examples/acados_python/pendulum_on_cart/solution_sensitivities/policy_gradient_example.py)
- [MATLAB/Octave](https://github.com/acados/acados/blob/main/examples/acados_matlab_octave/getting_started/extensive_example_ocp.m)


### Store iterates
One can always get the last iterate using `solver.get()`, see
- [Python](https://github.com/acados/acados/blob/main/examples/acados_python/linear_mass_model/linear_mass_test_problem.py)
- [MATLAB/Octave](https://github.com/acados/acados/blob/main/examples/acados_matlab_octave/getting_started/extensive_example_ocp.m)

In addition, one can set the solver option [`store_iterates`](https://docs.acados.org/python_interface/index.html#acados_template.acados_ocp_options.AcadosOcpOptions.store_iterates) to store all intermediate NLP solver iterates and get them after a solver call.
- [Python](https://github.com/acados/acados/blob/main/examples/acados_python/convex_ocp_with_onesided_constraints/main_convex_onesided.py)
- [MATLAB/Octave](https://github.com/acados/acados/blob/main/examples/acados_matlab_octave/getting_started/extensive_example_ocp.m)

