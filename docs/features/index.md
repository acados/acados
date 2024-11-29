# Features by Example
This page showcases how specific `acados` features can be used, by pointing to the relevant examples.
If you are new to `acados` we highly recommend you to start with the `getting_started` examples:
- for Python: [`examples/acados_python/getting_started`](https://github.com/acados/acados/blob/master/examples/acados_python/getting_started)
- for MATLAB, Octave: [`examples/acados_matlab_octave/getting_started`](https://github.com/acados/acados/blob/master/examples/acados_matlab_octave/getting_started)


## Simulation and Sensitivity propagation via `AcadosSimSolver`
- for Python: [`examples/acados_python/pendulum_on_cart/sim/extensive_example_sim.py`](https://github.com/acados/acados/blob/master/examples/acados_python/pendulum_on_cart/sim/extensive_example_sim.py)
- for MATLAB, Octave, Simulink: [`examples/acados_matlab_octave/getting_started/minimal_example_sim.m`](https://github.com/acados/acados/blob/master/examples/acados_matlab_octave/getting_started/minimal_example_sim.m)


## Optimal Control Problem Formulation
### Parameter updates
- for Python: [`examples/acados_python/tests/test_parametric_nonlinear_constraint_h.py`](https://github.com/acados/acados/blob/master/examples/acados_python/tests/test_parametric_nonlinear_constraint_h.py)
- for MATLAB, Octave, Simulink: [`examples/acados_matlab_octave/test/param_test.m`](https://github.com/acados/acados/blob/master/examples/acados_matlab_octave/test/param_test.m)

### Solver time-out
- for Python: [`examples/acados_python/furuta_pendulum/main_closed_loop.py`](https://github.com/acados/acados/blob/master/examples/acados_python/furuta_pendulum/main_closed_loop.py)
- for MATLAB, Octave, Simulink: [`examples/acados_matlab_octave/control_rates/main.m`](https://github.com/acados/acados/blob/master/examples/acados_matlab_octave/control_rates/main.m)

### Convex-over-nonlinear constraints
- for Python: [`examples/acados_python/convex_ocp_with_onesided_constraints/main_convex_onesided.py`](https://github.com/acados/acados/blob/master/examples/acados_python/convex_ocp_with_onesided_constraints/main_convex_onesided.py)

### Cost formulations
- for Python: [`examples/acados_python/pendulum_on_cart/ocp/ocp_example_cost_formulations.py`](examples/acados_python/pendulum_on_cart/ocp/ocp_example_cost_formulations.py)

### Multi-phase OCP
- for Python: [`examples/acados_python/mocp_transition_example`](https://github.com/acados/acados/blob/master/examples/acados_python/mocp_transition_example)
- for MATLAB, Octave, Simulink: [`examples/acados_matlab_octave/mocp_transition_example`](https://github.com/acados/acados/blob/master/examples/acados_matlab_octave/mocp_transition_example)

### Rate constraints & cost
- for MATLAB, Octave, Simulink: [`examples/acados_matlab_octave/control_rates/main.m`](https://github.com/acados/acados/blob/master/examples/acados_matlab_octave/control_rates/main.m)

### Soft constraints
- for Python: `examples/acados_python/linear_mass_model/solve_marathos_ocp.py`(https://github.com/acados/acados/blob/master/examples/acados_python/linear_mass_model/solve_marathos_ocp.py)
- for MATLAB, Octave: [`examples/acados_matlab_octave/test/create_slacked_ocp_qp_solver_formulation.m`](https://github.com/acados/acados/blob/master/examples/acados_matlab_octave/test/create_slacked_ocp_qp_solver_formulation.m)


### Moving horizon estimation (MHE)
- for Python: `examples/acados_python/pendulum_on_cart/mhe/closed_loop_mhe_ocp.py`(https://github.com/acados/acados/blob/master/examples/acados_python/pendulum_on_cart/mhe/closed_loop_mhe_ocp.py)
- for MATLAB, Octave: [`examples/acados_matlab_octave/lorentz/example_mhe.m`](https://github.com/acados/acados/blob/master/examples/acados_matlab_octave/lorentz/example_mhe.m)

### Nonuniform grid OCP discretization
- [`examples/acados_python/furuta_pendulum/main_closed_loop.py`](https://github.com/acados/acados/blob/master/examples/acados_python/furuta_pendulum/main_closed_loop.py)
- [`examples/acados_matlab_octave/getting_started/extensive_example_ocp.m`](https://github.com/acados/acados/blob/master/examples/acados_matlab_octave/getting_started/extensive_example_ocp.m)

### Time varying reference tracking
- for MATLAB, Octave: [`examples/acados_matlab_octave/control_rates/main.m`](https://github.com/acados/acados/blob/master/examples/acados_matlab_octave/control_rates/main.m)
- for Python: `examples/acados_python/pendulum_on_cart/mhe/closed_loop_mhe_ocp.py`(https://github.com/acados/acados/blob/master/examples/acados_python/pendulum_on_cart/mhe/closed_loop_mhe_ocp.py)


### One-sided constraints
- for Python: [`examples/acados_python/convex_ocp_with_onesided_constraints/main_convex_onesided.py`](https://github.com/acados/acados/blob/master/examples/acados_python/convex_ocp_with_onesided_constraints/main_convex_onesided.py)
- for MATLAB, Octave: [`examples/acados_matlab_octave/linear_mpc/main.m`](https://github.com/acados/acados/blob/master/examples/acados_matlab_octave/linear_mpc/main.m)



## Algorithmic features and solver options

### Real-time Iterations (RTI)
- for Python: [`examples/acados_python/furuta_pendulum/main_closed_loop.py`](https://github.com/acados/acados/blob/master/examples/acados_python/furuta_pendulum/main_closed_loop.py)
- for MATLAB, Octave: [`examples/acados_matlab_octave/masses_chain_model/example_closed_loop.m`](https://github.com/acados/acados/blob/master/examples/acados_matlab_octave/masses_chain_model/example_closed_loop.m)

### Advanced-Step Real-Time Iterations (AS-RTI)
Relevant publications: [Frey2024a](https://publications.syscop.de/Frey2024a.pdf), [Nurkanovic2019a](https://publications.syscop.de/Nurkanovic2019a.pdf)
- for Python: [`examples/acados_python/pendulum_on_cart/as_rti`](https://github.com/acados/acados/blob/master/examples/acados_python/pendulum_on_cart/as_rti)


### Cost integration
- for Python: [`examples/acados_python/tests/test_ggn_cost_integration.py`](examples/acados_python/tests/test_ggn_cost_integration.py)

### Globalization
- for Python: [`examples/acados_python/convex_problem_globalization_needed/convex_problem_globalization_necessary.py`](examples/acados_python/convex_problem_globalization_needed/convex_problem_globalization_necessary.py)

### Differential Dynamic Programming (DDP)
- for Python: [`examples/acados_python/unconstrained_ocps/hour_glass_p2p_motion/hour_glass_time_optimal_p2p_motion.py`](examples/acados_python/unconstrained_ocps/hour_glass_p2p_motion/hour_glass_time_optimal_p2p_motion.py)

### Partial condensing
Use a `qp_solver` starting with `PARTIAL_CONDENSING`, use `qp_solver_cond_N` to set the horizon of the partially condensed QP.
Additionally, one can use `qp_solver_cond_block_size` to specify how many blocks are condensed into one block in partial condensing.
- for Python: [`examples/acados_python/pendulum_on_cart/ocp/nonuniform_discretization_example.py`](examples/acados_python/pendulum_on_cart/ocp/nonuniform_discretization_example.py)
- for MATLAB, Octave, Simulink: [`examples/acados_matlab_octave/p_global_example/set_solver_options.m`](https://github.com/acados/acados/blob/master/examples/acados_matlab_octave/p_global_example/set_solver_options.m)


### zero-order robust optimization (zoRO)
Relevant publications:
[Frey2024](https://publications.syscop.de/Frey2024.pdf)
[Zanelli2021](https://publications.syscop.de/Zanelli2021.pdf)

Examples
- for Python: [`examples/acados_python/zoRO_example`](https://github.com/acados/acados/blob/master/examples/acados_python/zoRO_example)
- for MATLAB, Octave, Simulink: [`examples/acados_matlab_octave/pendulum_on_cart_model/zoro_example.m`](https://github.com/acados/acados/blob/master/examples/acados_matlab_octave/pendulum_on_cart_model/zoro_example.m)


### Solution sensitivities
- for Python: [`examples/acados_python/chain_mass/solution_sensitivity_example.py`](examples/acados_python/chain_mass/solution_sensitivity_example.py)

### Adjoiont solution sensitivities
- for Python: [`examples/acados_python/pendulum_on_cart/solution_sensitivities/forw_vs_adj_param_sens.py`](examples/acados_python/pendulum_on_cart/solution_sensitivities/forw_vs_adj_param_sens.py)




## Trouble shooting
Check Solver status: see [documentation of acados status for NLP solver status.](https://docs.acados.org/python_interface/index.html#acados_template.acados_ocp_solver.AcadosOcpSolver.get_status)
In order to asses the QP solver status, one has to check the corresponding QP solver status definitions.


### QP diagnostics
- for Python: [`examples/acados_python/pendulum_on_cart/solution_sensitivities/policy_gradient_example.py`](examples/acados_python/pendulum_on_cart/solution_sensitivities/policy_gradient_example.py)
- for MATLAB, Octave: [`examples/acados_matlab_octave/getting_started/extensive_example_ocp.m`](https://github.com/acados/acados/blob/master/examples/acados_matlab_octave/getting_started/extensive_example_ocp.m)


### Store iterates
One can always get the last iterate using `ocp_solver.get()`, see
- for Python: [`examples/acados_python/linear_mass_model/linear_mass_test_problem.py`](examples/acados_python/linear_mass_model/linear_mass_test_problem.py)
- for MATLAB, Octave: [`examples/acados_matlab_octave/getting_started/extensive_example_ocp.m`](https://github.com/acados/acados/blob/master/examples/acados_matlab_octave/getting_started/extensive_example_ocp.m)

In addition, one can set the solver option [`store_iterates`](https://docs.acados.org/python_interface/index.html#acados_template.acados_ocp_options.AcadosOcpOptions.store_iterates) to store all intermediate NLP solver iterates and get them after a solver call.
- for Python: [`examples/acados_python/convex_ocp_with_onesided_constraints/main_convex_onesided.py`](examples/acados_python/convex_ocp_with_onesided_constraints/main_convex_onesided.py)
- for MATLAB, Octave: [`examples/acados_matlab_octave/getting_started/extensive_example_ocp.m`](https://github.com/acados/acados/blob/master/examples/acados_matlab_octave/getting_started/extensive_example_ocp.m)

