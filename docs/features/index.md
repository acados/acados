# Features
This page showcases how specific `acados` features can be used, by pointing to the relevant examples.
If you are new to `acados` we highly recommend you to start with the `getting_started` examples:
- for Python: [`examples/acados_python/getting_started`](https://github.com/acados/acados/blob/master/examples/acados_python/getting_started)
- for MATLAB, Octave, Simulink: [`examples/acados_matlab_octave/getting_started`](https://github.com/acados/acados/blob/master/examples/acados_matlab_octave/getting_started)


## Simulation via AcadosSimSolver & Sensitivity propagation


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

### Moving horizon estimation (MHE)
- for Python: `examples/acados_python/pendulum_on_cart/mhe/closed_loop_mhe_ocp.py`(https://github.com/acados/acados/blob/master/examples/acados_python/pendulum_on_cart/mhe/closed_loop_mhe_ocp.py)
- for MATLAB / Octave: [`examples/acados_matlab_octave/lorentz/example_mhe.m`](https://github.com/acados/acados/blob/master/examples/acados_matlab_octave/lorentz/example_mhe.m)

### Nonuniform grid OCP discretization
- [`examples/acados_python/furuta_pendulum/main_closed_loop.py`](https://github.com/acados/acados/blob/master/examples/acados_python/furuta_pendulum/main_closed_loop.py)
- [`examples/acados_matlab_octave/getting_started/extensive_example_ocp.m`](https://github.com/acados/acados/blob/master/examples/acados_matlab_octave/getting_started/extensive_example_ocp.m)

### Time varying reference tracking

### One-sided constraints
- for Python: [`examples/acados_python/convex_ocp_with_onesided_constraints/main_convex_onesided.py`](https://github.com/acados/acados/blob/master/examples/acados_python/convex_ocp_with_onesided_constraints/main_convex_onesided.py)
- for MATLAB / Octave: [`examples/acados_matlab_octave/linear_mpc/main.m`](https://github.com/acados/acados/blob/master/examples/acados_matlab_octave/linear_mpc/main.m)


## Algorithmic features and solver options

### Real-time Iterations (RTI)
- for Python: [`examples/acados_python/furuta_pendulum/main_closed_loop.py`](https://github.com/acados/acados/blob/master/examples/acados_python/furuta_pendulum/main_closed_loop.py)
- for MATLAB / Octave: [`examples/acados_matlab_octave/masses_chain_model/example_closed_loop.m`](https://github.com/acados/acados/blob/master/examples/acados_matlab_octave/masses_chain_model/example_closed_loop.m)

### Advanced-Step Real-Time Iterations (AS-RTI)
Relevant publications: [Frey2024a](https://publications.syscop.de/Frey2024a.pdf), [Nurkanovic2019a](https://publications.syscop.de/Nurkanovic2019a.pdf)
- for Python: [`examples/acados_python/pendulum_on_cart/as_rti`](https://github.com/acados/acados/blob/master/examples/acados_python/pendulum_on_cart/as_rti)


### Cost integration
- for Python: [`examples/acados_python/tests/test_ggn_cost_integration.py`](examples/acados_python/tests/test_ggn_cost_integration.py)

### Solution sensitivities

### Globalization
- for Python: [`examples/acados_python/convex_problem_globalization_needed/convex_problem_globalization_necessary.py`](examples/acados_python/convex_problem_globalization_needed/convex_problem_globalization_necessary.py)

### DDP

### Partial condensing

### zero-order robust optimization (zoRO)
Relevant publications:
[Frey2024](https://publications.syscop.de/Frey2024.pdf)
[Zanelli2021](https://publications.syscop.de/Zanelli2021.pdf)

Examples
- for Python: [`examples/acados_python/zoRO_example`](https://github.com/acados/acados/blob/master/examples/acados_python/zoRO_example)
- for MATLAB, Octave, Simulink: [`examples/acados_matlab_octave/pendulum_on_cart_model/zoro_example.m`](https://github.com/acados/acados/blob/master/examples/acados_matlab_octave/pendulum_on_cart_model/zoro_example.m)

## Trouble shooting
Solver status.

### Store iterates

### QP diagnostics





