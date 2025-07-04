# Features by Example

This page showcases how specific `acados` features can be used by pointing to the relevant examples.
If you are new to `acados`, we highly recommend you to start with the `getting_started` examples:
- [Python](https://github.com/acados/acados/blob/main/examples/acados_python/getting_started)
- [MATLAB/Octave](https://github.com/acados/acados/blob/main/examples/acados_matlab_octave/getting_started)

### Getting started with Model Predictive Control
In particular, the **closed-loop** examples are a great starting point to develop model predictive control (MPC) in a simulation.
Here, an OCP solver (`AcadosOcpSolver`) is used to compute the control inputs and an integrator (`AcadosSimSolver`) is used to simulate the real system.
- [Python](https://github.com/acados/acados/blob/main/examples/acados_python/getting_started/minimal_example_closed_loop.py)
- [MATLAB/Octave](https://github.com/acados/acados/blob/main/examples/acados_matlab_octave/getting_started/minimal_example_closed_loop.m)

## Simulation and Sensitivity propagation
- [Python](https://github.com/acados/acados/blob/main/examples/acados_python/pendulum_on_cart/sim/extensive_example_sim.py)
- [MATLAB/Octave and Simulink](https://github.com/acados/acados/blob/main/examples/acados_matlab_octave/getting_started/minimal_example_sim.m)


## Optimal Control Problem Formulation

### Parameter updates
- [Python](https://github.com/acados/acados/blob/main/examples/acados_python/tests/test_parametric_nonlinear_constraint_h.py)
- [MATLAB/Octave and Simulink](https://github.com/acados/acados/blob/main/examples/acados_matlab_octave/test/param_test.m)


### Cost formulations

`acados` supports general nonlinear cost, but can also exploit particular cost structures such as (non)linear least-squares costs and convex-over-nonlinear costs.

- [Python](https://github.com/acados/acados/blob/main/examples/acados_python/pendulum_on_cart/ocp/ocp_example_cost_formulations.py)


### Soft constraints
- [Python](https://github.com/acados/acados/blob/main/examples/acados_python/tests/soft_constraint_test.py)
- [MATLAB/Octave](https://github.com/acados/acados/blob/main/examples/acados_matlab_octave/test/create_slacked_ocp_qp_solver_formulation.m)


### Convex-over-nonlinear constraints
- [Python](https://github.com/acados/acados/blob/main/examples/acados_python/convex_ocp_with_onesided_constraints/main_convex_onesided.py)
- currently not supported for MATLAB

### Multi-phase OCP
- [Python](https://github.com/acados/acados/blob/main/examples/acados_python/mocp_transition_example)
- [MATLAB/Octave and Simulink](https://github.com/acados/acados/blob/main/examples/acados_matlab_octave/mocp_transition_example)

### Constraints and cost on control rate
- [MATLAB/Octave and Simulink](https://github.com/acados/acados/blob/main/examples/acados_matlab_octave/control_rates/main.m)


### Moving horizon estimation (MHE)
- [Python](https://github.com/acados/acados/blob/main/examples/acados_python/pendulum_on_cart/mhe/closed_loop_mhe_ocp.py)
- [MATLAB/Octave](https://github.com/acados/acados/blob/main/examples/acados_matlab_octave/lorentz/example_mhe.m)

### Discretization with a nonuniform grid
- [Python](https://github.com/acados/acados/blob/main/examples/acados_python/furuta_pendulum/main_closed_loop.py)
- [MATLAB/Octave and Simulink](https://github.com/acados/acados/blob/main/examples/acados_matlab_octave/getting_started/extensive_example_ocp.m)

### Time-varying reference tracking
- [Python](https://github.com/acados/acados/blob/main/examples/acados_python/pendulum_on_cart/mhe/closed_loop_mhe_ocp.py)
- [MATLAB/Octave](https://github.com/acados/acados/blob/main/examples/acados_matlab_octave/control_rates/main.m)


### One-sided constraints
- [Python](https://github.com/acados/acados/blob/main/examples/acados_python/convex_ocp_with_onesided_constraints/main_convex_onesided.py)
- [MATLAB/Octave](https://github.com/acados/acados/blob/main/examples/acados_matlab_octave/linear_mpc/main.m)



## Algorithmic features and solver options

### Real-time iterations (RTI)
- [Python](https://github.com/acados/acados/blob/main/examples/acados_python/furuta_pendulum/main_closed_loop.py)
- [MATLAB/Octave](https://github.com/acados/acados/blob/main/examples/acados_matlab_octave/masses_chain_model/example_closed_loop.m)

### Advanced-step real-time iterations (AS-RTI)
Relevant publications: [Frey2024a](https://publications.syscop.de/Frey2024a.pdf), [Nurkanovic2019a](https://publications.syscop.de/Nurkanovic2019a.pdf)
- [Python](https://github.com/acados/acados/blob/main/examples/acados_python/pendulum_on_cart/as_rti/as_rti_closed_loop_example.py)


### Solver timeout
- [Python](https://github.com/acados/acados/blob/main/examples/acados_python/furuta_pendulum/main_closed_loop.py)
- [MATLAB/Octave and Simulink](https://github.com/acados/acados/blob/main/examples/acados_matlab_octave/control_rates/main.m)


### Cost integration
- [Python](https://github.com/acados/acados/blob/main/examples/acados_python/tests/test_ggn_cost_integration.py)

### Globalization
- [Python](https://github.com/acados/acados/blob/main/examples/acados_python/convex_problem_globalization_needed/convex_problem_globalization_necessary.py)

### Differential dynamic programming (DDP)
- [Python](https://github.com/acados/acados/blob/main/examples/acados_python/unconstrained_ocps/hour_glass_p2p_motion/hour_glass_time_optimal_p2p_motion.py)

### Partial condensing

Use a `qp_solver` starting with `PARTIAL_CONDENSING`, use `qp_solver_cond_N` to set the horizon of the partially condensed QP.
Additionally, one can use `qp_solver_cond_block_size` to specify how many blocks are condensed into one block in partial condensing.
- [Python](https://github.com/acados/acados/blob/main/examples/acados_python/pendulum_on_cart/ocp/nonuniform_discretization_example.py)
- [MATLAB/Octave](https://github.com/acados/acados/blob/main/examples/acados_matlab_octave/p_global_example/set_solver_options.m)


### Zero-order robust optimization (zoRO)
Relevant publications: [Frey2024](https://publications.syscop.de/Frey2024.pdf), [Zanelli2021](https://publications.syscop.de/Zanelli2021.pdf)

- [Python](https://github.com/acados/acados/blob/main/examples/acados_python/zoRO_example)
- [MATLAB/Octave](https://github.com/acados/acados/blob/main/examples/acados_matlab_octave/pendulum_on_cart_model/zoro_example.m)


### Solution sensitivities
- [Python](https://github.com/acados/acados/blob/main/examples/acados_python/chain_mass/solution_sensitivity_example.py)
- currently not supported for MATLAB/Octave

### Adjoint solution sensitivities
- [Python](https://github.com/acados/acados/blob/main/examples/acados_python/pendulum_on_cart/solution_sensitivities/forw_vs_adj_param_sens.py)
- currently not supported for MATLAB/Octave

