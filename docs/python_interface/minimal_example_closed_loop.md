# Minimal Closed-Loop MPC Example - Walkthrough

This walkthrough provides a detailed explanation of the [minimal_example_closed_loop.py](https://github.com/acados/acados/blob/main/examples/acados_python/getting_started/minimal_example_closed_loop.py) example, which demonstrates how to implement a closed-loop Model Predictive Control (MPC) simulation using `acados` in Python.

## Overview

The example simulates a **pendulum on a cart** system using MPC. In closed-loop control:
1. An OCP solver computes optimal control inputs based on the current state
2. These control inputs are applied to the system
3. An integrator simulates the system's response to get the next state
4. The process repeats for multiple time steps

This creates a feedback loop that mimics real-world MPC applications.

## Problem Setup

The system is an inverted pendulum on a cart with:
- **States** (4): cart position `x`, pendulum angle `θ`, cart velocity `v`, angular velocity `dθ/dt`
- **Control** (1): horizontal force `F` applied to the cart
- **Goal**: Swing up the pendulum from hanging down (θ = π) to upright (θ = 0) and stabilize it

### System Parameters
```python
x0 = np.array([0.0, np.pi, 0.0, 0.0])  # Initial state: cart at origin, pendulum hanging down
Fmax = 80                               # Maximum control force [N]
Tf = 0.8                                # Prediction horizon time [s]
N_horizon = 40                          # Number of shooting intervals
Nsim = 100                              # Number of simulation steps
```

## Step-by-Step Walkthrough

### 1. Define the Dynamic Model

The file `pendulum_model.py` exports the pendulum dynamics using CasADi symbolic expressions:

```python
def export_pendulum_ode_model() -> AcadosModel:
    # System parameters
    M = 1.0    # mass of cart [kg]
    m = 0.1    # mass of pendulum [kg]
    g = 9.81   # gravity [m/s^2]
    l = 0.8    # length of pendulum rod [m]
    
    # States and controls
    x1, theta, v1, dtheta = SX.sym('x1'), SX.sym('theta'), SX.sym('v1'), SX.sym('dtheta')
    x = vertcat(x1, theta, v1, dtheta)
    F = SX.sym('F')
    u = vertcat(F)
    
    # Dynamics equations (derived from physics)
    f_expl = vertcat(
        v1,
        dtheta,
        (-m*l*sin(theta)*dtheta**2 + m*g*cos(theta)*sin(theta) + F) / denominator,
        (-m*l*cos(theta)*sin(theta)*dtheta**2 + F*cos(theta) + (M+m)*g*sin(theta)) / (l*denominator)
    )
    
    model = AcadosModel()
    model.f_expl_expr = f_expl
    model.x = x
    model.u = u
    return model
```

### 2. Setup the OCP Solver

The function `setup_ocp_solver()` configures the optimal control problem:

```python
def setup_ocp_solver(x0, Fmax, N_horizon, Tf):
    ocp = AcadosOcp()
    
    # Set the model
    model = export_pendulum_ode_model()
    ocp.model = model
```

#### Cost Function
The cost penalizes deviations from the target state (upright pendulum at origin) and control effort:

```python
Q_mat = 2*np.diag([1e3, 1e3, 1e-2, 1e-2])  # State weights
R_mat = 2*np.diag([1e-2])                   # Control weight

ocp.cost.cost_type = 'NONLINEAR_LS'         # Nonlinear least squares
ocp.cost.W = scipy.linalg.block_diag(Q_mat, R_mat)  # Stage cost weights
ocp.cost.W_e = Q_mat                         # Terminal cost weights

ocp.model.cost_y_expr = vertcat(model.x, model.u)    # Stage cost terms
ocp.model.cost_y_expr_e = model.x                     # Terminal cost terms
ocp.cost.yref = np.zeros((nx + nu,))         # Target is zero (upright, at rest)
ocp.cost.yref_e = np.zeros((nx,))
```

The cost function minimizes: `|| x - 0 ||²_Q + || u - 0 ||²_R` at each stage, plus `|| x - 0 ||²_Q` at the terminal stage.

#### Constraints
Control input constraints limit the force:

```python
ocp.constraints.lbu = np.array([-Fmax])  # Lower bound
ocp.constraints.ubu = np.array([+Fmax])  # Upper bound
ocp.constraints.idxbu = np.array([0])    # Apply to control input 0
ocp.constraints.x0 = x0                  # Initial state constraint
```

#### Solver Options
```python
ocp.solver_options.N_horizon = N_horizon   # Discretization
ocp.solver_options.tf = Tf                 # Prediction horizon
ocp.solver_options.hessian_approx = 'GAUSS_NEWTON'  # Fast approximation
ocp.solver_options.qp_tol = 1e-8          # Solver tolerance
```

Finally, create the solver:
```python
ocp.code_export_directory = 'c_generated_code_ocp'
acados_ocp_solver = AcadosOcpSolver(ocp)
return acados_ocp_solver
```

This generates C code and compiles a fast solver specific to this problem.

### 3. Setup the Integrator

The function `setup_integrator()` creates a simulator for the "real" system:

```python
def setup_integrator(dt):
    sim = AcadosSim()
    sim.model = export_pendulum_ode_model()  # Same dynamics as OCP
    
    sim.solver_options.T = dt           # Integration time step
    sim.solver_options.num_steps = 2    # Integration sub-steps (for accuracy)
    
    sim.code_export_directory = 'c_generated_code_sim'
    acados_integrator = AcadosSimSolver(sim)
    return acados_integrator
```

The integrator is more accurate than the OCP's internal integrator, simulating model mismatch.

### 4. Closed-Loop Simulation

The `main()` function runs the closed-loop simulation:

```python
def main():
    # Initialize solvers
    ocp_solver = setup_ocp_solver(x0, Fmax, N_horizon, Tf)
    integrator = setup_integrator(Tf/N_horizon)
    
    # Dimensions
    nx = ocp_solver.acados_ocp.dims.nx  # Number of states
    nu = ocp_solver.acados_ocp.dims.nu  # Number of controls
    
    # Allocate storage
    simX = np.zeros((Nsim+1, nx))  # State trajectory
    simU = np.zeros((Nsim, nu))    # Control trajectory
    simX[0,:] = x0                 # Initial state
    t = np.zeros((Nsim))           # Computation times
```

#### The Closed-Loop
The main loop implements the MPC feedback:

```python
for i in range(Nsim):
    # 1. Solve OCP for current state
    simU[i,:] = ocp_solver.solve_for_x0(x0_bar=simX[i, :])
    
    # 2. Record computation time
    t[i] = ocp_solver.get_stats('time_tot')
    
    # 3. Simulate system with computed control
    simX[i+1, :] = integrator.simulate(x=simX[i, :], u=simU[i,:])
```

**Key points:**
- `solve_for_x0()` updates the initial constraint `x0` and solves the OCP, returning the first control input
- The OCP predicts `N_horizon` steps ahead, but only the first control is applied
- The integrator simulates the **actual** system response to that control
- The loop continues with the new measured state

This is called **receding horizon control** - the optimization horizon moves forward at each step.

### 5. Analyze Results

After the simulation, the example evaluates performance:

```python
# Timing statistics
t *= 1000  # Convert to milliseconds
print(f'Computation time in ms: min {np.min(t):.3f} median {np.median(t):.3f} max {np.max(t):.3f}')

# Plot trajectories
model = ocp_solver.acados_ocp.model
plot_pendulum(np.linspace(0, (Tf/N_horizon)*Nsim, Nsim+1), 
              Fmax, simU, simX, 
              latexify=False,
              time_label=model.t_label, 
              x_labels=model.x_labels, 
              u_labels=model.u_labels)
```

The plots show:
- State trajectories: cart position, pendulum angle, velocities
- Control input: force over time
- The pendulum successfully swings up and stabilizes

## Key Concepts

### Receding Horizon Principle
At each time step:
1. Measure current state
2. Solve an OCP over a future horizon
3. Apply only the first control action
4. Move one step forward and repeat

This provides feedback to handle disturbances and model mismatch.

### Separation of OCP and Simulation
The example uses two separate components:
- **OCP Solver**: Assumes a model to predict and optimize
- **Integrator**: Simulates the "true" system

In real applications, the integrator is replaced by the actual physical system. Using separate solvers in simulation helps test robustness.

### Code Generation
`acados` generates optimized C code for the specific problem:
- The OCP structure is exploited for efficiency
- Generated code can be deployed on embedded systems
- Python provides a convenient interface to the generated solver

## Running the Example

1. Ensure `acados` is installed with Python interface
2. Navigate to the example directory:
   ```bash
   cd <acados_root>/examples/acados_python/getting_started
   ```
3. Run the script:
   ```bash
   python minimal_example_closed_loop.py
   ```

The script will:
- Generate and compile C code (first run takes longer)
- Run the closed-loop simulation
- Display computation time statistics
- Show plots of the trajectories

## Next Steps

After understanding this example, you can:
- Modify the cost function weights to change behavior
- Adjust the prediction horizon and see the effect on performance
- Add state constraints (e.g., limit cart position)
- Try different initial conditions
- Explore the [advanced-step RTI example](https://github.com/acados/acados/blob/main/examples/acados_python/pendulum_on_cart/as_rti/as_rti_closed_loop_example.py) for real-time implementations
- Check out other examples in the [features documentation](../features/index.md)

## Related Examples

- **Open-loop OCP**: [minimal_example_ocp.py](https://github.com/acados/acados/blob/main/examples/acados_python/getting_started/minimal_example_ocp.py) - solve a single OCP without feedback
- **Integrator**: [minimal_example_sim.py](https://github.com/acados/acados/blob/main/examples/acados_python/getting_started/minimal_example_sim.py) - standalone simulation
- **Real-time MPC**: [as_rti_closed_loop_example.py](https://github.com/acados/acados/blob/main/examples/acados_python/pendulum_on_cart/as_rti/as_rti_closed_loop_example.py) - advanced real-time iteration scheme

## Further Reading

- [Python Interface Documentation](./index.md)
- [Problem Formulation](../problem_formulation/index.md)
- [acados Paper](https://www.sciencedirect.com/science/article/pii/S0005109819303620)
