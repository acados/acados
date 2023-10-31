# About
This folder contains examples of a fast implementation of the zero-order robust optimization (zoRO) algorithm [1].
The disturbance propagation and the constraint tightening of the zoRO algorithm are in a template based `C` function.
The function can be called through the *custom_update* interface in acados.

In Python, one can define the disturbance/uncertainty matrix and set which constraints to be tightened. The robustification of the optimal control problem against disturbance/uncertainty is achieved by calling the *custom_update()* function before solving the MPC problem.

## How to use it
In the Python file, define the filenames of the custom update function:
```
ocp.solver_options.custom_update_filename = 'custom_update_function.c'
    ocp.solver_options.custom_update_header_filename = 'custom_update_function.h'

ocp.solver_options.custom_update_copy = False
ocp.solver_options.custom_templates = [
    ('custom_update_function_zoro_template.in.c', 'custom_update_function.c'),
    ('custom_update_function_zoro_template.in.h', 'custom_update_function.h'),
]
```

Then setup the `ZoroDescription`. For all properties, check the documentation of this class.
These settings include the initial uncertainty matrix $P_0$ (P0_mat), the feedback matrix $K$ (fdbk_K_mat), the $W$ matrix (W_mat), and the sensitivity of the discretized dynamics with respect to the noise, i.e. the $G$ matrix (unc_jac_G_mat), the scaling of the backoff terms, which relates to the probability level of constraint satisfaction.

For example:
```
zoro_description = ZoroDescription()
zoro_description.fdbk_K_mat = np.zeros((nu, nx))
zoro_description.P0_mat = np.zeros((nx, nx))
zoro_description.W_mat = np.eye(nx)
```
Setting the constraints to be tightened can be done using the index sets idx_lbx_t (the index of the lower bound state constraints to be tightened), idx_lbx_e_t, idx_lbu, idx_lg_t, idx_lg_e_t, idx_lh_t, idx_lh_e_t, etc.
E.g.
```
zoro_description.idx_lbu_t = [0]
```

To perform the uncertainty propagation and update the backoff terms of the OCP, run
```
ocp_solver.custom_update([])
```

The initial disturbance matrix $P_0$ can be passed as a function argument:
```
ocp_solver.custom_update([P0_mat.flatten()])
```

## Examples
A minimal example can be found in *pendulum_on_cart/minimal_example_zoro.py*.
Other examples include the continuous stirred-tank reactor, and the differential drive robot and are also in this folder.

[1] Andrea Zanelli, Jonathan Frey, Florian Messerer, Moritz Diehl, Zero-Order Robust Nonlinear Model Predictive Control with Ellipsoidal Uncertainty Sets, IFAC-PapersOnLine,
Volume 54, Issue 6, 2021, Pages 50-57, ISSN 2405-8963, https://doi.org/10.1016/j.ifacol.2021.08.523.
