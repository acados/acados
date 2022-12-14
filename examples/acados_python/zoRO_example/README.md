# About
Here is the implementation of the zero-order robust optimization (zoRO) algorithm [1]. The disturbance propagation and the constraint tightening of the zoRO algorithm are in the template c function. The function can be called through the *custom_update* interface in Acados.

In Python, one can define the disturbance/uncertainty matrix and set which constraints to be tightened. The robustification of the optimal control problem against disturbance/uncertainty is achieved by calling the *custom_update()* function before solving the MPC problem.

## Disturbance/Uncertainty Propagation
The disturbance/uncertainty propagation is modelled as
$$P_{k+1} = (A_k + B_kK)P_k(A_k + B_kK)^\top + GWG^\top$$

## How to use it



In the python file, define the filenames of the custom update function:
```
ocp.solver_options.custom_update_filename = 'custom_update_function.c'
    ocp.solver_options.custom_update_header_filename = 'custom_update_function.h'

ocp.solver_options.custom_update_copy = False
ocp.solver_options.custom_templates = [
    ('custom_update_function_zoro_template.in.c', 'custom_update_function.c'),
    ('custom_update_function_zoro_template.in.h', 'custom_update_function.h'),
]
```

Then do the zoRO settings. The settings of the disturbance propagation include the initial disturbance matrix $P_0$ (P0_mat), the feedback K matrix $K$ (fdbk_K_mat), the $W$ matrix (W_mat), and the $G$ matrix (unc_jac_G_mat).


```
zoro_description = ZoroDescription()
zoro_description.fdbk_K_mat = np.zeros((nu, nx))
zoro_description.P0_mat = np.zeros((nx, nx))
zoro_description.W_mat = np.eye(nx)
```

The settings of constraints tightening include idx_lbx_t (the index of the lower bound state constraints to be tightened), idx_lbx_e_t, idx_lbu, idx_lg_t, idx_lg_e_t, idx_lh_t, idx_lh_e_t, etc.
```
zoro_description.idx_lbu_t = [0]
```

After all the settings are done, call the *process_zoro_description()* function:
```
ocp.zoro_description = process_zoro_description(zoro_description)
```

Everytime before solving the optimal control problem, run
```
ocp_solver.custom_update([])
```

## Examples
The minimum example can be found in *pendulum_on_cart/minimum_example_zoro.py*. Other examples of the mass chain, the continuous stirred-tank reactor, and the differential drive robot are also in this folder.



[1] Andrea Zanelli, Jonathan Frey, Florian Messerer, Moritz Diehl, Zero-Order Robust Nonlinear Model Predictive Control with Ellipsoidal Uncertainty Sets, IFAC-PapersOnLine,
Volume 54, Issue 6, 2021, Pages 50-57, ISSN 2405-8963, https://doi.org/10.1016/j.ifacol.2021.08.523.
