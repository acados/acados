# Custom Update function feature in acados

This example shows how to implement a custom C function to update solver data and parameters in between solver calls.

## Modify it for your purpose

1. Replace `pendulum_ode` in this example with your `model.name`.

2. Think about what memory you need. You are responsible for handling it!
The acados developer guide contains some info on memory management:
https://docs.acados.org/developer_guide/index.html#memory-management-in-acados

3. Implement your custom update function.

4. Run your fast custom acados application.
