Integrator Interface:
=============

TODO: Complete integrator interface
==

BEFORE MERGE
- [T] output dict/struct
- [] Cleanup & Refactor
- [] add support for parameters in cpp function generation

After MERGE
- [] generalize integrator class
- [T] try to compile on Windows




Done
=====
- [x] fix dependencies
    - [x] travis
- [x] extend to non native model, i.e. IRK with expl model
- [X] 2 integrators at the same time
- [X] minimal MX support
- [X] IRK support
- [X] Matlab support
- [x] Fix ocp_nlp
- [x] .get_settings
- [x] add output check
- [x] move function generation up

- [x] add cpp generation for all integrator function
    - [x] ERK
        - [x] expl_ode_fun
        - [x] expl_ode_hes
        - [x] expl_vde_for
        - [x] expl_vde_adj
    - [x] IRK
        - [x] impl_ode_fun;
        - [x] impl_ode_fun_jac_x_xdot_z;
        - [x] impl_ode_jac_x_xdot_u_z;
        - [x] impl_ode_hess
    - [x] LIFTED_IRK
        - [x] impl_ode_fun;
        - [x] impl_ode_fun_jac_x_xdot_u
    - [x] LIFTED IRK

- [x] add logic, which functions to generate for which settings
    - [x] ERK
    - [x] IRK
    - [x] LIFTED IRK

- [x] TEST!
    - [x] ERK
    - [x] IRK
