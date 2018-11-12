Integrator Interface:
=============

Complete integrator interface
==
- [] add cpp generation for all integrator function
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
    - [] GNSF?!
- [] add logic, which functions to generate for which settings
    - [x] ERK
    - [x] IRK
    - [] LIFTED IRK
- [] TEST!
    - [] ERK
    - [] IRK
    - [] LIFTED IRK

CURRENT BUGS:
    - ERK with hessian, might be C problem
    - 

FOR LATER
- [] extend to non native model, i.e. IRK with expl model
- [] add support for parameters in cpp function generation


After DEMO
- [] output dict/struct
- [] Cleanup & Refactor
- [x] .get_settings

- [T] try to compile on Windows

Done
=====
- [X] 2 integrators at the same time
- [X] minimal MX support
- [X] IRK support
- [X] Matlab support
- [x] Fix ocp_nlp