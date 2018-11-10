Integrator Interface:
=============

Complete integrator interface
==
- [] add cpp generation for all integrator function
    - [] ERK
        - [x] expl_ode_fun
        - [close] expl_ode_hes -- how to transpose SX matrix
        - [x] expl_vde_for
        - [x] expl_vde_adj
    - [] IRK
        - [x] impl_ode_fun;
        - [x] impl_ode_fun_jac_x_xdot_z;
        - [x] impl_ode_jac_x_xdot_u_z;
        - [] impl_ode_hess; - how to transpose?!
    - [x] LIFTED_IRK
        - [x] impl_ode_fun;
        - [x] impl_ode_fun_jac_x_xdot_u
    - [] GNSF?!
- [] add logic, which functions to generate for which settings
    - [x] ERK
    - [x] IRK
- [] extend to non native model, i.e. IRK with expl model
- [] add support for parameters in cpp function generation


After DEMO
- [J] IRK with EXPLICIT
- [] output dict/struct
- [x] Fix ocp_nlp
- [JT] Cleanup & Refactor
- [] .get_settings

Done
=====
TODOs:
- [X] 2 integrators at the same time
- [X] minimal MX support
- [X] IRK support
- [X] Matlab support
- [T] try to compile on Windows

For Presenting
- [X] Minimal Python
- [X] Mighty Python example with timings
- [X] Matlab example with timings
- [T] Design choices