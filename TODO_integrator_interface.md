Integrator Interface:
=============

Complete integrator interface
==
- [] add cpp generation for all integrator function
    - [] ERK
        - [x] expl_ode_fun
        - [x] expl_ode_jac - REMOVE, it is not used within the implementation!
        - [] expl_ode_hes
        - [x] expl_vde_for
        - [x] expl_vde_adj
    - [] IRK
    - [] NEW_LIFTED_IRK; oj: suggest to skip and hide lifted_irk
    - [] GNSF?!
- [] add logic, which functions to generate for which settings
- [] extend to non native model, i.e. IRK with expl model



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