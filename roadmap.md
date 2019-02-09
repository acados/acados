## Roadmap

#### `ocp_nlp`
- [x] Gauss-Newton SQP
- [ ] exact Hessian SQP
- [ ] partial tightening
- [ ] HPNMPC
- [ ] blockSQP

#### `ocp_qp`
- [x] qpOASES v1.0
- [ ] qpOASES v2.0 (no dynamic memory allocation, same structure as rest of QP solvers)
- [ ] block condensing (should maybe get its own category)
- [x] OOQP
- [x] qpDUNES
- [x] HPMPC
- [ ] FORCES
- [x] OSQP

#### `sim`
- [x] explicit Runge-Kutta
- [x] lifted IRK
- [x] collocation integrators GL
- [ ] collocation integrators Radau
- [x] discrete-time systems
- [x] second order sensitivities

#### `general`
- [ ] allow for models with varying state and control dimensions
- [ ] unit tests in Python (`unit_test` or `nosetests`) or MATLAB (built-in framework)
- [ ] 'code generation' of driver files such that users don't have to code them manually
- [ ] short Python doctests to exemplify usage of interface functions
- [ ] installation via conda and/or pip
- [ ] RTI implementation similar to ACADO
- [ ] Provide Simulink blocks as an interface (self-written or Legacy Code Tool)
- [ ] HPMPC and HPIPM for MSVC / dSpace
- [ ] Proper installation instructions/Troubleshooting on the Wiki/Readme
- [ ] Installation instructions for MacOS
- [ ] Provide separate documentation page
