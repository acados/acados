## Roadmap

#### documentation
- [ ] provide OCP NLP formulation that is handled by `ocp_nlp` as a formula in docs
    - [ ] closely stick to setter names!
- [ ] Set up and document binary workflow
    - [ ] Windows Matlab
    - [ ] MacOS Matlab

#### `matlab interface`
- [ ] separate `acados_ocp()` into generating the C object and setting the numerical data
- [ ] code generation workflow!
- [ ] support nonuniform grids?
- [ ] OCP with DAEs

#### `templating`
- [ ] add GNSF

#### `general`
- [ ] RTI implementation similar to ACADO

#### `ocp_nlp`
- [x] Gauss-Newton SQP
- [x] exact Hessian SQP
- [ ] partial tightening <!-- - [ ] HPNMPC (what?!) -->
- [ ] blockSQP

#### `ocp_qp`
- [x] qpOASES v3.1
- [x] OOQP
- [x] qpDUNES
- [x] HPMPC
- [x] OSQP

#### `sim`
- [x] explicit Runge-Kutta
- [x] lifted IRK
- [x] collocation integrators GL
- [ ] collocation integrators Radau
- [ ] GNSF Hessians