## Roadmap

#### core
- [ ] propagate cost in integrator
- [ ] restore download and compilation of OOQP

#### documentation
- [ ] provide OCP NLP formulation that is handled by `ocp_nlp` as a formula in docs
    - [ ] closely stick to setter names!
- [ ] Set up and document binary workflow
    - [ ] Windows Matlab
    - [ ] MacOS Matlab

#### `matlab interface`
- [ ] code generation workflow! # Mex templating
- [x] separate `acados_ocp()` into generating the C object and setting the numerical data
- [x] support nonuniform grids
- [x] OCP with DAEs

#### build
- [ ] cmake: add openmp parallelization

#### `templating`
- [ ] add GNSF

#### `ocp_nlp`
- [x] Gauss-Newton SQP
- [x] exact Hessian SQP
- [ ] partial tightening <!-- - [ ] HPNMPC (what?!) -->
- [ ] blockSQP
- [ ] RTI implementation similar to ACADO

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