## Roadmap

#### core
- [ ] propagate cost in integrator
- [ ] restore download and compilation of OOQP
- [ ] split ocp solve into prepare and feedback
- [ ] stage transition functions for changing model dimensions

#### documentation
- [ ] provide OCP NLP formulation that is handled by `ocp_nlp` as a formula in docs
    - [ ] closely stick to setter names!
- [ ] Set up and document binary workflow
    - [ ] Windows Matlab
    - [ ] MacOS Matlab

#### `matlab interface`
- [x] detect dimensions
- [x] structure detections for constraints
- [ ] code generation workflow! # Mex templating
    - 1) nonlinear least-squares 
    - 2) Vz (already implemented?)
    - 3) exact Hessian
    - 4) external cost
    - 5) GNSF
    - 6) discrete dynamics
- [x] separate `acados_ocp()` into generating the C object and setting the numerical data
- [x] support nonuniform grids
- [x] OCP with DAEs

#### build
- [ ] cmake: add openmp parallelization

#### `templating`
- [ ] remove duplicated code to generate external functions
- [ ] add documentation
  - [ ] how to set up trenderer
  - [ ] ...
- [ ] explore possibility to interact with generated C code
- [ ] set up CI tests - matlab needed..
- [ ] add GNSF
- [ ] make compatible with octave (if possible) - jsonencode alternative?!

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
