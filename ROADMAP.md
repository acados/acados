## Roadmap

#### core
- [ ] propagate cost in integrator
- [ ] restore download and compilation of OOQP
- [x] split ocp solve into prepare and feedback
- [ ] stage transition functions for changing model dimensions

#### Python interface
- [ ] lifted IRK
- [x] exact hessian
- [x] regularization
- [x] discrete dynamics

#### matlab interface
- [x] detect dimensions
- [x] structure detections for constraints
- [x] getting started folder
- [x] add Mex templating support for: ( in prioritized order )
    - [x] nonlinear least-squares
    - [ ] Vz (already implemented?)
    - [ ] exact Hessian
    - [x] external cost
    - [ ] GNSF
    - [x] discrete dynamics
- [x] separate `acados_ocp()` into generating the C object and setting the numerical data
- [x] support nonuniform grids
- [x] OCP with DAEs

#### documentation
- [x] provide OCP NLP formulation that is handled by `ocp_nlp` as a formula
    - [x] closely stick to setter names!
- [ ] Set up and document binary workflow
    - [x] Windows Matlab
    - [ ] MacOS Matlab

#### build
- [ ] cmake: add openmp parallelization

#### `ocp_nlp`
- [x] Gauss-Newton SQP
- [x] exact Hessian SQP
- [ ] partial tightening <!-- - [ ] HPNMPC (what?!) -->
- [ ] blockSQP
- [ ] RTI implementation similar to ACADO

#### `sim`
- [ ] collocation integrators Radau
- [ ] GNSF Hessians
