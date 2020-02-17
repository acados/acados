## Roadmap

#### core
- [ ] propagate cost in integrator
- [ ] restore download and compilation of OOQP
- [ ] split ocp solve into prepare and feedback
- [ ] stage transition functions for changing model dimensions


#### documentation
- [x] provide OCP NLP formulation that is handled by `ocp_nlp` as a formula
    - [x] closely stick to setter names!
- [ ] Set up and document binary workflow
    - [x] Windows Matlab
    - [ ] MacOS Matlab

#### `matlab interface`
- [x] detect dimensions
- [x] structure detections for constraints
- [x] getting started folder
- [ ] add Mex templating support for: ( in prioritized order )
    - [ ] nonlinear least-squares
    - [ ] Vz (already implemented?)
    - [ ] exact Hessian
    - [ ] external cost
    - [ ] GNSF
    - [ ] discrete dynamics
- [x] separate `acados_ocp()` into generating the C object and setting the numerical data
- [x] support nonuniform grids
- [x] OCP with DAEs

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
