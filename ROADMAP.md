## Roadmap
- [ ] Templates: avoid global memory!

#### core
- [ ] propagate cost in integrator
- [ ] restore download and compilation of OOQP
- [ ] stage transition functions for changing model dimensions

#### interfaces
- [ ] Lifted integrators
- [ ] add support for manual model functions

#### documentation
- [x] provide OCP NLP formulation that is handled by `ocp_nlp` as a formula
    - [x] closely stick to setter names!
- [ ] Set up and document binary workflow
    - [ ] Windows Matlab, reiterate, look into Visual C
    - [ ] MacOS Matlab

#### `ocp_nlp`
- [ ] partial tightening <!-- - [ ] HPNMPC (what?!) -->
- [ ] blockSQP (https://github.com/djanka2/blockSQP)
- [ ] RTI implementation similar to ACADO

#### `sim`
- [ ] collocation integrators Radau
- [ ] GNSF Hessians



## DONE
- [x] closed loop example MPC + MHE

#### C
- [x] split ocp solve into prepare and feedback

#### build
- [x] cmake: add openmp parallelization

#### matlab interface
- [x] detect dimensions
- [ ] detect slack dimensions
- [x] structure detections for constraints
- [x] getting started folder
- [x] add Mex templating support for: ( in prioritized order )
    - [x] nonlinear least-squares
    - [x] Vz
    - [x] exact Hessian
    - [x] external cost
    - [x] GNSF
    - [x] discrete dynamics
- [x] separate `acados_ocp()` into generating the C object and setting the numerical data
- [x] support nonuniform grids
- [x] OCP with DAEs

#### Python interface
- [x] exact hessian
- [x] regularization
- [x] discrete dynamics
