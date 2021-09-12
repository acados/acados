## Roadmap
- [ ] Test 32 bit
- [ ] Test OOQP

#### core
- [ ] propagate cost in integrator
    - or: add support for quadrature state, separate dimension in integrator and OCP solver
- [ ] stage transition functions for changing model dimensions


#### documentation
- [x] provide OCP NLP formulation that is handled by `ocp_nlp` as a formula
    - [x] closely stick to setter names!
- [x] Windows Matlab, reiterate, look into Visual C
- [ ] Windows Python
- [ ] 
- [ ] MacOS Matlab

#### `ocp_nlp`
- [ ] partial tightening <!-- - [ ] HPNMPC (what?!) -->
- [ ] blockSQP (https://github.com/djanka2/blockSQP)
- [x] RTI implementation similar to ACADO
- [ ] support cost on z for external, NLS

#### `sim`
- [ ] collocation integrators Radau
    - NOTE: currently always Gauss(-Legendre) Butcher tables
        - A-stable, but not L-stable
        - order is 2 * num_stages
    - implement also Radau IIA collocation methods
        - L-stable
        - order is 2 * num_stages - 1
- [ ] GNSF Hessians



## DONE
- [x] closed loop example MPC + MHE
- [x] Templates: avoid global memory
- [x] add support for manual model functions -- partly done: external cost and discrete dynamics

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
