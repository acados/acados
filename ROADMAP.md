## Roadmap
- [ ] Test 32 bit
- [ ] Test OOQP
- [ ] get_optimal_value_hessian() - at least for HPIPM and exact Hessians
- [ ] remove old layer in MATLAB interface

#### core
- [x] propagate cost in integrator for NLS+IRK
    - or: add support for quadrature state, separate dimension in integrator and OCP solver
- [ ] faster workspace memory casting


#### `sim`
- [ ] GNSF Hessians
- [x] propagate cost in integrator for CONL+IRK
- [ ] time in integrator + time dependent model functions


## DONE
- [x] closed loop example MPC + MHE
- [x] Templates: avoid global memory
- [x] add support for manual model functions -- partly done: external cost and discrete dynamics
- [x] More flexible solution sensitivities

#### `ocp_nlp`
- [x] partial tightening -> implemented via multi-phase
- [x] RTI implementation similar to ACADO
- [x] support cost on z for external, NLS

#### C
- [x] split ocp solve into prepare and feedback

#### build
- [x] cmake: add openmp parallelization

#### matlab interface
- [x] detect dimensions
- [x] detect slack dimensions
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
