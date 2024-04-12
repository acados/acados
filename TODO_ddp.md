## TODO List for DDP solver in acados:

1. Miscellaneous:
- [x] Sanity check for unconstrained OCP, HPIPM no condensing (OJ)
- [] Check iteration counter and statistics size, if they are correctly implemented

2. Globalization:
- [X] Check feasibility of initial iterate and make it feasible if necessary
- [X] Merit function only with objective
- [X] Line search
- [X] Regularization with Levenberg-Marquardt term
- [ ] Check if NAN in forward dynamics simulation
- [X] For the moment remove time step weighing in Levenberg-Marquardt

3. Feasibility Problem Translation:
- [ ] Parse rockit problems to acados directly in feasibility form (D)
- [x] Parse acados problems to feasibility problem, i.e., move constraints into objective (OJ)
- [] Check parser such that all combinations of constraints: bx, bu, bxe, NL constraints are correctly parsed and no errors are thrown
- [X] Weights for constraints in objective per N

4. Testing and ideally same OR better results as python Implementation
- [x] Add QP problem, then SQP and DDP with full step should be equivalent (D)
- [X] Add NL_LS problem with linear dynamics, then SQP and DDP with full step using GN Hessian should be equivalent (D)
- [ ] add test in CI https://github.com/FreyJo/acados/blob/8f5ee7cf4c8d5a688c83f5c705fd5bbbf733f9ec/interfaces/CMakeLists.txt
- [x] swingup example (OJ)
