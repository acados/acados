## TODO List for DDP solver in acados:

1. Miscellaneous:
- [x] Sanity check for unconstrained OCP, HPIPM no condensing (OJ)

2. Globalization:
- [ ] Check feasibility of initial iterate and make it feasible if necessary
- [ ] Merit function only with objective
- [ ] Line search
- [ ] Regularization with Levenberg-Marquardt term
- [ ] Check if NAN in forward dynamics simulation

3. Feasibility Problem Translation:
- [ ] Parse rockit problems to acados directly in feasibility form (D)
- [ ] Parse acados problems to feasibility problem, i.e., move constraints into objective (OJ)
- [ ] Weights for constraints in objective per N

4. Testing and ideally same OR better results as python Implementation
- [ ] Add QP problem, then SQP and DDP with full step should be equivalent (D)
- [ ] add test in CI https://github.com/FreyJo/acados/blob/8f5ee7cf4c8d5a688c83f5c705fd5bbbf733f9ec/interfaces/CMakeLists.txt
- [ ] swingup example (OJ)
