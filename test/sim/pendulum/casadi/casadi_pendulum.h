#ifndef TEST_SIM_PENDULUM_CASADI_CASADI_PENDULUM_H_
#define TEST_SIM_PENDULUM_CASADI_CASADI_PENDULUM_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "acados/utils/types.h"

int vde_forw_pendulum(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);

void VDE_forw_pendulum(const real_t* in, real_t* out, int (*vde)(const real_t**, real_t**, int*, real_t*, int));

void VDE_adj_pendulum(const real_t* in, real_t* out);

void VDE_hess_pendulum(const real_t* in, real_t* out);

void jac_fun_pendulum(const real_t* in, real_t* out);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // TEST_SIM_PENDULUM_CASADI_CASADI_PENDULUM_H_
