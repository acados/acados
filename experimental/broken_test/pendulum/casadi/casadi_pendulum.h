#ifndef TEST_SIM_PENDULUM_CASADI_CASADI_PENDULUM_H_
#define TEST_SIM_PENDULUM_CASADI_CASADI_PENDULUM_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "acados/utils/types.h"

int vde_forw_pendulum(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);

int jac_pendulum(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);

int vde_adj_pendulum(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);

int vde_hess_pendulum(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // TEST_SIM_PENDULUM_CASADI_CASADI_PENDULUM_H_
