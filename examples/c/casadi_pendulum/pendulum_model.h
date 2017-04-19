#ifndef EXAMPLES_CASADI_PENDULUM_PENDULUM_MODEL_H_
#define EXAMPLES_CASADI_PENDULUM_PENDULUM_MODEL_H_

#include "acados/utils/types.h"

void VDE_fun_pendulum(const real_t* in, real_t* out, 
    int (*vde)(const real_t**, real_t**, int*, real_t*, int));

void jac_fun_pendulum(const real_t* in, real_t* out);

#endif  // EXAMPLES_CASADI_PENDULUM_PENDULUM_MODEL_H_
