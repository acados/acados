#ifndef EXAMPLES_C_QUADCOPTER_MODEL_QUADCOPTER_MODEL_H_
#define EXAMPLES_C_QUADCOPTER_MODEL_QUADCOPTER_MODEL_H_

#include "acados/utils/types.h"

int vdeFun(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);

void VDE_fun_quadcopter(const real_t* in, real_t* out,
                      int (*vde)(const real_t**, real_t**, int*, real_t*, int));

void jac_fun_quadcopter(const real_t* in, real_t* out);

#endif  // EXAMPLES_C_QUADCOPTER_MODEL_QUADCOPTER_MODEL_H_
