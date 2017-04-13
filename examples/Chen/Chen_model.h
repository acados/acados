#ifndef EXAMPLES_CHEN_CHEN_MODEL_H_
#define EXAMPLES_CHEN_CHEN_MODEL_H_

#include "acados/utils/types.h"

void VDE_fun(const real_t* in, real_t* out, 
    int (*vde)(const real_t**, real_t**, int*, real_t*, int));

#endif  // EXAMPLES_CHEN_CHEN_MODEL_H_
