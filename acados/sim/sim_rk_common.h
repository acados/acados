#ifndef ACADOS_SIM_SIM_RK_COMMON_H_
#define ACADOS_SIM_SIM_RK_COMMON_H_

#include "acados/utils/types.h"
#include "acados/sim/sim_common.h"

typedef struct sim_RK_opts_ {
    int_t num_stages;
    real_t *A_mat;
    real_t *c_vec;
    real_t *b_vec;
} sim_RK_opts;

#endif  // ACADOS_SIM_SIM_RK_COMMON_H_
