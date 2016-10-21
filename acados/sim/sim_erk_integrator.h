#ifndef ACADOS_SIM_SIM_ERK_INTEGRATOR_H_
#define ACADOS_SIM_SIM_ERK_INTEGRATOR_H_

#include "acados/utils/types.h"
#include "acados/sim/sim_rk_common.h"

typedef struct sim_erk_workspace_ {
    real_t *rhs_in;
    real_t *K_tmp;
    real_t *out_tmp;
} sim_erk_workspace;


void sim_erk(const sim_in *in, sim_out *out, const sim_RK_opts *opts, sim_erk_workspace *work);

void sim_erk_create_workspace(const sim_in *in, sim_RK_opts *opts, sim_erk_workspace *work);

void sim_erk_create_opts(int_t num_stages, sim_RK_opts *opts);

#endif  // ACADOS_SIM_SIM_ERK_INTEGRATOR_H_
