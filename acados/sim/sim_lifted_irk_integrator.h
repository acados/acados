#ifndef ACADOS_SIM_SIM_LIFTED_IRK_INTEGRATOR_H_
#define ACADOS_SIM_SIM_LIFTED_IRK_INTEGRATOR_H_

#include "acados/utils/types.h"
#include "acados/sim/sim_rk_common.h"

typedef struct sim_lifted_irk_workspace_ {
    real_t *rhs_in;
    real_t *VDE_tmp;
    real_t *out_tmp;
    int *ipiv;
    real_t *sys_mat;
    real_t *sys_sol;
    struct d_strmat *str_mat;
    struct d_strmat *str_sol;
} sim_lifted_irk_workspace;

typedef struct sim_lifted_irk_memory_ {
    real_t *K_traj;
    real_t *DK_traj;
    real_t *x;
    real_t *u;
} sim_lifted_irk_memory;


void sim_lifted_irk(const sim_in *in, sim_out *out, const sim_RK_opts *opts,
        sim_lifted_irk_memory *mem, sim_lifted_irk_workspace *work);

void sim_lifted_irk_create_workspace(const sim_in *in, sim_RK_opts *opts,
        sim_lifted_irk_workspace *work);

void sim_lifted_irk_create_memory(const sim_in *in, sim_RK_opts *opts,
        sim_lifted_irk_memory *mem);

void sim_irk_create_opts(int_t num_stages, const char* name, sim_RK_opts *opts);

#endif  // ACADOS_SIM_SIM_LIFTED_IRK_INTEGRATOR_H_
