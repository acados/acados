/*
 *    This file is part of acados.
 *
 *    acados is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    acados is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with acados; if not, write to the Free Software Foundation,
 *    Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */
/* OLD implementation removed from acados kernel, since there is a new implementation,
    which is easier to maintain - FreyJo */
#ifndef ACADOS_SIM_SIM_LIFTED_IRK_INTEGRATOR_H_
#define ACADOS_SIM_SIM_LIFTED_IRK_INTEGRATOR_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "acados/sim/sim_collocation_utils.h"
#include "acados/sim/sim_common.h"
#include "acados/utils/types.h"

#define TRIPLE_LOOP 1
#define CODE_GENERATION 0

typedef struct
{
    int nx;
    int nu;
    int nz;
} sim_lifted_irk_dims;

typedef struct
{
    /* external functions */
    // jacobian explicit ode
    external_function_generic *expl_ode_jac;
    // forward explicit vde
    external_function_generic *expl_vde_for;

} lifted_irk_model;

typedef struct
{
    real_t *rhs_in;
    real_t *jac_tmp;
    real_t **VDE_tmp;
    real_t *out_tmp;
    int_t *ipiv;

    real_t *sys_mat;
    real_t *sys_sol;
    real_t *sys_sol_trans;

    real_t *trans;
    struct blasfeo_dmat *str_mat;
    struct blasfeo_dmat *str_sol;

    real_t *out_adj_tmp;
} sim_lifted_irk_workspace;

typedef struct
{
    double *grad_correction;
    double *grad_K;  // gradient correction

    real_t *K_traj;
    real_t *DK_traj;
    real_t *mu_traj;

    real_t *x;
    real_t *u;

    real_t *delta_DK_traj;
    real_t *adj_traj;
    real_t **jac_traj;

    real_t **sys_mat2;
    int_t **ipiv2;
    real_t **sys_sol2;

    struct blasfeo_dmat **str_mat2;
    struct blasfeo_dmat **str_sol2;

} sim_lifted_irk_memory;

// get & set functions
void sim_lifted_irk_set_nx(void *dims_, int nx);
void sim_lifted_irk_set_nu(void *dims_, int nu);
void sim_lifted_irk_set_nz(void *dims_, int nz);

void sim_lifted_irk_get_nx(void *dims_, int *nx);
void sim_lifted_irk_get_nu(void *dims_, int *nu);
void sim_lifted_irk_get_nz(void *dims_, int *nz);
//
int sim_lifted_irk_dims_calculate_size();
//
void *sim_lifted_irk_dims_assign(void *config_, void *raw_memory);
//
int sim_lifted_irk_model_calculate_size(void *config, void *dims);
//
void *sim_lifted_irk_model_assign(void *config, void *dims, void *raw_memory);
//
int sim_lifted_irk_model_set_function(void *model_, sim_function_t fun_type, void *fun);
//
int sim_lifted_irk_opts_calculate_size(void *config, void *dims);
//
void *sim_lifted_irk_opts_assign(void *config, void *dims, void *raw_memory);
//
void sim_lifted_irk_opts_initialize_default(void *config, void *dims, void *opts_);
//
void sim_lifted_irk_opts_update(void *config_, void *dims, void *opts_);
//
int sim_lifted_irk_memory_calculate_size(void *config, void *dims, void *opts);
//
void *sim_lifted_irk_memory_assign(void *config, void *dims, void *opts_, void *raw_memory);
//
int sim_lifted_irk_workspace_calculate_size(void *config, void *in, void *args);
//
int sim_lifted_irk(void *config, sim_in *in, sim_out *out, void *args, void *mem, void *work);
//
void sim_lifted_irk_config_initialize_default(void *config);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_SIM_SIM_LIFTED_IRK_INTEGRATOR_H_
