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

#ifndef ACADOS_SIM_SIM_LIFTED_IRK_INTEGRATOR_H_
#define ACADOS_SIM_SIM_LIFTED_IRK_INTEGRATOR_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "acados/sim/sim_collocation.h"
#include "acados/sim/sim_rk_common.h"
#include "acados/utils/types.h"

#define TRIPLE_LOOP 1
#define CODE_GENERATION 0

typedef struct {
    real_t *rhs_in;
    real_t *jac_tmp;
    real_t **VDE_tmp;
    real_t *out_tmp;
    int_t *ipiv;

    real_t *sys_mat;
    real_t *sys_sol;
    real_t *sys_sol_trans;

    real_t *trans;
    struct d_strmat *str_mat;
    struct d_strmat *str_sol;

    real_t *out_adj_tmp;
} sim_lifted_irk_workspace;

typedef struct {
    real_t *K_traj;
    real_t *DK_traj;
    real_t *delta_DK_traj;
    real_t *mu_traj;

    real_t **sys_mat2;
    real_t **sys_sol2;
    struct d_strmat **str_mat2;
    struct d_strmat **str_sol2;
    int_t **ipiv2;
    real_t *adj_traj;

    real_t **jac_traj;

    real_t *x;
    real_t *u;
} sim_lifted_irk_memory;

int_t sim_lifted_irk(const sim_in *in, sim_out *out, void *args, void *mem,
                     void *work);

int_t sim_lifted_irk_calculate_workspace_size(const sim_in *in, void *args);

void sim_lifted_irk_create_memory(const sim_in *in, void *args,
                                  sim_lifted_irk_memory *mem);
void sim_lifted_irk_free_memory(void *mem_);

void sim_irk_create_arguments(void *args, const int_t num_stages, const char* name);

void sim_lifted_irk_initialize(const sim_in *in, void *args_, void *mem_,
                               void **work);
void sim_lifted_irk_destroy(void *mem, void *work);

void sim_irk_control_collocation(void *args, int_t num_stages,
                                 const char *name);

void sim_irk_create_Newton_scheme(void *args, int_t num_stages,
                                  const char *name,
                                  enum Newton_type_collocation type);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_SIM_SIM_LIFTED_IRK_INTEGRATOR_H_
