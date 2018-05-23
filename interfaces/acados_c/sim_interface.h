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

#ifndef INTERFACES_ACADOS_C_SIM_INTERFACE_H_
#define INTERFACES_ACADOS_C_SIM_INTERFACE_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "acados/sim/sim_common.h"

typedef enum { ERK, LIFTED_IRK, IRK, GNSF, NEW_LIFTED_IRK } sim_solver_t;

typedef struct
{
    sim_solver_t sim_solver;
} sim_solver_plan;

typedef struct
{
    sim_solver_config *config;
    void *dims;
    void *opts;
    void *mem;
    void *work;
} sim_solver;

//
sim_solver_config *sim_config_create(sim_solver_plan plan);
//
void *sim_dims_create(void *config_);
//
sim_in *sim_in_create(sim_solver_config *config, void *dims);
//
int sim_set_model(sim_solver_config *config, sim_in *in, const char *fun_type, void *fun_ptr);
//
int sim_set_model_internal(sim_solver_config *config, void *model, const char *fun_type,
                           void *fun_ptr);
//
sim_out *sim_out_create(sim_solver_config *config, void *dims);
//
void *sim_opts_create(sim_solver_config *config, void *dims);
//
int sim_calculate_size(sim_solver_config *config, void *dims, void *opts_);
//
sim_solver *sim_assign(sim_solver_config *config, void *dims, void *opts_, void *raw_memory);
//
sim_solver *sim_create(sim_solver_config *config, void *dims, void *opts_);
//
int sim_solve(sim_solver *solver, sim_in *qp_in, sim_out *qp_out);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // INTERFACES_ACADOS_C_SIM_INTERFACE_H_
