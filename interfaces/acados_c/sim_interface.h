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

typedef enum { ERK, IRK, GNSF, LIFTED_IRK } sim_solver_t;

typedef struct
{
    sim_solver_t sim_solver;
} sim_solver_plan;

typedef struct
{
    sim_config *config;
    void *dims;
    void *opts;
    void *mem;
    void *work;
} sim_solver;

/* config */
//
sim_config *sim_config_create(sim_solver_plan plan);
//
void sim_config_destroy(void *config);

/* dims */
//
void *sim_dims_create(void *config_);
//
void sim_dims_destroy(void *dims);
//
void sim_dims_set(sim_config *config, void *dims, const char *field, const int* value);
//
void sim_dims_get(sim_config *config, void *dims, const char *field, int* value);

/* in */
//
sim_in *sim_in_create(sim_config *config, void *dims);
//
void sim_in_destroy(void *out);
//
int sim_in_set(void *config_, void *dims_, sim_in *in, const char *field, void *value);


/* out */
//
sim_out *sim_out_create(sim_config *config, void *dims);
//
void sim_out_destroy(void *out);
//
int sim_out_get(void *config, void *dims, sim_out *out, const char *field, void *value);

/* opts */
//
void *sim_opts_create(sim_config *config, void *dims);
//
void sim_opts_destroy(void *opts);
//
int sim_opts_set(sim_config *config, void *opts, const char *field,
                           void *value);
/* solver */
//
int sim_calculate_size(sim_config *config, void *dims, void *opts_);
//
sim_solver *sim_assign(sim_config *config, void *dims, void *opts_, void *raw_memory);
//
sim_solver *sim_solver_create(sim_config *config, void *dims, void *opts_);
//
void sim_solver_destroy(void *solver);
//
int sim_solve(sim_solver *solver, sim_in *in, sim_out *out);
//
int sim_precompute(sim_solver *solver, sim_in *in, sim_out *out);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // INTERFACES_ACADOS_C_SIM_INTERFACE_H_
