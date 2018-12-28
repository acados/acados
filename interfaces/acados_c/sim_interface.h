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
    sim_solver_config *config;
    void *dims;
    void *opts;
    void *mem;
    void *work;
} sim_solver;

/* config */
//
sim_solver_config *sim_config_create(sim_solver_plan plan);
//
void sim_config_free(void *config);

/* dims */
//
void *sim_dims_create(void *config_);
//
void sim_dims_free(void *dims);
//
void sim_dims_set(sim_solver_config *config, void *dims, const char *field, const int* value);

/* in */
//
sim_in *sim_in_create(sim_solver_config *config, void *dims);
//
void sim_in_free(void *out);
//
void sim_in_set_T(sim_solver_config *config, double T, sim_in *in);
//
int sim_model_set(sim_solver_config *config, sim_in *in, const char *fun_type, void *fun_ptr);
//
int sim_model_set_internal(sim_solver_config *config, void *model, const char *fun_type,
                           void *fun_ptr);
//
void sim_in_set_x(sim_solver_config *config, void *dims, double *x, sim_in *in);
//
void sim_in_set_xdot(sim_solver_config *config, void *dims, double *xdot, sim_in *in);
//
void sim_in_set_u(sim_solver_config *config, void *dims, double *u, sim_in *in);
//
void sim_in_set_Sx(sim_solver_config *config, void *dims, double *Sx, sim_in *in);
//
void sim_in_set_Su(sim_solver_config *config, void *dims, double *Su, sim_in *in);

/* out */
//
sim_out *sim_out_create(sim_solver_config *config, void *dims);
//
void sim_out_free(void *out);
//
void sim_out_get_xn(sim_solver_config *config, void *dims, sim_out *out, double *xn);
//
void sim_out_get_Sxn(sim_solver_config *config, void *dims, sim_out *out, double *Sxn);
//
void sim_out_get_Sun(sim_solver_config *config, void *dims, sim_out *out, double *Sun);

/* opts */
//
void *sim_opts_create(sim_solver_config *config, void *dims);
//
void sim_opts_free(void *opts);
//
int sim_opts_set(sim_solver_config *config, void *opts, const char *field,
                           void *value);
/* solver */
//
int sim_calculate_size(sim_solver_config *config, void *dims, void *opts_);
//
sim_solver *sim_assign(sim_solver_config *config, void *dims, void *opts_, void *raw_memory);
//
sim_solver *sim_create(sim_solver_config *config, void *dims, void *opts_);
//
void sim_free(void *solver);
//
int sim_solve(sim_solver *solver, sim_in *in, sim_out *out);
//
int sim_precompute(sim_solver *solver, sim_in *in, sim_out *out);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // INTERFACES_ACADOS_C_SIM_INTERFACE_H_
