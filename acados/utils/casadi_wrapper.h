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

#ifndef ACADOS_UTILS_CASADI_WRAPPER_H_
#define ACADOS_UTILS_CASADI_WRAPPER_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "acados/utils/types.h"

typedef struct {
    const real_t *x;
    const real_t *u;
    const real_t *p;

    bool compute_jac;
    bool compute_hess;
} casadi_wrapper_in;

typedef struct {
    real_t *y;
    real_t *jac_y;
    real_t *hess_y;
} casadi_wrapper_out;

typedef struct {
    int_t (*fun)(const real_t **, real_t **, int_t *, real_t *, int_t);
    int_t (*dims)(int_t *, int_t *, int_t *, int_t *);
    const int_t *(*sparsity)(int_t);
} casadi_wrapper_args;

typedef struct {
    const real_t **arg;
    real_t **res;
    int_t *iw;
    real_t *w;
    real_t **sparse_res;
} casadi_wrapper_workspace;

casadi_wrapper_args *casadi_wrapper_create_arguments();

int_t casadi_wrapper_calculate_workspace_size(const casadi_wrapper_in *cw_in,
                                              casadi_wrapper_args *args);

char *casadi_wrapper_assign_workspace(const casadi_wrapper_in *cw_in,
                                      casadi_wrapper_args *args,
                                      casadi_wrapper_workspace **work,
                                      void *raw_memory);

casadi_wrapper_workspace *casadi_wrapper_create_workspace(
    const casadi_wrapper_in *cw_in, casadi_wrapper_args *args);

int_t casadi_wrapper(const casadi_wrapper_in *cw_in, casadi_wrapper_out *cw_out,
                     casadi_wrapper_args *args, casadi_wrapper_workspace *work);

void casadi_wrapper_initialize(const casadi_wrapper_in *cw_in,
                               casadi_wrapper_args *args,
                               casadi_wrapper_workspace **work);

void casadi_wrapper_destroy(casadi_wrapper_workspace *work);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_UTILS_CASADI_WRAPPER_H_
