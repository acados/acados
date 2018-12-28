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

#ifndef ACADOS_DENSE_QP_DENSE_QP_COMMON_H_
#define ACADOS_DENSE_QP_DENSE_QP_COMMON_H_

#ifdef __cplusplus
extern "C" {
#endif

// hpipm
#include "hpipm/include/hpipm_d_dense_qp.h"
#include "hpipm/include/hpipm_d_dense_qp_res.h"
#include "hpipm/include/hpipm_d_dense_qp_sol.h"
// acados
#include "acados/utils/types.h"

typedef struct d_dense_qp_dim dense_qp_dims;
typedef struct d_dense_qp dense_qp_in;
typedef struct d_dense_qp_sol dense_qp_out;
typedef struct d_dense_qp_res dense_qp_res;
typedef struct d_dense_qp_res_workspace dense_qp_res_ws;

#ifndef QP_SOLVER_CONFIG_
#define QP_SOLVER_CONFIG_

typedef struct
{
    // TODO(dimitris): pass dims to evaluate?
    void (*dims_set)(void *config_, void *dims_, const char *field, const int* value);
    int (*evaluate)(void *config, void *qp_in, void *qp_out, void *args, void *mem, void *work);
    int (*opts_calculate_size)(void *config, void *dims);
    void *(*opts_assign)(void *config, void *dims, void *raw_memory);
    void (*opts_initialize_default)(void *config, void *dims, void *args);
    void (*opts_update)(void *config, void *dims, void *args);
    int (*memory_calculate_size)(void *config, void *dims, void *args);
    void *(*memory_assign)(void *config, void *dims, void *args, void *raw_memory);
    int (*workspace_calculate_size)(void *config, void *dims, void *args);
} qp_solver_config;

#endif

typedef struct
{
    double solve_QP_time;
    double interface_time;
    double total_time;
    int num_iter;
    int t_computed;
} dense_qp_info;

/* config */
//
int dense_qp_solver_config_calculate_size();
//
qp_solver_config *dense_qp_solver_config_assign(void *raw_memory);

/* dims */
//
int dense_qp_dims_calculate_size();
//
dense_qp_dims *dense_qp_dims_assign(void *raw_memory);
//
void dense_qp_dims_set(void *config_, void *dims_, const char *field, const int* value);
//

/* in */
//
int dense_qp_in_calculate_size(void *config, dense_qp_dims *dims);
//
dense_qp_in *dense_qp_in_assign(void *config, dense_qp_dims *dims, void *raw_memory);
//
int dense_qp_out_calculate_size(void *config, dense_qp_dims *dims);
//
dense_qp_out *dense_qp_out_assign(void *config, dense_qp_dims *dims, void *raw_memory);
//
int dense_qp_res_calculate_size(dense_qp_dims *dims);
//
dense_qp_res *dense_qp_res_assign(dense_qp_dims *dims, void *raw_memory);
//
int dense_qp_res_workspace_calculate_size(dense_qp_dims *dims);
//
dense_qp_res_ws *dense_qp_res_workspace_assign(dense_qp_dims *dims, void *raw_memory);
//
void dense_qp_compute_t(dense_qp_in *qp_in, dense_qp_out *qp_out);
//
void dense_qp_res_compute(dense_qp_in *qp_in, dense_qp_out *qp_out, dense_qp_res *qp_res,
                          dense_qp_res_ws *res_ws);
//
void dense_qp_res_compute_nrm_inf(dense_qp_res *qp_res, double res[4]);
//
void dense_qp_stack_slacks_dims(dense_qp_dims *in, dense_qp_dims *out);
//
void dense_qp_stack_slacks(dense_qp_in *in, dense_qp_in *out);
//
void dense_qp_unstack_slacks(dense_qp_out *in, dense_qp_in *qp_out, dense_qp_out *out);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_DENSE_QP_DENSE_QP_COMMON_H_
