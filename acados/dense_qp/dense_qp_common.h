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
#include "hpipm/include/hpipm_d_dense_qp_sol.h"
#include "hpipm/include/hpipm_d_dense_qp_res.h"
// acados
#include "acados/utils/types.h"

typedef struct d_dense_qp_dim dense_qp_dims;
typedef struct d_dense_qp dense_qp_in;
typedef struct d_dense_qp_sol dense_qp_out;
typedef struct d_dense_qp_res dense_qp_res;
typedef struct d_dense_qp_res_workspace dense_qp_res_ws;



typedef struct {
    int (*fun)(dense_qp_in *qp_in, dense_qp_out *qp_out, void *args, void *mem, void *work);
    int (*calculate_args_size)(dense_qp_dims *dims);
    void *(*assign_args)(dense_qp_dims *dims, void *raw_memory);
    void (*initialize_default_args)(void *args);
    int (*calculate_memory_size)(dense_qp_dims *dims, void *args);
    void *(*assign_memory)(dense_qp_dims *dims, void *args, void *raw_memory);
    int (*calculate_workspace_size)(dense_qp_dims *dims, void *args);
} dense_qp_solver;


typedef struct {
    double solve_QP_time;
    double interface_time;
    double total_time;
} dense_qp_info;


//
int dense_qp_in_calculate_size(dense_qp_dims *dims);
//
dense_qp_in *assign_dense_qp_in(dense_qp_dims *dims, void *raw_memory);
//
int dense_qp_out_calculate_size(dense_qp_dims *dims);
//
dense_qp_out *assign_dense_qp_out(dense_qp_dims *dims, void *raw_memory);
//
int dense_qp_res_calculate_size(dense_qp_dims *dims);
//
dense_qp_res *assign_dense_qp_res(dense_qp_dims *dims, void *raw_memory);
//
int dense_qp_res_ws_calculate_size(dense_qp_dims *dims);
//
dense_qp_res_ws *assign_dense_qp_res_ws(dense_qp_dims *dims, void *raw_memory);
//
void compute_dense_qp_res(dense_qp_in *qp_in, dense_qp_out *qp_out, dense_qp_res *qp_res, dense_qp_res_ws *res_ws);
//
void compute_dense_qp_res_nrm_inf(dense_qp_res *qp_res, double res[4]);
//

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_DENSE_QP_DENSE_QP_COMMON_H_
