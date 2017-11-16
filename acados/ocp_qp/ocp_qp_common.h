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

#ifndef ACADOS_OCP_QP_OCP_QP_COMMON_H_
#define ACADOS_OCP_QP_OCP_QP_COMMON_H_

#ifdef __cplusplus
extern "C" {
#endif

// hpipm
#include "hpipm/include/hpipm_d_ocp_qp.h"
#include "hpipm/include/hpipm_d_ocp_qp_sol.h"
#include "hpipm/include/hpipm_d_ocp_qp_dim.h"
// acados
#include "acados/utils/types.h"


typedef struct d_ocp_qp_dim ocp_qp_dims;
typedef struct d_ocp_qp ocp_qp_in;
typedef struct d_ocp_qp_sol ocp_qp_out;




// NOTE(dimitris): contains both ocp_qp solvers and condensing with dense solvers
typedef enum {
    HPIPM,
    CONDENSING_HPIPM,
    CONDENSING_QPOASES
} qp_solver_t;


typedef struct {
    int (*fun)(ocp_qp_in *qp_in, ocp_qp_out *qp_out, void *args, void *mem);
    int (*calculate_args_size)(ocp_qp_dims *dims);
    void *(*assign_args)(ocp_qp_dims *dims, void *raw_memory);
    void (*initialize_default_args)(void *args);
    int (*calculate_memory_size)(ocp_qp_dims *dims, void *args);
    void *(*assign_memory)(ocp_qp_dims *dims, void *args, void *raw_memory);
} ocp_qp_solver;


typedef struct {
    int (*fun)(ocp_qp_in *qp_in, ocp_qp_out *qp_out, void *args, void *mem);
    int (*calculate_args_size)(ocp_qp_dims *dims, qp_solver_t solver_name);
    void *(*assign_args)(ocp_qp_dims *dims, qp_solver_t solver_name, void *raw_memory);
    void (*initialize_default_args)(void *args);
    int (*calculate_memory_size)(ocp_qp_dims *dims, void *args);
    void *(*assign_memory)(ocp_qp_dims *dims, void *args, void *raw_memory);
} ocp_qp_xcond_solver;


//
int ocp_qp_in_calculate_size(ocp_qp_dims *dims);
//
ocp_qp_in *assign_ocp_qp_in(ocp_qp_dims *dims, void *raw_memory);
//
int ocp_qp_out_calculate_size(ocp_qp_dims *dims);
//
ocp_qp_out *assign_ocp_qp_out(ocp_qp_dims *dims, void *raw_memory);
//
// TODO(dimitris):remove
ocp_qp_solver initialize_ocp_qp_solver(qp_solver_t qp_solver_name);
//
void set_qp_solver_fun_ptrs(qp_solver_t qp_solver_name, void *qp_solver);
//
void set_xcond_qp_solver_fun_ptrs(qp_solver_t qp_solver_name, ocp_qp_xcond_solver *qp_solver);

// TODO TEMP
// void form_nbu_nbx_rev(int N, int *nbu, int *nbx, int *nb, int* nx, int *nu, int **idxb_rev);
// void ocp_qp_in_copy_dynamics(const real_t *A, const real_t *B, const real_t *b, ocp_qp_in *qp_in,
//     int_t stage);
// void ocp_qp_in_copy_objective(const real_t *Q, const real_t *S, const real_t *R, const real_t *q,
//      const real_t *r, ocp_qp_in *qp_in, int_t stage);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif // ACADOS_OCP_QP_OCP_QP_COMMON_H_
