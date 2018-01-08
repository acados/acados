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
#include "hpipm/include/hpipm_d_ocp_qp_res.h"
// acados
#include "acados/utils/types.h"


typedef struct d_ocp_qp_dim ocp_qp_dims;
typedef struct d_ocp_qp ocp_qp_in;
typedef struct d_ocp_qp_sol ocp_qp_out;
typedef struct d_ocp_qp_res ocp_qp_res;
typedef struct d_ocp_qp_res_workspace ocp_qp_res_ws;




// NOTE(dimitris): contains both ocp_qp solvers and condensing with dense solvers
typedef enum {
    HPIPM,
    HPMPC,
    CONDENSING_HPIPM,
    CONDENSING_QPOASES,
    CONDENSING_QORE
} qp_solver_t;


typedef struct {
    int (*fun)(ocp_qp_in *qp_in, ocp_qp_out *qp_out, void *args, void *mem, void *work);
    int (*calculate_args_size)(ocp_qp_dims *dims);
    void *(*assign_args)(ocp_qp_dims *dims, void *raw_memory);
    void (*initialize_default_args)(void *args);
    int (*calculate_memory_size)(ocp_qp_dims *dims, void *args);
    void *(*assign_memory)(ocp_qp_dims *dims, void *args, void *raw_memory);
    int (*calculate_workspace_size)(ocp_qp_dims *dims, void *args);
} ocp_qp_solver;


typedef struct {
    int (*fun)(ocp_qp_in *qp_in, ocp_qp_out *qp_out, void *args, void *mem, void *work);
    int (*calculate_args_size)(ocp_qp_dims *dims, void *solver_funs);
    void *(*assign_args)(ocp_qp_dims *dims, void *solver_funs, void *raw_memory);
    void (*initialize_default_args)(void *args);
    int (*calculate_memory_size)(ocp_qp_dims *dims, void *args);
    void *(*assign_memory)(ocp_qp_dims *dims, void *args, void *raw_memory);
    int (*calculate_workspace_size)(ocp_qp_dims *dims, void *args);
    void *qp_solver_funs;
} ocp_qp_xcond_solver;



typedef struct {
    double solve_QP_time;
    double condensing_time;
    double interface_time;
    double total_time;
} ocp_qp_info;


//
int ocp_qp_in_calculate_size(ocp_qp_dims *dims);
//
ocp_qp_in *assign_ocp_qp_in(ocp_qp_dims *dims, void *raw_memory);
//
int ocp_qp_out_calculate_size(ocp_qp_dims *dims);
//
ocp_qp_out *assign_ocp_qp_out(ocp_qp_dims *dims, void *raw_memory);
//
int ocp_qp_res_calculate_size(ocp_qp_dims *dims);
//
ocp_qp_res *assign_ocp_qp_res(ocp_qp_dims *dims, void *raw_memory);
//
int ocp_qp_res_ws_calculate_size(ocp_qp_dims *dims);
//
ocp_qp_res_ws *assign_ocp_qp_res_ws(ocp_qp_dims *dims, void *raw_memory);
//
void compute_ocp_qp_res(ocp_qp_in *qp_in, ocp_qp_out *qp_out, ocp_qp_res *qp_res, ocp_qp_res_ws *res_ws);
//
void compute_ocp_qp_res_nrm_inf(ocp_qp_res *qp_res, double res[4]);
//
int set_qp_solver_fun_ptrs(qp_solver_t qp_solver_name, void *qp_solver);
//
void set_xcond_qp_solver_fun_ptrs(qp_solver_t qp_solver_name, ocp_qp_xcond_solver *qp_solver);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif // ACADOS_OCP_QP_OCP_QP_COMMON_H_
