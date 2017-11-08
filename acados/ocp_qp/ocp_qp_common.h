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
#include "hpipm/include/hpipm_d_ocp_qp_size.h"
// acados
#include "acados/utils/types.h"


typedef struct d_ocp_qp_size ocp_qp_dims;


typedef struct d_ocp_qp ocp_qp_in;


typedef struct d_ocp_qp_sol ocp_qp_out;


typedef enum {
    HPIPM,
    CONDENSING_HPIPM,
    CONDENSING_QPOASES
} qp_solver_t;


typedef struct {
    int (*calculate_args_size)(ocp_qp_dims *dims);
    void *(*assign_args)(ocp_qp_dims *dims, void *mem);
    void (*initialize_default_args)(void *args);
    //...
    int (*solve)(ocp_qp_in *qp_in, ocp_qp_out *qp_out, void *args, void *mem, void *work);
    ocp_qp_in *qp_in;
    ocp_qp_out *qp_out;
    void *args;
    void *mem;
    void *work;
} new_ocp_qp_solver;



typedef struct {
    int (*fun)(ocp_qp_in *qp_in, ocp_qp_out *qp_out, void *args, void *mem);
	// TODO remove ???
    void (*initialize)(ocp_qp_in *qp_in, void *args, void **mem, void **work);
	// TODO remove ???
    void (*destroy)(void *mem, void *work);
	// TODO add calculate_size and assign instead ???
	int (*calculate_memory)(ocp_qp_in *qp_in, void *args);
	char* (*assign_memory)(ocp_qp_in *qp_in, void *args, void **mem, void *raw_mem);
    ocp_qp_in *qp_in;
    ocp_qp_out *qp_out;
    void *args;
    void *mem;
    void *work;
} ocp_qp_solver;



//
int ocp_qp_in_calculate_size(ocp_qp_dims *dims);
// TODO make name consistent !!! (ocp_qp_out_assign)
char *assign_ocp_qp_in(ocp_qp_dims *dims, ocp_qp_in **qp_in, void *mem);
//
int ocp_qp_out_calculate_size(ocp_qp_dims *dims);
// TODO make name consistent !!! (ocp_qp_out_assign)
char *assign_ocp_qp_out(ocp_qp_dims *dims, ocp_qp_out **qp_out, void *mem);
//
new_ocp_qp_solver initialize_ocp_qp_solver(qp_solver_t qp_solver_name);

// TODO TEMP

void form_nbu_nbx_rev(int N, int *nbu, int *nbx, int *nb, int* nx, int *nu, int **idxb_rev);

void ocp_qp_in_copy_dynamics(const real_t *A, const real_t *B, const real_t *b, ocp_qp_in *qp_in,
    int_t stage);

void ocp_qp_in_copy_objective(const real_t *Q, const real_t *S, const real_t *R, const real_t *q,
     const real_t *r, ocp_qp_in *qp_in, int_t stage);

ocp_qp_solver *create_ocp_qp_solver(const ocp_qp_in *qp_in, const char *name, void *options);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif // ACADOS_OCP_QP_OCP_QP_COMMON_H_
