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

#ifndef ACADOS_OCP_QP_OCP_QP_HPMPC_H_
#define ACADOS_OCP_QP_OCP_QP_HPMPC_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/utils/types.h"

typedef enum hpmpc_options_t_ {
    HPMPC_DEFAULT_ARGUMENTS,
    HPMPC_PARTIAL_TIGHTENING  
} hpmpc_options_t;

typedef struct ocp_qp_hpmpc_args_ {
    double tol;
    int max_iter;
    //  double min_step;
    double mu0;
    //  double sigma_min;
    int warm_start;
    int N2;  // horizion length of the partially condensed problem
    double **ux0;
    double **pi0;
    double **lam0;
    double **t0;
    int out_iter;          // number of performed iterations
    double *inf_norm_res;  // array of size 5, returning inf norm res

    // partial tightening
    double sigma_mu;
    int N;
    int M;
} ocp_qp_hpmpc_args;

typedef void ocp_qp_hpmpc_memory;  // HPMPC does not have a memory struct

typedef void ocp_qp_hpmpc_workspace;  // // HPMPC does not have a workspace struct

ocp_qp_hpmpc_args *ocp_qp_hpmpc_create_arguments(const ocp_qp_in *qp_in, hpmpc_options_t opts);

int_t ocp_qp_hpmpc_calculate_memory_size(const ocp_qp_in *in, void *args_);

void *ocp_qp_hpmpc_create_memory(const ocp_qp_in *input, void *args_);

void ocp_qp_hpmpc_free_memory(void *mem);

int_t ocp_qp_hpmpc_calculate_workspace_size(const ocp_qp_in *in, void *args);

int_t ocp_qp_hpmpc(const ocp_qp_in *qp_in, ocp_qp_out *qp_out, void *args_, void *mem_,
                   void *workspace_);

// int ocp_qp_hpmpc_libstr_pt(ocp_qp_in *qp_in, ocp_qp_out *qp_out,
//   ocp_qp_hpmpc_args *qp_args, int M, double sigma_mu, void *workspace);

// TODO(Andrea): need to merge hpmpc in order to use this... (Body is ready)
// int ocp_qp_hpnmpc(ocp_qp_in *qp_in, ocp_qp_out *qp_out, ocp_qp_hpmpc_args
// *qp_args,
//   void *workspace);

void ocp_qp_hpmpc_initialize(const ocp_qp_in *qp_in, void *args_, void **mem, void **work);

void ocp_qp_hpmpc_destroy(void *mem, void *work);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_OCP_QP_OCP_QP_HPMPC_H_
