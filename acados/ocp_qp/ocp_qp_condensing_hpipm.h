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

#ifndef ACADOS_OCP_QP_OCP_QP_CONDENSING_HPIPM_H_
#define ACADOS_OCP_QP_OCP_QP_CONDENSING_HPIPM_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/utils/types.h"

// struct of arguments to the solver
typedef struct ocp_qp_condensing_hpipm_args_ {
    real_t alpha_min;
    real_t res_g_max;
    real_t res_b_max;
    real_t res_d_max;
    real_t res_m_max;
    real_t mu0;
    void *scrapspace;
    int_t iter_max;
} ocp_qp_condensing_hpipm_args;

// struct of the solver memory
typedef struct ocp_qp_condensing_hpipm_memory_ {
    struct d_ocp_qp *qp;
    struct d_ocp_qp_sol *qp_sol;
    struct d_dense_qp *qpd;
    struct d_dense_qp_sol *qpd_sol;
    struct d_cond_qp_ocp2dense_workspace *cond_workspace;
    struct d_dense_qp_ipm_arg *ipm_arg;
    struct d_dense_qp_ipm_workspace *ipm_workspace;
    real_t **hlam_lb;
    real_t **hlam_ub;
    real_t **hlam_lg;
    real_t **hlam_ug;
    int_t **hidxb_rev;
    real_t inf_norm_res[5];
    int_t iter;
} ocp_qp_condensing_hpipm_memory;

ocp_qp_condensing_hpipm_args *ocp_qp_condensing_hpipm_create_arguments(const ocp_qp_in *qp_in);

int_t ocp_qp_condensing_hpipm_calculate_memory_size(const ocp_qp_in *qp_in, void *args_);

char *ocp_qp_condensing_hpipm_assign_memory(const ocp_qp_in *qp_in, void *args_, void **mem_,
                                            void *raw_memory);
//
ocp_qp_condensing_hpipm_memory *ocp_qp_condensing_hpipm_create_memory(const ocp_qp_in *qp_in,
                                                                      void *args_);

int_t ocp_qp_condensing_hpipm_calculate_workspace_size(const ocp_qp_in *qp_in, void *args_);

int_t ocp_qp_condensing_hpipm(const ocp_qp_in *qp_in, ocp_qp_out *qp_out, void *args_,
                              void *memory_, void *workspace_);

void ocp_qp_condensing_hpipm_initialize(const ocp_qp_in *qp_in, void *args_, void **mem,
                                        void **work);

void ocp_qp_condensing_hpipm_destroy(void *mem, void *work);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_OCP_QP_OCP_QP_CONDENSING_HPIPM_H_
