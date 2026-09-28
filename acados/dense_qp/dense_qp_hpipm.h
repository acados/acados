/*
 * Copyright (c) The acados authors.
 *
 * This file is part of acados.
 *
 * Licensed under the 2-Clause BSD License.
 */


#ifndef ACADOS_DENSE_QP_DENSE_QP_HPIPM_H_
#define ACADOS_DENSE_QP_DENSE_QP_HPIPM_H_

#ifdef __cplusplus
extern "C" {
#endif

// hpipm
#include "hpipm/include/hpipm_d_dense_qp.h"
#include "hpipm/include/hpipm_d_dense_qp_ipm.h"
#include "hpipm/include/hpipm_d_dense_qp_sol.h"
// acados
#include "acados/dense_qp/dense_qp_common.h"
#include "acados/utils/types.h"



typedef struct dense_qp_hpipm_opts_
{
    struct d_dense_qp_ipm_arg *hpipm_opts;
    double m_relax;
    int print_level;
} dense_qp_hpipm_opts;



typedef struct dense_qp_hpipm_memory_
{
    struct d_dense_qp_ipm_ws *hpipm_workspace;
    double time_qp_solver_call;
    int iter;

} dense_qp_hpipm_memory;



//
acados_size_t dense_qp_hpipm_opts_calculate_size(void *config, void *dims);
//
void *dense_qp_hpipm_opts_assign(void *config, void *dims, void *raw_memory);
//
void dense_qp_hpipm_opts_initialize_default(void *config, void *dims, void *opts_);
//
void dense_qp_hpipm_opts_update(void *config, void *dims, void *opts_);
//
acados_size_t dense_qp_hpipm_memory_calculate_size(void *config_, void *dims_, void *opts_);
//
void *dense_qp_hpipm_memory_assign(void *config_, void *dims_, void *opts_, void *raw_memory);
//
acados_size_t dense_qp_hpipm_calculate_workspace_size(void *dims, void *opts_);
//
int dense_qp_hpipm(void *config, void *qp_in, void *qp_out, void *opts_, void *mem_, void *work_);
//
void dense_qp_hpipm_config_initialize_default(void *config_);
//
void dense_qp_hpipm_memory_reset(void *config, void *qp_in, void *qp_out, void *opts, void *mem, void *work);
//
void dense_qp_hpipm_solver_get(void *config_, void *qp_in_, void *qp_out_, void *opts_, void *mem_, const char *field, int stage, void* value, int size1, int size2);


#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_DENSE_QP_DENSE_QP_HPIPM_H_
