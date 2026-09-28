/*
 * Copyright (c) The acados authors.
 *
 * This file is part of acados.
 *
 * Licensed under the 2-Clause BSD License.
 */


#ifndef ACADOS_OCP_QP_OCP_QP_HPIPM_H_
#define ACADOS_OCP_QP_OCP_QP_HPIPM_H_

#ifdef __cplusplus
extern "C" {
#endif

// hpipm
#include "hpipm/include/hpipm_d_ocp_qp_ipm.h"
// acados
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/utils/types.h"



// struct of arguments to the solver
typedef struct ocp_qp_hpipm_opts_
{
    struct d_ocp_qp_ipm_arg *hpipm_opts;
    int print_level;
    double m_relax;
} ocp_qp_hpipm_opts;



// struct of the solver memory
typedef struct ocp_qp_hpipm_memory_
{
    struct d_ocp_qp_ipm_ws *hpipm_workspace;
    double time_qp_solver_call;
    int iter;
    int status;

} ocp_qp_hpipm_memory;



//
acados_size_t ocp_qp_hpipm_opts_calculate_size(void *config, void *dims);
//
void *ocp_qp_hpipm_opts_assign(void *config, void *dims, void *raw_memory);
//
void ocp_qp_hpipm_opts_initialize_default(void *config, void *dims, void *opts_);
//
void ocp_qp_hpipm_opts_update(void *config, void *dims, void *opts_);
//
void ocp_qp_hpipm_opts_set(void *config_, void *opts_, const char *field, void *value);
//
acados_size_t ocp_qp_hpipm_memory_calculate_size(void *config, void *dims, void *opts_);
//
void *ocp_qp_hpipm_memory_assign(void *config, void *dims, void *opts_, void *raw_memory);
//
acados_size_t ocp_qp_hpipm_workspace_calculate_size(void *config, void *dims, void *opts_);
//
int ocp_qp_hpipm(void *config, void *qp_in, void *qp_out, void *opts_, void *mem_, void *work_);
//
void ocp_qp_hpipm_memory_reset(void *config_, void *qp_in_, void *qp_out_, void *opts_, void *mem_, void *work_);
//
void ocp_qp_hpipm_solver_get(void *config_, void *qp_in_, void *qp_out_, void *opts_, void *mem_, const char *field, int stage, void* value, int size1, int size2);
//
void ocp_qp_hpipm_config_initialize_default(void *config);



#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_OCP_QP_OCP_QP_HPIPM_H_
