/*
 * Copyright (c) The acados authors.
 *
 * This file is part of acados.
 *
 * Licensed under the 2-Clause BSD License.
 */


#ifndef ACADOS_OCP_QP_OCP_QP_CLARABEL_H_
#define ACADOS_OCP_QP_OCP_QP_CLARABEL_H_

#ifdef __cplusplus
extern "C" {
#endif

// clarabel
// #include "clarabel/Clarabel.h"
#include "Clarabel.cpp/include/clarabel.h"

// acados
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/utils/types.h"

typedef struct ocp_qp_clarabel_opts_
{
    // settings *clarabel_opts;
    ClarabelDefaultSettings *clarabel_opts;
    int print_level;
    int first_run;

} ocp_qp_clarabel_opts;


typedef struct ocp_qp_clarabel_memory_
{
    ClarabelCscMatrix P; // just upper triangular is enough
    uintptr_t P_nnzmax;
    uintptr_t *P_col_ptr;
    uintptr_t *P_rowval;
    ClarabelFloat *P_nzval;
    uintptr_t P_nnz;

    ClarabelCscMatrix A;
    uintptr_t A_nnzmax;
    uintptr_t *A_col_ptr;
    uintptr_t *A_rowval;
    ClarabelFloat *A_nzval;
    uintptr_t A_nnz;

    ClarabelFloat *q;
    uintptr_t q_nnz;
    ClarabelFloat *b;
    uintptr_t b_nnz;

    ClarabelSupportedConeT cones[2];

    ClarabelDefaultSolver *solver;
    ClarabelDefaultSolution solution;

    double time_qp_solver_call;
    int iter;
    int status;

} ocp_qp_clarabel_memory;

acados_size_t ocp_qp_clarabel_opts_calculate_size(void *config, void *dims);
//
void *ocp_qp_clarabel_opts_assign(void *config, void *dims, void *raw_memory);
//
void ocp_qp_clarabel_opts_initialize_default(void *config, void *dims, void *opts_);
//
void ocp_qp_clarabel_opts_update(void *config, void *dims, void *opts_);
//
acados_size_t ocp_qp_clarabel_memory_calculate_size(void *config, void *dims, void *opts_);
//
void *ocp_qp_clarabel_memory_assign(void *config, void *dims, void *opts_, void *raw_memory);
//
acados_size_t ocp_qp_clarabel_workspace_calculate_size(void *config, void *dims, void *opts_);
//
int ocp_qp_clarabel(void *config, void *qp_in, void *qp_out, void *opts_, void *mem_, void *work_);
//
void ocp_qp_clarabel_memory_reset(void *config_, void *qp_in_, void *qp_out_, void *opts_, void *mem_, void *work_);
//
void ocp_qp_clarabel_solver_get(void *config_, void *qp_in_, void *qp_out_, void *opts_, void *mem_, const char *field, int stage, void* value, int size1, int size2);
//
void ocp_qp_clarabel_config_initialize_default(void *config);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_OCP_QP_OCP_QP_CLARABEL_H_
