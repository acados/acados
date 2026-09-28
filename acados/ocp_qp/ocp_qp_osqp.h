/*
 * Copyright (c) The acados authors.
 *
 * This file is part of acados.
 *
 * Licensed under the 2-Clause BSD License.
 */


#ifndef ACADOS_OCP_QP_OCP_QP_OSQP_H_
#define ACADOS_OCP_QP_OCP_QP_OSQP_H_

#ifdef __cplusplus
extern "C" {
#endif

// osqp
#include "osqp/include/public/osqp.h"

// acados
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/utils/types.h"

typedef struct ocp_qp_osqp_opts_
{
    OSQPSettings *osqp_opts;
    int print_level;
} ocp_qp_osqp_opts;


typedef struct ocp_qp_osqp_memory_
{
    OSQPInt first_run;

    OSQPFloat *q;
    OSQPFloat *l;
    OSQPFloat *u;

    OSQPInt P_nnzmax;
    OSQPInt *P_i;
    OSQPInt *P_p;
    OSQPFloat *P_x;

    OSQPInt A_nnzmax;
    OSQPInt *A_i;
    OSQPInt *A_p;
    OSQPFloat *A_x;

    OSQPCscMatrix *P;
    OSQPCscMatrix *A;
    OSQPSolver *osqp_solver;

    double time_qp_solver_call;
    int iter;
    int status;

} ocp_qp_osqp_memory;

acados_size_t ocp_qp_osqp_opts_calculate_size(void *config, void *dims);
//
void *ocp_qp_osqp_opts_assign(void *config, void *dims, void *raw_memory);
//
void ocp_qp_osqp_opts_initialize_default(void *config, void *dims, void *opts_);
//
void ocp_qp_osqp_opts_update(void *config, void *dims, void *opts_);
//
acados_size_t ocp_qp_osqp_memory_calculate_size(void *config, void *dims, void *opts_);
//
void *ocp_qp_osqp_memory_assign(void *config, void *dims, void *opts_, void *raw_memory);
//
acados_size_t ocp_qp_osqp_workspace_calculate_size(void *config, void *dims, void *opts_);
//
int ocp_qp_osqp(void *config, void *qp_in, void *qp_out, void *opts_, void *mem_, void *work_);
//
void ocp_qp_osqp_memory_reset(void *config_, void *qp_in_, void *qp_out_, void *opts_, void *mem_, void *work_);
//
void ocp_qp_osqp_solver_get(void *config_, void *qp_in_, void *qp_out_, void *opts_, void *mem_, const char *field, int stage, void* value, int size1, int size2);
//
void ocp_qp_osqp_config_initialize_default(void *config);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_OCP_QP_OCP_QP_OSQP_H_
