/*
 * Copyright (c) The acados authors.
 *
 * This file is part of acados.
 *
 * Licensed under the 2-Clause BSD License.
 */


#ifndef ACADOS_DENSE_QP_DENSE_QP_DAQP_H_
#define ACADOS_DENSE_QP_DENSE_QP_DAQP_H_

#ifdef __cplusplus
extern "C" {
#endif

// blasfeo
#include "blasfeo_common.h"

// daqp
#include "daqp/include/types.h"

// acados
#include "acados/dense_qp/dense_qp_common.h"
#include "acados/utils/types.h"


typedef struct dense_qp_daqp_opts_
{
    DAQPSettings* daqp_opts;
    int warm_start;
    int print_level;
} dense_qp_daqp_opts;


typedef struct dense_qp_daqp_memory_
{
    double* blower;
    double* bupper;
    int* idxs;

    double* Zl;
    double* Zu;
    double* zl;
    double* zu;
    double* d_ls;
    double* d_us;

    double time_qp_solver_call;
    int iter;
    int matrices_initialized;
    DAQPWorkspace * daqp_work;
    struct blasfeo_dmat *H_factor;
    struct blasfeo_dmat *M_factor;
    struct blasfeo_dvec *rhs_factor;
    struct blasfeo_dvec *v_factor;
    struct blasfeo_dvec *constraint_value;
    double *ldp_rows;  // packed rows of Rinv for the simple bounds, followed by M

} dense_qp_daqp_memory;

// opts
acados_size_t dense_qp_daqp_opts_calculate_size(void *config, dense_qp_dims *dims);
//
void *dense_qp_daqp_opts_assign(void *config, dense_qp_dims *dims, void *raw_memory);
//
void dense_qp_daqp_opts_initialize_default(void *config, dense_qp_dims *dims, void *opts_);
//
void dense_qp_daqp_opts_update(void *config, dense_qp_dims *dims, void *opts_);
//
// memory
acados_size_t dense_qp_daqp_workspace_calculate_size(void *config, dense_qp_dims *dims, void *opts_);
//
void *dense_qp_daqp_workspace_assign(void *config, dense_qp_dims *dims, void *raw_memory);
//
acados_size_t dense_qp_daqp_memory_calculate_size(void *config, dense_qp_dims *dims, void *opts_);
//
void *dense_qp_daqp_memory_assign(void *config, dense_qp_dims *dims, void *opts_, void *raw_memory);
//
// functions
int dense_qp_daqp(void *config, dense_qp_in *qp_in, dense_qp_out *qp_out, void *opts_, void *memory_, void *work_);
//
void dense_qp_daqp_memory_reset(void *config_, void *qp_in, void *qp_out, void *opts_, void *mem_, void *work_);
//
void dense_qp_daqp_config_initialize_default(void *config_);
//
void dense_qp_daqp_memory_reset(void *config, void *qp_in, void *qp_out, void *opts, void *mem, void *work);
//
void dense_qp_daqp_solver_get(void *config_, void *qp_in_, void *qp_out_, void *opts_, void *mem_, const char *field, int stage, void* value, int size1, int size2);


#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_DENSE_QP_DENSE_QP_DAQP_H_
