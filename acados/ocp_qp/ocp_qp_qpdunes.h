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

#ifdef ACADOS_WITH_QPDUNES

#ifndef ACADOS_OCP_QP_OCP_QP_QPDUNES_H_
#define ACADOS_OCP_QP_OCP_QP_QPDUNES_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "qpDUNES.h"

#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/utils/types.h"

typedef enum qpdunes_options_t_ {
    QPDUNES_DEFAULT_ARGUMENTS,
    QPDUNES_LINEAR_MPC,    // TODO(dimitris): implement
    QPDUNES_NONLINEAR_MPC  // TODO(dimitris): implement
} qpdunes_options_t;


typedef enum {
    QPDUNES_WITH_QPOASES,
    QPDUNES_WITH_CLIPPING
} qpdunes_stage_qp_solver_t;


typedef struct ocp_qp_qpdunes_args_ {
    qpOptions_t options;
    qpdunes_stage_qp_solver_t stageQpSolver;
    bool isLinearMPC;
} ocp_qp_qpdunes_args;


typedef struct ocp_qp_qpdunes_memory_ {
    int firstRun;
    int nx;
    int nu;
    int nz;
    int nDmax;  // max(dims->ng)
    qpData_t qpData;
} ocp_qp_qpdunes_memory;


typedef struct ocp_qp_qpdunes_workspace_ {
    double *H;
    double *Q;
    double *R;
    double *S;
    double *g;
    double *ABt;
    double *b;
    double *Ct;
    double *lc;
    double *uc;
    double *zLow;
    double *zUpp;
} ocp_qp_qpdunes_workspace;


//
int ocp_qp_qpdunes_calculate_args_size(ocp_qp_dims *dims);
//
void *ocp_qp_qpdunes_assign_args(ocp_qp_dims *dims, void *raw_memory);
//
void ocp_qp_qpdunes_initialize_default_args(void *args_);
//
int ocp_qp_qpdunes_calculate_memory_size(ocp_qp_dims *dims, void *args_);
//
void *ocp_qp_qpdunes_assign_memory(ocp_qp_dims *dims, void *args_, void *raw_memory);
//
int ocp_qp_qpdunes_calculate_workspace_size(ocp_qp_dims *dims, void *args_);
//
int ocp_qp_qpdunes(ocp_qp_in *qp_in, ocp_qp_out *qp_out, void *args_, void *memory_, void *work_);
//
void ocp_qp_qpdunes_free_memory(void *mem_);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_OCP_QP_OCP_QP_QPDUNES_H_

#endif  // ACADOS_WITH_QPDUNES