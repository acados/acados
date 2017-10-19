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

#ifndef ACADOS_OCP_QP_OCP_QP_QPDUNES_H_
#define ACADOS_OCP_QP_OCP_QP_QPDUNES_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <stdbool.h>

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
    bool isLinearMPC;
} ocp_qp_qpdunes_args;

typedef struct ocp_qp_qpdunes_memory_ {
    int_t firstRun;
    int_t dimA;
    int_t dimB;
    int_t dimC;    // maximum number of elements of matrix: [Cx Cu]
    int_t maxDim;  // max(dimA+dimB, dimC) needed for scrap (to transpose
                   // matrices)
    int_t dimz;
    int_t nDmax;
    qpData_t qpData;
    qpdunes_stage_qp_solver_t stageQpSolver;
} ocp_qp_qpdunes_memory;

typedef struct ocp_qp_qpdunes_workspace_ {
    real_t *H;
    real_t *g;
    real_t *ABt;
    real_t *Ct;
    real_t *scrap;
    real_t *zLow;
    real_t *zUpp;
    int tmp;  // TODO(dimitris): is this to make sizeof(struct) == 64??
} ocp_qp_qpdunes_workspace;

ocp_qp_qpdunes_args *ocp_qp_qpdunes_create_arguments(qpdunes_options_t opts);

ocp_qp_qpdunes_memory *ocp_qp_qpdunes_create_memory(const ocp_qp_in *in, void *args_);

void ocp_qp_qpdunes_free_memory(void *mem_);

int_t ocp_qp_qpdunes_calculate_workspace_size(const ocp_qp_in *in, void *args_);

int_t ocp_qp_qpdunes(const ocp_qp_in *input, ocp_qp_out *output, void *args_, void *memory_,
                     void *work_);

void ocp_qp_qpdunes_initialize(const ocp_qp_in *qp_in, void *args_, void **mem, void **work);

void ocp_qp_qpdunes_destroy(void *mem, void *work);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_OCP_QP_OCP_QP_QPDUNES_H_
