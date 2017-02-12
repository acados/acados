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

#include "acados/utils/types.h"
#include "acados/ocp_qp/ocp_qp_common.h"
#include "qpDUNES/include/qpDUNES.h"

typedef enum qpdunes_options_t_ {
    QPDUNES_DEFAULT_ARGUMENTS,
    QPDUNES_LINEAR_MPC,
    QPDUNES_NONLINEAR_MPC
} qpdunes_options_t;

typedef enum qpdunes_stage_qp_solver_t_ {
    QPDUNES_WITH_QPOASES,
    QPDUNES_WITH_CLIPPING
} qpdunes_stage_qp_solver_t;

typedef struct ocp_qp_qpdunes_args_ {
    qpOptions_t options;
} ocp_qp_qpdunes_args;

typedef struct ocp_qp_qpdunes_workspace_ {
    real_t *At;
    real_t *Bt;
    real_t *scrap;
    real_t *zLow;
    real_t *zUpp;
    real_t *g;
    int tmp;
} ocp_qp_qpdunes_workspace;

typedef struct ocp_qp_qpdunes_memory_ {
    int_t firstRun;
    int_t dimA;
    int_t dimB;
    int_t maxDim;
    int_t dimz;
    qpData_t qpData;
    qpdunes_stage_qp_solver_t stageQpSolver;
} ocp_qp_qpdunes_memory;

int_t ocp_qp_qpdunes_create_arguments(void *args_, int_t opts_);
int_t ocp_qp_qpdunes_calculate_workspace_size(const ocp_qp_in *in, void *args_);
int_t ocp_qp_qpdunes_create_memory(const ocp_qp_in *input, void *args_, void *memory_);

void ocp_qp_qpdunes_free_memory(void *mem_);
// void ocp_qp_qpdunes_free_workspace(void *work_);

int_t ocp_qp_qpdunes(ocp_qp_in *input, ocp_qp_out *output,
    void *args_, void *memory_, void *work_);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_OCP_QP_OCP_QP_QPDUNES_H_
