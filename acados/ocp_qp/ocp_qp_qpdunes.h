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

#include "qpDUNES.h"

#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/utils/types.h"

typedef enum qpdunes_options_t_ {
    QPDUNES_DEFAULT_ARGUMENTS,
    QPDUNES_LINEAR_MPC,     // TODO(dimitris): partly implemented
    QPDUNES_NONLINEAR_MPC,  // TODO(dimitris): not implemented yet
    QPDUNES_ACADO_SETTINGS
} qpdunes_options_t;

typedef enum { QPDUNES_WITH_QPOASES, QPDUNES_WITH_CLIPPING } qpdunes_stage_qp_solver_t;

typedef struct ocp_qp_qpdunes_opts_
{
    qpOptions_t options;
    qpdunes_stage_qp_solver_t stageQpSolver;
    int warmstart;  // warmstart = 0: all multipliers set to zero, warmstart = 1: use previous mult.
    bool isLinearMPC;
} ocp_qp_qpdunes_opts;

typedef struct ocp_qp_qpdunes_memory_
{
    int firstRun;
    int nx;
    int nu;
    int nz;
    int nDmax;  // max(dims->ng)
    qpData_t qpData;
} ocp_qp_qpdunes_memory;

typedef struct ocp_qp_qpdunes_workspace_
{
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
int ocp_qp_qpdunes_opts_calculate_size(void *config_, ocp_qp_dims *dims);
//
void *ocp_qp_qpdunes_opts_assign(void *config_, ocp_qp_dims *dims, void *raw_memory);
//
void ocp_qp_qpdunes_opts_initialize_default(void *config_, ocp_qp_dims *dims, void *opts_);
//
void ocp_qp_qpdunes_opts_update(void *config_, ocp_qp_dims *dims, void *opts_);
//
int ocp_qp_qpdunes_memory_calculate_size(void *config_, ocp_qp_dims *dims, void *opts_);
//
void *ocp_qp_qpdunes_memory_assign(void *config_, ocp_qp_dims *dims, void *opts_, void *raw_memory);
//
int ocp_qp_qpdunes_workspace_calculate_size(void *config_, ocp_qp_dims *dims, void *opts_);
//
int ocp_qp_qpdunes(void *config_, ocp_qp_in *qp_in, ocp_qp_out *qp_out, void *opts_, void *memory_,
                   void *work_);
//
void ocp_qp_qpdunes_free_memory(void *mem_);
//
void ocp_qp_qpdunes_config_initialize_default(void *config_);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_OCP_QP_OCP_QP_QPDUNES_H_
