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

#ifndef ACADOS_OCP_QP_OCP_QP_OOQP_H_
#define ACADOS_OCP_QP_OCP_QP_OOQP_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "acados/utils/types.h"
#include "acados/ocp_qp/ocp_qp_common.h"

typedef struct ocp_qp_ooqp_args_ {
    int_t printLevel;
    int_t workspaceMode;        // EXPERIMENTAL: 1 for structs, 2 for chunk of memory
    int_t fixHessian;
    int_t fixHessianSparsity;
    int_t useDiagonalWeights;   // TODO(dimitris): implement option
} ocp_qp_ooqp_args;

typedef struct ocp_qp_ooqp_workspace_ {
    real_t *x;
    real_t *gamma;
    real_t *phi;
    real_t *y;
    real_t *z;
    real_t *lambda;
    real_t *pi;
    real_t objectiveValue;
    // int_t ierr;
} ocp_qp_ooqp_workspace;

typedef struct ocp_qp_ooqp_memory_ {
    int_t firstRun;
    real_t *c;
    int_t nx;
    int_t *irowQ;
    int_t nnzQ;
    int_t *jcolQ;
    int_t *orderQ;
    real_t *dQ;
    real_t *xlow;
    char *ixlow;
    real_t *xupp;
    char *ixupp;
    int_t *irowA;
    int_t nnzA;
    int_t *jcolA;
    int_t *orderA;
    real_t *dA;
    real_t *bA;
    int_t my;
    int_t *irowC;
    int_t nnzC;
    int_t *jcolC;
    int_t *orderC;
    real_t *dC;
    real_t *clow;
    int_t mz;
    char *iclow;
    real_t *cupp;
    char *icupp;
} ocp_qp_ooqp_memory;

int_t ocp_qp_ooqp_create_memory(const ocp_qp_in *input, void *args_, void *memory_);
int_t ocp_qp_ooqp_create_workspace(const ocp_qp_in *input,  void *args_,  void *work_);
int_t ocp_qp_ooqp_calculate_workspace_size(const ocp_qp_in *in, void *args_);

void ocp_qp_ooqp_free_memory(void *mem_);
void ocp_qp_ooqp_free_workspace(void *work_);

int_t ocp_qp_ooqp(ocp_qp_in *input, ocp_qp_out *output,
    void *args_, void *memory_, void *work_);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_OCP_QP_OCP_QP_OOQP_H_
