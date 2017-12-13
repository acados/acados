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

#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/utils/types.h"

typedef struct ocp_qp_ooqp_args_ {
    int_t printLevel;
    int_t useDiagonalWeights;  // TODO(dimitris): implement option
    int_t fixHessian;
    int_t fixHessianSparsity;
    int_t fixDynamics;
    int_t fixDynamicsSparsity;
    int_t fixInequalities;
    int_t fixInequalitiesSparsity;
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
    int_t *tmpInt;    // temporary vector to sort indicies sparse matrices
    real_t *tmpReal;  // temporary vector to sort data of sparse matrices
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
    int_t nnz;  // max(nnzQ, nnzA, nnzC)
} ocp_qp_ooqp_memory;

ocp_qp_ooqp_args *ocp_qp_ooqp_create_arguments();

ocp_qp_ooqp_memory *ocp_qp_ooqp_create_memory(const ocp_qp_in *input, void *args_);

void ocp_qp_ooqp_free_memory(void *mem_);

int_t ocp_qp_ooqp_calculate_workspace_size(const ocp_qp_in *in, void *args_);

int_t ocp_qp_ooqp(const ocp_qp_in *input, ocp_qp_out *output, void *args_, void *memory_,
                  void *work_);

void ocp_qp_ooqp_initialize(const ocp_qp_in *qp_in, void *args_, void **mem, void **work);

void ocp_qp_ooqp_destroy(void *mem_, void *work);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_OCP_QP_OCP_QP_OOQP_H_
