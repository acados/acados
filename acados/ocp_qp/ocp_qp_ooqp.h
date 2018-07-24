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

enum ocp_qp_ooqp_termination_code
{
  SPARSE_SUCCESSFUL_TERMINATION = 0,
  SPARSE_NOT_FINISHED,
  SPARSE_MAX_ITS_EXCEEDED,
  SPARSE_INFEASIBLE,
  SPARSE_UNKNOWN
};

typedef struct ocp_qp_ooqp_opts_
{
    int printLevel;
    int useDiagonalWeights;  // TODO(dimitris): implement option
    int fixHessian;
    int fixHessianSparsity;
    int fixDynamics;
    int fixDynamicsSparsity;
    int fixInequalities;
    int fixInequalitiesSparsity;
} ocp_qp_ooqp_opts;

typedef struct ocp_qp_ooqp_workspace_
{
    double *x;
    double *gamma;
    double *phi;
    double *y;
    double *z;
    double *lambda;
    double *pi;
    double objectiveValue;
    int *tmpInt;    // temporary vector to sort indicies sparse matrices
    double *tmpReal;  // temporary vector to sort data of sparse matrices
    // int ierr;
} ocp_qp_ooqp_workspace;

typedef struct ocp_qp_ooqp_memory_
{
    int firstRun;
    double *c;
    int nx;
    int *irowQ;
    int nnzQ;
    int *jcolQ;
    int *orderQ;
    double *dQ;
    double *xlow;
    char *ixlow;
    double *xupp;
    char *ixupp;
    int *irowA;
    int nnzA;
    int *jcolA;
    int *orderA;
    double *dA;
    double *bA;
    int my;
    int *irowC;
    int nnzC;
    int *jcolC;
    int *orderC;
    double *dC;
    double *clow;
    int mz;
    char *iclow;
    double *cupp;
    char *icupp;
    int nnz;  // max(nnzQ, nnzA, nnzC)
} ocp_qp_ooqp_memory;

//
int ocp_qp_ooqp_opts_calculate_size(void *config_, ocp_qp_dims *dims);
//
void *ocp_qp_ooqp_opts_assign(void *config_, ocp_qp_dims *dims, void *raw_memory);
//
void ocp_qp_ooqp_opts_initialize_default(void *config_, ocp_qp_dims *dims, void *opts_);
//
void ocp_qp_ooqp_opts_update(void *config_, ocp_qp_dims *dims, void *opts_);
//
int ocp_qp_ooqp_memory_calculate_size(void *config_, ocp_qp_dims *dims, void *opts_);
//
void *ocp_qp_ooqp_memory_assign(void *config_, ocp_qp_dims *dims, void *opts_, void *raw_memory);
//
int ocp_qp_ooqp_workspace_calculate_size(void *config_, ocp_qp_dims *dims, void *opts_);
//
int ocp_qp_ooqp(void *config_, ocp_qp_in *qp_in, ocp_qp_out *qp_out, void *opts_, void *memory_,
                void *work_);
//
void ocp_qp_ooqp_destroy(void *mem_, void *work);
//
void ocp_qp_ooqp_config_initialize_default(void *config_);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_OCP_QP_OCP_QP_OOQP_H_
