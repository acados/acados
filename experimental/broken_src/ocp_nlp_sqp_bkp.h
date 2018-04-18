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

#ifndef ACADOS_OCP_NLP_OCP_NLP_SQP_H_
#define ACADOS_OCP_NLP_OCP_NLP_SQP_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "acados/ocp_nlp/ocp_nlp_common.h"
#include "acados/ocp_nlp/ocp_nlp_sm_common.h"
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/utils/types.h"

typedef struct {
    int_t maxIter;
    ocp_qp_solver *qp_solver;
    ocp_nlp_sm *sensitivity_method;
    // char qp_solver_name[MAX_STR_LEN];
    // char sm_method_name[MAX_STR_LEN];
} ocp_nlp_sqp_args;

typedef struct {
    ocp_nlp_memory *common;

    // TODO(nielsvd): qp solver and sensitivity method
    //       could be automatically initialized for user
    //       convenience, look into this!
    ocp_nlp_sm_in *sm_in;
    ocp_nlp_sm_out *sm_out;
    // ocp_qp_solver *qp_solver;
    // ocp_nlp_sm *sensitivity_method;
} ocp_nlp_sqp_memory;

typedef struct {
    int_t dummy;
} ocp_nlp_sqp_workspace;

ocp_nlp_sqp_args *ocp_nlp_sqp_create_arguments();

int_t ocp_nlp_sqp_calculate_memory_size(const ocp_nlp_in *nlp_in, void *args_);

char *ocp_nlp_sqp_assign_memory(const ocp_nlp_in *nlp_in, void *args_, void **mem_,
                                void *raw_memory);

ocp_nlp_sqp_memory *ocp_nlp_sqp_create_memory(const ocp_nlp_in *nlp_in, void *args_);

int_t ocp_nlp_sqp_calculate_workspace_size(const ocp_nlp_in *nlp_in, void *args_);

char *ocp_nlp_sqp_assign_workspace(const ocp_nlp_in *nlp_in, void *args_, void **work_,
                                   void *raw_memory);

ocp_nlp_sqp_workspace *ocp_nlp_sqp_create_workspace(const ocp_nlp_in *nlp_in, void *args_);

int_t ocp_nlp_sqp(const ocp_nlp_in *nlp_in, ocp_nlp_out *nlp_out, void *args_, void *memory_,
                  void *workspace_);

void ocp_nlp_sqp_initialize(const ocp_nlp_in *nlp_in, void *args_, void **mem_, void **work_);

void ocp_nlp_sqp_destroy(void *mem_, void *work_);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_OCP_NLP_OCP_NLP_SQP_H_
