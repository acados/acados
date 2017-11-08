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

#ifndef ACADOS_OCP_NLP_OCP_NLP_GN_SQP_H_
#define ACADOS_OCP_NLP_OCP_NLP_GN_SQP_H_

#ifdef __cplusplus
extern "C" {
#endif

// acados
#include "acados/ocp_nlp/ocp_nlp_common.h"
#include "acados/sim/sim_collocation.h"
#include "acados/sim/sim_rk_common.h"
#include "acados/utils/types.h"

typedef struct {
    ocp_nlp_args *common;
    void *qp_solver_args;
} ocp_nlp_gn_sqp_args;

typedef struct {
    ocp_nlp_memory *common;
    ocp_qp_solver *qp_solver;
    ocp_nlp_dims *dims;
} ocp_nlp_gn_sqp_memory;

typedef struct {
    // nlp workspace
    ocp_nlp_work *common;
    // sqp workspace
    ocp_qp_in *qp_in;
    ocp_qp_out *qp_out;
    void *qp_mem;
    void *qp_work;
} ocp_nlp_gn_sqp_work;


int_t ocp_nlp_gn_sqp_calculate_args_size(ocp_nlp_dims *dims, ocp_qp_solver *qp_solver);
//
ocp_nlp_gn_sqp_args *ocp_nlp_gn_sqp_assign_args(ocp_nlp_dims *dims, ocp_qp_solver *qp_solver, void *mem);
//
#if defined(EXT_DEPS)
ocp_nlp_gn_sqp_args *ocp_nlp_gn_sqp_create_args(ocp_nlp_dims *dims, ocp_qp_solver *qp_solver);
#endif

int ocp_nlp_gn_sqp(ocp_nlp_in *nlp_in, ocp_nlp_out *nlp_out, ocp_nlp_gn_sqp_args *args, ocp_nlp_gn_sqp_memory *mem, void *work_);

void ocp_nlp_gn_sqp_create_memory(ocp_nlp_dims *dims, ocp_qp_solver *qp_solver, const ocp_nlp_in *in, void *args_, void *memory_);

void ocp_nlp_gn_sqp_free_memory(void *memory_);

//
int ocp_nlp_gn_sqp_calculate_workspace_size(ocp_nlp_dims *dims, ocp_qp_solver *qp_solver, ocp_nlp_gn_sqp_args *args);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_OCP_NLP_OCP_NLP_GN_SQP_H_
