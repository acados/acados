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


#ifndef ACADOS_UTILS_CREATE_H_
#define ACADOS_UTILS_CREATE_H_

#ifdef __cplusplus
extern "C" {
#endif


#include "acados/dense_qp/dense_qp_common.h"
#include "acados/dense_qp/dense_qp_qpoases.h"
#include "acados/dense_qp/dense_qp_hpipm.h"
#include "acados/ocp_nlp/ocp_nlp_common.h"
#include "acados/ocp_qp/ocp_qp_condensing.h"
#include "acados/ocp_qp/ocp_qp_condensing_hpipm.h"
#include "acados/ocp_qp/ocp_qp_condensing_qpoases.h"
#include "acados/ocp_qp/ocp_qp_partial_condensing.h"
#include "acados/ocp_qp/ocp_qp_sparse_solver.h"
#include "acados/ocp_qp/ocp_qp_condensing_solver.h"
#include "acados/ocp_qp/ocp_qp_hpipm.h"

ocp_qp_in *create_ocp_qp_in(ocp_qp_dims *dims);

ocp_qp_out *create_ocp_qp_out(ocp_qp_dims *dims);

ocp_qp_hpipm_args *ocp_qp_hpipm_create_arguments(ocp_qp_dims *dims);

ocp_qp_hpipm_memory *ocp_qp_hpipm_create_memory(ocp_qp_dims *dims, void *args_);

ocp_qp_condensing_args *ocp_qp_condensing_create_arguments(ocp_qp_dims *dims);

ocp_qp_condensing_memory *ocp_qp_condensing_create_memory(ocp_qp_dims *dims, ocp_qp_condensing_args *args);

ocp_qp_partial_condensing_args *ocp_qp_partial_condensing_create_arguments(ocp_qp_dims *dims);

ocp_qp_partial_condensing_memory *ocp_qp_partial_condensing_create_memory(ocp_qp_dims *dims, ocp_qp_partial_condensing_args *args);

ocp_qp_sparse_solver_args *ocp_qp_sparse_solver_create_arguments(ocp_qp_dims *dims, qp_solver_t solver_name);

ocp_qp_sparse_solver_memory *ocp_qp_sparse_solver_create_memory(ocp_qp_dims *dims, void *args_);

ocp_qp_condensing_solver_args *ocp_qp_condensing_solver_create_arguments(ocp_qp_dims *dims, qp_solver_t solver_name);

ocp_qp_condensing_solver_memory *ocp_qp_condensing_solver_create_memory(ocp_qp_dims *dims, void *args_);

ocp_qp_condensing_hpipm_args *ocp_qp_condensing_hpipm_create_arguments(ocp_qp_dims *dims);

ocp_qp_condensing_hpipm_memory *ocp_qp_condensing_hpipm_create_memory(ocp_qp_dims *dims, void *args_);

ocp_qp_condensing_qpoases_args *ocp_qp_condensing_qpoases_create_arguments(ocp_qp_dims *dims);

ocp_qp_condensing_qpoases_memory *ocp_qp_condensing_qpoases_create_memory(ocp_qp_dims *dims, void *args_);

dense_qp_in *create_dense_qp_in(dense_qp_dims *dims);

dense_qp_out *create_dense_qp_out(dense_qp_dims *dims);

dense_qp_hpipm_args *dense_qp_hpipm_create_arguments(dense_qp_dims *dims);

dense_qp_hpipm_memory *dense_qp_hpipm_create_memory(dense_qp_dims *dims, void *args_);

dense_qp_qpoases_args *dense_qp_qpoases_create_arguments(dense_qp_dims *dims);

dense_qp_qpoases_memory *dense_qp_qpoases_create_memory(dense_qp_dims *dims, void *args_);

ocp_nlp_in *create_ocp_nlp_in(ocp_nlp_dims *dims, int num_stages);



#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_UTILS_CREATE_H_
