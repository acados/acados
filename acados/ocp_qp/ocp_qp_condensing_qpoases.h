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

#ifndef ACADOS_OCP_QP_OCP_QP_CONDENSING_QPOASES_H_
#define ACADOS_OCP_QP_OCP_QP_CONDENSING_QPOASES_H_

#ifdef __cplusplus
extern "C" {
#endif

// acados
#include "acados/ocp_qp/ocp_qp_condensing.h"
#include "acados/dense_qp/dense_qp_qpoases.h"
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/dense_qp/dense_qp_common.h"
#include "acados/utils/types.h"


// struct of solver arguments
typedef struct ocp_qp_condensing_qpoases_args_ {
    ocp_qp_condensing_args *cond_args;
    dense_qp_qpoases_args *solver_args;
} ocp_qp_condensing_qpoases_args;


// struct of the solver memory
typedef struct ocp_qp_condensing_qpoases_memory_ {
    ocp_qp_condensing_memory *condensing_memory;
    dense_qp_qpoases_memory *solver_memory;
    dense_qp_in *qpd_in;
    dense_qp_out *qpd_out;
} ocp_qp_condensing_qpoases_memory;

//
int ocp_qp_condensing_qpoases_calculate_args_size(ocp_qp_dims *dims);
//
char *ocp_qp_condensing_qpoases_assign_args(ocp_qp_dims *dims, ocp_qp_condensing_qpoases_args **args, void *mem);
//
void ocp_qp_condensing_qpoases_initialize_default_args(ocp_qp_condensing_qpoases_args *args);
//
int ocp_qp_condensing_qpoases_calculate_memory_size(ocp_qp_dims *dims, ocp_qp_condensing_qpoases_args *args);
//
char *ocp_qp_condensing_qpoases_assign_memory(ocp_qp_in *qp_in, ocp_qp_condensing_qpoases_args *args, void **mem_, void *raw_memory);
//
int ocp_qp_condensing_qpoases(ocp_qp_in *qp_in, ocp_qp_out *qp_out, void *args_, void *mem_);



#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_OCP_QP_OCP_QP_CONDENSING_QPOASES_H_
