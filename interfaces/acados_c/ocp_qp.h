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

#ifndef ACADOS_C_OCP_QP_H_
#define ACADOS_C_OCP_QP_H_

#ifdef __cplusplus
extern "C" {
#endif

// acados
#include <acados/ocp_qp/ocp_qp_common.h>
#include <acados/utils/types.h>

// NOTE(dimitris): contains both ocp_qp solvers and condensing with dense
// solvers
typedef enum {
    HPIPM,
    CONDENSING_HPIPM,
    CONDENSING_QPOASES,
    CONDENSING_QORE
} qp_solver_t;

//
int set_qp_solver_fun_ptrs(qp_solver_t qp_solver_name, void *qp_solver);
//
void set_xcond_qp_solver_fun_ptrs(qp_solver_t qp_solver_name, ocp_qp_xcond_solver *qp_solver);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_C_OCP_QP_H_
