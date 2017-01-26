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

#ifndef ACADOS_UTILS_ALLOCATE_OCP_QP_IN_H_
#define ACADOS_UTILS_ALLOCATE_OCP_QP_IN_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "acados/utils/types.h"
#include "acados/ocp_qp/ocp_qp_common.h"

// allocate a full QP
void allocate_QP_full(int_t N, int_t *nx, int_t *nu, int_t *nb, int_t *nc, ocp_qp_in *const qp);
void free_QP_full(ocp_qp_in *const qp);

// allocate only equality constraints and upper/lower bound on x0
void allocate_QP_unconstrained(int_t N, int_t *nx, int_t *nu, ocp_qp_in *const qp);
void free_QP_unconstrained(ocp_qp_in *const qp);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_UTILS_ALLOCATE_OCP_QP_IN_H_
