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

// allocate objective and equality constraints
void allocate_ocp_qp_in_basic(int_t N, int_t *nx, int_t *nu, ocp_qp_in *const qp);
void free_ocp_qp_in_basic(ocp_qp_in *const qp);

// allocate upper/lower bounds only
void allocate_ocp_qp_in_bounds(int_t N, int_t *nb, ocp_qp_in *const qp);
void free_ocp_qp_in_bounds(ocp_qp_in *const qp);

// allocate polyhedral inequality constraints only
void allocate_ocp_qp_in_polyhedral(int_t N, int_t *nc, ocp_qp_in *const qp);
void free_ocp_qp_in_polyhedral(ocp_qp_in *const qp);

// allocate upper/lower bounds of first stage for constraint on x0
void allocate_ocp_qp_in_x0(int_t nx0, ocp_qp_in *const qp);
void free_ocp_qp_in_x0(ocp_qp_in *const qp);

// basic + bounds + polyhedral
void allocate_ocp_qp_in_full(int_t N, int_t *nx, int_t *nu, int_t *nb, int_t *nc,
    ocp_qp_in *const qp);
void free_ocp_qp_in_full(ocp_qp_in *const qp);

// basic + x0
void allocate_ocp_qp_in_unconstrained(int_t N, int_t *nx, int_t *nu, ocp_qp_in *const qp);
void free_ocp_qp_in_unconstrained(ocp_qp_in *const qp);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_UTILS_ALLOCATE_OCP_QP_IN_H_
