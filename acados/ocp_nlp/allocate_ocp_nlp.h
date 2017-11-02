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

#ifndef ACADOS_OCP_NLP_ALLOCATE_OCP_NLP_H_
#define ACADOS_OCP_NLP_ALLOCATE_OCP_NLP_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "acados/ocp_nlp/ocp_nlp_common.h"
#include "acados/ocp_nlp/ocp_nlp_sm_gn.h"
#include "acados/utils/types.h"

void allocate_ocp_nlp_in(int_t N, int_t *nx, int_t *nu, int_t *nb, int_t *ng,
                         ocp_nlp_in *const nlp);
void free_ocp_nlp_in(ocp_nlp_in *const nlp);

void allocate_ocp_nlp_out(ocp_nlp_in *const in, ocp_nlp_out *out);

void free_ocp_nlp_out(int_t N, ocp_nlp_out *out);

void allocate_ls_cost(int_t N, int_t *nx, int_t *nu, int_t *ny, ocp_nlp_ls_cost *ls_cost);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_OCP_NLP_ALLOCATE_OCP_NLP_H_
