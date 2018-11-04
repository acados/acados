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

#ifndef ACADOS_OCP_NLP_OCP_NLP_REG_MIRROR_H_
#define ACADOS_OCP_NLP_OCP_NLP_REG_MIRROR_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "acados/ocp_nlp/ocp_nlp_reg_common.h"

typedef struct {
    double *reg_hess;
    double *V;
    double *d;
} ocp_nlp_reg_mirror_memory;

int ocp_nlp_reg_mirror_memory_calculate_size(ocp_nlp_reg_dims *dims);

void *ocp_nlp_reg_mirror_memory_assign(ocp_nlp_reg_dims *dims, void *raw_memory);

void ocp_nlp_reg_mirror_config_initialize_default(ocp_nlp_reg_config *config);

#ifdef __cplusplus
}
#endif

#endif  // ACADOS_OCP_NLP_OCP_NLP_REG_MIRROR_H_
