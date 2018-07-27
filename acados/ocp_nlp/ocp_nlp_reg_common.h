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

#ifndef ACADOS_OCP_NLP_OCP_NLP_REG_COMMON_H_
#define ACADOS_OCP_NLP_OCP_NLP_REG_COMMON_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "acados/ocp_qp/ocp_qp_common.h"

typedef ocp_qp_dims ocp_nlp_reg_dims;

typedef ocp_qp_in ocp_nlp_reg_in;

typedef ocp_qp_out ocp_nlp_reg_out;

typedef struct {
    double delta;
    double gamma;
} ocp_nlp_reg_opts;

typedef struct {
    int (*opts_calculate_size)(void);
    void *(*opts_assign)(void *raw_memory);

    int (*memory_calculate_size)(ocp_nlp_reg_dims *dims);
    void *(*memory_assign)(ocp_nlp_reg_dims *dims, void *raw_memory);

    void (*evaluate)(void *config, ocp_nlp_reg_dims *dims, ocp_nlp_reg_in *in, ocp_nlp_reg_out *out,
                    ocp_nlp_reg_opts *opts, void *mem_);
} ocp_nlp_reg_config;

int ocp_nlp_reg_opts_calculate_size(void);

void *ocp_nlp_reg_opts_assign(void *raw_memory);

int ocp_nlp_reg_config_calculate_size(void);

void *ocp_nlp_reg_config_assign(void *raw_memory);

#ifdef __cplusplus
}
#endif

#endif  // ACADOS_OCP_NLP_OCP_NLP_REG_COMMON_H_
