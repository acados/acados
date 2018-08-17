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

#ifndef ACADOS_OCP_NLP_OCP_NLP_REG_CONV_H_
#define ACADOS_OCP_NLP_OCP_NLP_REG_CONV_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "acados/ocp_nlp/ocp_nlp_reg_common.h"

#include "blasfeo/include/blasfeo_common.h"

typedef struct {
    double *R;
    double *V;
    double *d;
    double *reg_hess;

    struct blasfeo_dmat Q_tilde;
    struct blasfeo_dmat Q_bar;
    struct blasfeo_dmat BAQ;
    struct blasfeo_dmat L;
    struct blasfeo_dmat delta_eye;
    struct blasfeo_dmat St_copy;

    struct blasfeo_dmat *original_RSQrq;

    struct blasfeo_dvec grad;
    struct blasfeo_dvec b;

} ocp_nlp_reg_conv_memory;

int ocp_nlp_reg_conv_calculate_memory_size(ocp_nlp_reg_dims *dims);

void *ocp_nlp_reg_conv_assign_memory(ocp_nlp_reg_dims *dims, void *raw_memory);

void ocp_nlp_reg_conv_config_initialize_default(ocp_nlp_reg_config *config);

#ifdef __cplusplus
}
#endif

#endif  // ACADOS_OCP_NLP_OCP_NLP_REG_CONV_H_
