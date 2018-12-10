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

#ifndef ACADOS_OCP_NLP_OCP_NLP_COST_COMMON_H_
#define ACADOS_OCP_NLP_OCP_NLP_COST_COMMON_H_

#ifdef __cplusplus
extern "C" {
#endif

// acados
#include "acados/utils/external_function_generic.h"
#include "acados/utils/types.h"



/************************************************
 * config
 ************************************************/

typedef struct
{
    int (*dims_calculate_size)(void *config);
    void *(*dims_assign)(void *config, void *raw_memory);
    void (*dims_initialize)(void *config, void *dims, int nx, int nu, int ny, int ns);
    void (*dims_set)(void *config_, void *dims_, const char *field, int *value);
    int (*model_calculate_size)(void *config, void *dims);
    void *(*model_assign)(void *config, void *dims, void *raw_memory);
    int (*model_set)(void *config_, void *dims_, void *model_, const char *field, void *value_);
    int (*opts_calculate_size)(void *config, void *dims);
    void *(*opts_assign)(void *config, void *dims, void *raw_memory);
    void (*opts_initialize_default)(void *config, void *dims, void *opts);
    void (*opts_update)(void *config, void *dims, void *opts);
    int (*memory_calculate_size)(void *config, void *dims, void *opts);
    struct blasfeo_dvec *(*memory_get_grad_ptr)(void *memory);
    void (*memory_set_ux_ptr)(struct blasfeo_dvec *ux, void *memory);
    void (*memory_set_RSQrq_ptr)(struct blasfeo_dmat *RSQrq, void *memory);
    void (*memory_set_Z_ptr)(struct blasfeo_dvec *Z, void *memory);
    void *(*memory_assign)(void *config, void *dims, void *opts, void *raw_memory);
    int (*workspace_calculate_size)(void *config, void *dims, void *opts);
    void (*initialize)(void *config_, void *dims, void *model_, void *opts_, void *mem_,
                       void *work_);
    void (*update_qp_matrices)(void *config_, void *dims, void *model_, void *opts_, void *mem_,
                               void *work_);
    void (*config_initialize_default)(void *config);
} ocp_nlp_cost_config;

//
int ocp_nlp_cost_config_calculate_size();
//
ocp_nlp_cost_config *ocp_nlp_cost_config_assign(void *raw_memory);



#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_OCP_NLP_OCP_NLP_COST_COMMON_H_
