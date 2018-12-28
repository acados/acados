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

#ifndef ACADOS_OCP_NLP_OCP_NLP_CONSTRAINTS_COMMON_H_
#define ACADOS_OCP_NLP_OCP_NLP_CONSTRAINTS_COMMON_H_

#ifdef __cplusplus
extern "C" {
#endif

// acados
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/utils/external_function_generic.h"
#include "acados/utils/types.h"



/************************************************
 * config
 ************************************************/

typedef struct
{
    int (*dims_calculate_size)(void *config);
    void *(*dims_assign)(void *config, void *raw_memory);
    void (*dims_initialize)(void *config, void *dims, int nx, int nu, int nbx, int nbu, int ng,
                            int nh, int nq, int ns);
    int (*model_calculate_size)(void *config, void *dims);
    void *(*model_assign)(void *config, void *dims, void *raw_memory);
    int (*model_set)(void *config_, void *dims_, void *model_, const char *field, void *value);
    int (*opts_calculate_size)(void *config, void *dims);
    void *(*opts_assign)(void *config, void *dims, void *raw_memory);
    void (*opts_initialize_default)(void *config, void *dims, void *opts);
    void (*opts_update)(void *config, void *dims, void *opts);
    void (*opts_set)(void *config_, void *dims_, void *opts_, enum acados_opts name,
        void *ptr_value);
    int (*memory_calculate_size)(void *config, void *dims, void *opts);
    struct blasfeo_dvec *(*memory_get_fun_ptr)(void *memory);
    struct blasfeo_dvec *(*memory_get_adj_ptr)(void *memory);
    void (*memory_set_ux_ptr)(struct blasfeo_dvec *ux, void *memory);
    void (*memory_set_lam_ptr)(struct blasfeo_dvec *lam, void *memory);
    void (*memory_set_DCt_ptr)(struct blasfeo_dmat *DCt, void *memory);
    void (*memory_set_RSQrq_ptr)(struct blasfeo_dmat *RSQrq, void *memory);
    void (*memory_set_idxb_ptr)(int *idxb, void *memory);
    void (*memory_set_idxs_ptr)(int *idxs, void *memory);
    void *(*memory_assign)(void *config, void *dims, void *opts, void *raw_memory);
    int (*workspace_calculate_size)(void *config, void *dims, void *opts);
    void (*initialize)(void *config, void *dims, void *model, void *opts, void *mem, void *work);
    void (*update_qp_matrices)(void *config, void *dims, void *model, void *opts, void *mem,
                               void *work);
    void (*config_initialize_default)(void *config);
    // dimension setters
    void (*dims_set)(void *config_, void *dims_, const char *field, const int *value);
    void (*get_dims)(void *config_, void *dims_, const char *field, int* value);
} ocp_nlp_constraints_config;

//
int ocp_nlp_constraints_config_calculate_size();
//
ocp_nlp_constraints_config *ocp_nlp_constraints_config_assign(void *raw_memory);



#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_OCP_NLP_OCP_NLP_CONSTRAINTS_COMMON_H_
