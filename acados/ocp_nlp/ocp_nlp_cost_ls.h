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

#ifndef ACADOS_OCP_NLP_OCP_NLP_COST_LS_H_
#define ACADOS_OCP_NLP_OCP_NLP_COST_LS_H_

#ifdef __cplusplus
extern "C" {
#endif

// blasfeo
#include "blasfeo/include/blasfeo_target.h"
#include "blasfeo/include/blasfeo_common.h"
// acados
#include "acados/ocp_nlp/ocp_nlp_cost_common.h"
#include "acados/utils/types.h"
#include "acados/utils/external_function_generic.h"


typedef struct
{
	struct blasfeo_dmat Cyt;
	struct blasfeo_dmat W;
    struct blasfeo_dvec y_ref;
} ocp_nlp_cost_ls_model;

//
int ocp_nlp_cost_ls_model_calculate_size(void *config, ocp_nlp_cost_dims *dims);
//
void *ocp_nlp_cost_ls_model_assign(void *config, ocp_nlp_cost_dims *dims, void *raw_memory);
//
void ocp_nlp_cost_ls_config_initialize_default(void *config);



/* options */

typedef struct
{
    bool gauss_newton_hess; // gauss-newton hessian approximation
} ocp_nlp_cost_ls_opts;

//
int ocp_nlp_cost_ls_opts_calculate_size(void *config, ocp_nlp_cost_dims *dims);
//
void *ocp_nlp_cost_ls_opts_assign(void *config, ocp_nlp_cost_dims *dims, void *raw_memory);
//
void ocp_nlp_cost_ls_opts_initialize_default(void *config, ocp_nlp_cost_dims *dims, void *opts);



/* memory */

typedef struct
{
	struct blasfeo_dmat W_chol; // cholesky factor of weight matrix
    struct blasfeo_dvec res; // ls residual r(x)
	struct blasfeo_dvec grad; // gradient of cost function
	struct blasfeo_dvec *ux; // pointer to ux in nlp_out
	struct blasfeo_dmat *RSQrq; // pointer to RSQrq in qp_in
} ocp_nlp_cost_ls_memory;

//
int ocp_nlp_cost_ls_memory_calculate_size(void *config, ocp_nlp_cost_dims *dims, void *opts);
//
void *ocp_nlp_cost_ls_memory_assign(void *config, ocp_nlp_cost_dims *dims, void *opts, void *raw_memory);
//
struct blasfeo_dvec *ocp_nlp_cost_ls_memory_get_grad_ptr(void *memory_);
//
void ocp_nlp_cost_ls_memory_set_RSQrq_ptr(struct blasfeo_dmat *RSQrq, void *memory);
//
void ocp_nlp_cost_ls_memory_set_ux_ptr(struct blasfeo_dvec *ux, void *memory_);



/* workspace */

typedef struct
{
	struct blasfeo_dmat tmp_nv_ny;
	struct blasfeo_dvec tmp_ny;
} ocp_nlp_cost_ls_workspace;

//
int ocp_nlp_cost_ls_workspace_calculate_size(void *config, ocp_nlp_cost_dims *dims, void *opts);


/* functions */

//
void ocp_nlp_cost_ls_initialize_qp(void *config_, ocp_nlp_cost_dims *dims, void *model_, void *opts_, void *mem_, void *work_);
//
void ocp_nlp_cost_ls_update_qp_matrices(void *config_, ocp_nlp_cost_dims *dims, void *model_, void *opts_, void *memory_, void *work_);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif // ACADOS_OCP_NLP_OCP_NLP_COST_LS_H_
