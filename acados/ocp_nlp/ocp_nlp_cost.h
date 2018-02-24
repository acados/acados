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

#ifndef ACADOS_OCP_NLP_OCP_NLP_COST_H_
#define ACADOS_OCP_NLP_OCP_NLP_COST_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/utils/types.h"
#include "acados/utils/external_function_generic.h"


/************************************************
* dims
************************************************/

typedef struct
{
	int nx; // number of states
	int nu; // number of inputs
	int ny; // number of outputs
} ocp_nlp_cost_dims;

//
int ocp_nlp_cost_dims_calculate_size();
//
ocp_nlp_cost_dims *ocp_nlp_cost_dims_assign(void *raw_memory);



/************************************************
* config
************************************************/

typedef struct
{
	int (*model_calculate_size) (void *config, ocp_nlp_cost_dims *dims);
	void *(*model_assign) (void *config, ocp_nlp_cost_dims *dims, void *raw_memory);
	void (*config_initialize_default) (void *config);
} ocp_nlp_cost_config;

//
int ocp_nlp_cost_config_calculate_size();
//
ocp_nlp_cost_config *ocp_nlp_cost_config_assign(void *raw_memory);


/************************************************
* least squares
************************************************/

/* least squares */

typedef struct
{
	ocp_nlp_cost_dims *dims;
	external_function_generic *nls_jac; // array of N+1 pointers; evaluation and jacobian of ls residuals
	struct blasfeo_dmat Cyt;
	struct blasfeo_dmat W;
    struct blasfeo_dvec y_ref;
	int nls_mask; // nonlinear least squares mask
} ocp_nlp_cost_ls_model;

//
int ocp_nlp_cost_ls_model_calculate_size(void *config, ocp_nlp_cost_dims *dims);
//
void *ocp_nlp_cost_ls_model_assign(void *config, ocp_nlp_cost_dims *dims, void *raw_memory);
//
void ocp_nlp_cost_ls_config_initialize_default(void *config_);



#endif // ACADOS_OCP_NLP_OCP_NLP_COST_H_
