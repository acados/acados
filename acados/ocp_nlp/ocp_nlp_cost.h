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

// blasfeo
#include "blasfeo/include/blasfeo_target.h"
#include "blasfeo/include/blasfeo_common.h"
// acados
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
	int (*opts_calculate_size) (void *config, ocp_nlp_cost_dims *dims);
	void *(*opts_assign) (void *config, ocp_nlp_cost_dims *dims, void *raw_memory);
	void (*opts_initialize_default) (void *config, ocp_nlp_cost_dims *dims, void *opts);
	int (*memory_calculate_size) (void *config, ocp_nlp_cost_dims *dims, void *opts);
	void *(*memory_assign) (void *config, ocp_nlp_cost_dims *dims, void *opts, void *raw_memory);
	int (*workspace_calculate_size) (void *config, ocp_nlp_cost_dims *dims, void *opts);
	void (*config_initialize_default) (void *config);
} ocp_nlp_cost_config;

//
int ocp_nlp_cost_config_calculate_size();
//
ocp_nlp_cost_config *ocp_nlp_cost_config_assign(void *raw_memory);



/************************************************
* linear least squares
************************************************/

/* model */

typedef struct
{
	ocp_nlp_cost_dims *dims;
	external_function_generic *nls_jac; // evaluation and jacobian of ls residuals
	struct blasfeo_dmat Cyt;
	struct blasfeo_dmat W;
    struct blasfeo_dvec y_ref;
	int nls_mask; // nonlinear least squares mask TODO lin and nonlin models instead
} ocp_nlp_cost_ls_model;

//
int ocp_nlp_cost_ls_model_calculate_size(void *config, ocp_nlp_cost_dims *dims);
//
void *ocp_nlp_cost_ls_model_assign(void *config, ocp_nlp_cost_dims *dims, void *raw_memory);
//
void ocp_nlp_cost_ls_config_initialize_default(void *config);


/************************************************
* nonlinear least squares
************************************************/

/* model */

typedef struct
{
	ocp_nlp_cost_dims *dims;
	external_function_generic *nls_jac; // evaluation and jacobian of ls residuals
	struct blasfeo_dmat Cyt;
	struct blasfeo_dmat W;
    struct blasfeo_dvec y_ref;
	int nls_mask; // nonlinear least squares mask TODO lin and nonlin models instead
} ocp_nlp_cost_nls_model;

//
int ocp_nlp_cost_nls_model_calculate_size(void *config, ocp_nlp_cost_dims *dims);
//
void *ocp_nlp_cost_nls_model_assign(void *config, ocp_nlp_cost_dims *dims, void *raw_memory);
//
void ocp_nlp_cost_nls_config_initialize_default(void *config);



/* options */

typedef struct
{
    int dummy; // XXX to make cmake happy
} ocp_nlp_cost_nls_opts;

//
int ocp_nlp_cost_nls_opts_calculate_size(void *config, ocp_nlp_cost_dims *dims);
//
void *ocp_nlp_cost_nls_opts_assign(void *config, ocp_nlp_cost_dims *dims, void *raw_memory);
//
void ocp_nlp_cost_nls_opts_initialize_default(void *config, ocp_nlp_cost_dims *dims, void *opts);



/* memory */

typedef struct
{
    int dummy; // XXX to make cmake happy
} ocp_nlp_cost_nls_memory;

//
int ocp_nlp_cost_nls_memory_calculate_size(void *config, ocp_nlp_cost_dims *dims, void *opts);
//
void *ocp_nlp_cost_nls_memory_assign(void *config, ocp_nlp_cost_dims *dims, void *opts, void *raw_memory);



/* workspace */

typedef struct
{
    int dummy; // XXX to make cmake happy
} ocp_nlp_cost_nls_workspace;

//
int ocp_nlp_cost_nls_workspace_calculate_size(void *config, ocp_nlp_cost_dims *dims, void *opts);


#ifdef __cplusplus
} /* extern "C" */
#endif

#endif // ACADOS_OCP_NLP_OCP_NLP_COST_H_
