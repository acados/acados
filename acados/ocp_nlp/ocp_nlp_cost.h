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
	struct blasfeo_dvec *(*memory_get_grad_ptr) (void *memory);
	void (*memory_set_ux_ptr) (struct blasfeo_dvec *ux, void *memory);
	void (*memory_set_RSQrq_ptr) (struct blasfeo_dmat *RSQrq, void *memory);
	void *(*memory_assign) (void *config, ocp_nlp_cost_dims *dims, void *opts, void *raw_memory);
	int (*workspace_calculate_size) (void *config, ocp_nlp_cost_dims *dims, void *opts);
	void (*initialize_qp) (void *config_, ocp_nlp_cost_dims *dims, void *model_, void *opts_, void *mem_, void *work_);
	void (*update_qp_matrices) (void *config_, ocp_nlp_cost_dims *dims, void *model_, void *opts_, void *mem_, void *work_);
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
	struct blasfeo_dmat W_chol; // cholesky factor of weight matrix
	struct blasfeo_dmat Jt; // jacobian of nls fun
    struct blasfeo_dvec res; // nls residual r(x)
	struct blasfeo_dvec grad; // gradient of cost function
	struct blasfeo_dvec *ux; // pointer to ux in nlp_out
	struct blasfeo_dmat *RSQrq; // pointer to RSQrq in qp_in
} ocp_nlp_cost_nls_memory;

//
int ocp_nlp_cost_nls_memory_calculate_size(void *config, ocp_nlp_cost_dims *dims, void *opts);
//
void *ocp_nlp_cost_nls_memory_assign(void *config, ocp_nlp_cost_dims *dims, void *opts, void *raw_memory);
//
struct blasfeo_dvec *ocp_nlp_cost_nls_memory_get_grad_ptr(void *memory_);
//
void ocp_nlp_cost_nls_memory_set_RSQrq_ptr(struct blasfeo_dmat *RSQrq, void *memory);
//
void ocp_nlp_cost_nls_memory_set_ux_ptr(struct blasfeo_dvec *ux, void *memory_);



/* workspace */

typedef struct
{
	struct blasfeo_dmat tmp_nv_ny;
	struct blasfeo_dvec tmp_ny;
	double *nls_jac_in;
	double *nls_jac_out;
} ocp_nlp_cost_nls_workspace;

//
int ocp_nlp_cost_nls_workspace_calculate_size(void *config, ocp_nlp_cost_dims *dims, void *opts);


/* functions */

//
void ocp_nlp_cost_nls_initialize_qp(void *config_, ocp_nlp_cost_dims *dims, void *model_, void *opts_, void *mem_, void *work_);
//
void ocp_nlp_cost_nls_update_qp_matrices(void *config_, ocp_nlp_cost_dims *dims, void *model_, void *opts_, void *memory_, void *work_);



#ifdef __cplusplus
} /* extern "C" */
#endif

#endif // ACADOS_OCP_NLP_OCP_NLP_COST_H_
