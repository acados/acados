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

#ifndef ACADOS_OCP_NLP_OCP_NLP_CONSTRAINTS_H_
#define ACADOS_OCP_NLP_OCP_NLP_CONSTRAINTS_H_

#ifdef __cplusplus
extern "C" {
#endif

// blasfeo
#include "blasfeo/include/blasfeo_target.h"
#include "blasfeo/include/blasfeo_common.h"
// acados
#include "acados/utils/types.h"
#include "acados/utils/external_function_generic.h"
#include "acados/ocp_qp/ocp_qp_common.h"



/************************************************
* dims
************************************************/

typedef struct
{
    int nx;
    int nu;
    int nb;  // nbx + nbu
    int nbx;
    int nbu;
    int ng;  // number of general linear constraints
    int ns;  // number of soft constraints
	int nh;  // number of nonlinear path constraints
} ocp_nlp_constraints_dims;

//
int ocp_nlp_constraints_dims_calculate_size();
//
ocp_nlp_constraints_dims *ocp_nlp_constraints_dims_assign(void *raw_memory);



/************************************************
* config
************************************************/

typedef struct
{
	int (*model_calculate_size) (void *config, ocp_nlp_constraints_dims *dims);
	void *(*model_assign) (void *config, ocp_nlp_constraints_dims *dims, void *raw_memory);
	int (*opts_calculate_size) (void *config, ocp_nlp_constraints_dims *dims);
	void *(*opts_assign) (void *config, ocp_nlp_constraints_dims *dims, void *raw_memory);
	void (*opts_initialize_default) (void *config, ocp_nlp_constraints_dims *dims, void *opts);
	int (*memory_calculate_size) (void *config, ocp_nlp_constraints_dims *dims, void *opts);
	struct blasfeo_dvec *(*memory_get_fun_ptr) (void *memory);
	struct blasfeo_dvec *(*memory_get_adj_ptr) (void *memory);
	void (*memory_set_ux_ptr) (struct blasfeo_dvec *ux, void *memory);
	void (*memory_set_lam_ptr) (struct blasfeo_dvec *lam, void *memory);
	void (*memory_set_DCt_ptr) (struct blasfeo_dmat *DCt, void *memory);
	void (*memory_set_idxb_ptr) (int *idxb, void *memory);
	void *(*memory_assign) (void *config, ocp_nlp_constraints_dims *dims, void *opts, void *raw_memory);
	int (*workspace_calculate_size) (void *config, ocp_nlp_constraints_dims *dims, void *opts);
	void (*initialize_qp) (void *config, ocp_nlp_constraints_dims *dims, void *model, void *opts, void *mem, void *work);
	void (*update_qp_matrices) (void *config, ocp_nlp_constraints_dims *dims, void *model, void *opts, void *mem, void *work);
	void (*config_initialize_default) (void *config);
} ocp_nlp_constraints_config;

//
int ocp_nlp_constraints_config_calculate_size();
//
ocp_nlp_constraints_config *ocp_nlp_constraints_config_assign(void *raw_memory);



/************************************************
* linear constraints
************************************************/

/* model */

typedef struct
{
	ocp_nlp_constraints_dims *dims;
    int *idxb;
	struct blasfeo_dvec d;
	struct blasfeo_dmat DCt;
}
ocp_nlp_constraints_linear_model;

//
int ocp_nlp_constraints_linear_model_calculate_size(void *config, ocp_nlp_constraints_dims *dims);
//
void *ocp_nlp_constraints_linear_model_assign(void *config, ocp_nlp_constraints_dims *dims, void *raw_memory);



/* options */

typedef struct
{
    int dummy; // so cmake is happy
} ocp_nlp_constraints_linear_opts;

//
int ocp_nlp_constraints_linear_opts_calculate_size(void *config, ocp_nlp_constraints_dims *dims);
//
void *ocp_nlp_constraints_linear_opts_assign(void *config, ocp_nlp_constraints_dims *dims, void *raw_memory);
//
void ocp_nlp_constraints_linear_opts_initialize_default(void *config, ocp_nlp_constraints_dims *dims, void *opts);



/* memory */

typedef struct
{
	struct blasfeo_dvec fun;
	struct blasfeo_dvec adj;
	struct blasfeo_dvec *ux; // pointer to ux in nlp_out
	struct blasfeo_dvec *lam; // pointer to lam in nlp_out
	struct blasfeo_dmat *DCt; // pointer to DCt in qp_in
	int *idxb; // pointer to idxb[ii] in qp_in
} ocp_nlp_constraints_linear_memory;

//
int ocp_nlp_constraints_linear_memory_calculate_size(void *config, ocp_nlp_constraints_dims *dims, void *opts);
//
void *ocp_nlp_constraints_linear_memory_assign(void *config, ocp_nlp_constraints_dims *dims, void *opts, void *raw_memory);
//
struct blasfeo_dvec *ocp_nlp_constraints_linear_memory_get_fun_ptr(void *memory_);
//
struct blasfeo_dvec *ocp_nlp_constraints_linear_memory_get_adj_ptr(void *memory_);
//
void ocp_nlp_constraints_linear_memory_set_ux_ptr(struct blasfeo_dvec *ux, void *memory_);
void ocp_nlp_constraints_linear_memory_set_lam_ptr(struct blasfeo_dvec *lam, void *memory_);
//
void ocp_nlp_constraints_linear_memory_set_DCt_ptr(struct blasfeo_dmat *DCt, void *memory);
//
void ocp_nlp_constraints_linear_memory_set_idxb_ptr(int *idxb, void *memory_);



/* workspace */

typedef struct
{
	struct blasfeo_dvec tmp_nbg;
} ocp_nlp_constraints_linear_workspace;

//
int ocp_nlp_constraints_linear_workspace_calculate_size(void *config, ocp_nlp_constraints_dims *dims, void *opts);



/* functions */

//
void ocp_nlp_constraints_linear_config_initialize_default(void *config);
//
void ocp_nlp_constraints_linear_initialize_qp(void *config, ocp_nlp_constraints_dims *dims, void *model, void *opts, void *mem, void *work);
//
void ocp_nlp_constraints_linear_update_qp_matrices(void *config_, ocp_nlp_constraints_dims *dims, void *model_, void *opts_, void *memory_, void *work_);


#ifdef __cplusplus
} /* extern "C" */
#endif

#endif // ACADOS_OCP_NLP_OCP_NLP_CONSTRAINTS_H_
