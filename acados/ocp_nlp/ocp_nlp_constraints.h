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
	int nh;  // number of nonlinear path constraints
	int nq;  // number of quadratic_over_nonlinear constraints TODO
    int ns;  // number of soft constraints TODO
} ocp_nlp_constraints_dims;

//
int ocp_nlp_constraints_dims_calculate_size(void *config);
//
void *ocp_nlp_constraints_dims_assign(void *config, void *raw_memory);
//
void ocp_nlp_constraints_dims_initialize(void *config, void *dims, int nx, int nu, int nbx, int nbu, int ng, int nh, int nq, int ns);



/************************************************
* config
************************************************/

typedef struct
{
	int (*dims_calculate_size) (void *config);
	void *(*dims_assign) (void *config, void *raw_memory);
	void (*dims_initialize) (void *config, void *dims, int nx, int nu, int nbx, int nbu, int ng, int nh, int nq, int ns);
	int (*model_calculate_size) (void *config, void *dims);
	void *(*model_assign) (void *config, void *dims, void *raw_memory);
	int (*opts_calculate_size) (void *config, void *dims);
	void *(*opts_assign) (void *config, void *dims, void *raw_memory);
	void (*opts_initialize_default) (void *config, void *dims, void *opts);
	void (*opts_update) (void *config, void *dims, void *opts);
	int (*memory_calculate_size) (void *config, void *dims, void *opts);
	struct blasfeo_dvec *(*memory_get_fun_ptr) (void *memory);
	struct blasfeo_dvec *(*memory_get_adj_ptr) (void *memory);
	void (*memory_set_ux_ptr) (struct blasfeo_dvec *ux, void *memory);
	void (*memory_set_lam_ptr) (struct blasfeo_dvec *lam, void *memory);
	void (*memory_set_DCt_ptr) (struct blasfeo_dmat *DCt, void *memory);
	void (*memory_set_RSQrq_ptr) (struct blasfeo_dmat *RSQrq, void *memory);
	void (*memory_set_idxb_ptr) (int *idxb, void *memory);
	void *(*memory_assign) (void *config, void *dims, void *opts, void *raw_memory);
	int (*workspace_calculate_size) (void *config, void *dims, void *opts);
	void (*initialize) (void *config, void *dims, void *model, void *opts, void *mem, void *work);
	void (*update_qp_matrices) (void *config, void *dims, void *model, void *opts, void *mem, void *work);
	void (*config_initialize_default) (void *config);
} ocp_nlp_constraints_config;

//
int ocp_nlp_constraints_config_calculate_size();
//
ocp_nlp_constraints_config *ocp_nlp_constraints_config_assign(void *raw_memory);



/************************************************
* (linear) constraints
************************************************/

/* model */

typedef struct
{
//	ocp_nlp_constraints_dims *dims;
    int *idxb;
	struct blasfeo_dvec d;
	struct blasfeo_dmat DCt;
	external_function_generic *h;
	external_function_generic *quadratic;
} ocp_nlp_constraints_model;

//
int ocp_nlp_constraints_calculate_size(void *config, void *dims);
//
void *ocp_nlp_constraints_assign(void *config, void *dims, void *raw_memory);



/* options */

typedef struct
{
    int dummy; // so cmake is happy
} ocp_nlp_constraints_opts;

//
int ocp_nlp_constraints_opts_calculate_size(void *config, void *dims);
//
void *ocp_nlp_constraints_opts_assign(void *config, void *dims, void *raw_memory);
//
void ocp_nlp_constraints_opts_initialize_default(void *config, void *dims, void *opts);
//
void ocp_nlp_constraints_opts_update(void *config, void *dims, void *opts);



/* memory */

typedef struct
{
	struct blasfeo_dvec fun;
	struct blasfeo_dvec adj;
	struct blasfeo_dvec *ux; // pointer to ux in nlp_out
	struct blasfeo_dvec *lam; // pointer to lam in nlp_out
	struct blasfeo_dmat *DCt; // pointer to DCt in qp_in
	struct blasfeo_dmat *RSQrq; // pointer to RSQrq in qp_in
	int *idxb; // pointer to idxb[ii] in qp_in
} ocp_nlp_constraints_memory;

//
int ocp_nlp_constraints_memory_calculate_size(void *config, void *dims, void *opts);
//
void *ocp_nlp_constraints_memory_assign(void *config, void *dims, void *opts, void *raw_memory);
//
struct blasfeo_dvec *ocp_nlp_constraints_memory_get_fun_ptr(void *memory_);
//
struct blasfeo_dvec *ocp_nlp_constraints_memory_get_adj_ptr(void *memory_);
//
void ocp_nlp_constraints_memory_set_ux_ptr(struct blasfeo_dvec *ux, void *memory_);
void ocp_nlp_constraints_memory_set_lam_ptr(struct blasfeo_dvec *lam, void *memory_);
//
void ocp_nlp_constraints_memory_set_DCt_ptr(struct blasfeo_dmat *DCt, void *memory);
//
void ocp_nlp_constraints_memory_set_idxb_ptr(int *idxb, void *memory_);



/* workspace */

typedef struct
{
	struct blasfeo_dvec tmp_nbgh;
//	struct blasfeo_dmat jacobian_quadratic;
} ocp_nlp_constraints_workspace;

//
int ocp_nlp_constraints_workspace_calculate_size(void *config, void *dims, void *opts);



/* functions */

//
void ocp_nlp_constraints_config_initialize_default(void *config);
//
void ocp_nlp_constraints_initialize(void *config, void *dims, void *model, void *opts, void *mem, void *work);
//
void ocp_nlp_constraints_update_qp_matrices(void *config_, void *dims, void *model_, void *opts_, void *memory_, void *work_);


#ifdef __cplusplus
} /* extern "C" */
#endif

#endif // ACADOS_OCP_NLP_OCP_NLP_CONSTRAINTS_H_
