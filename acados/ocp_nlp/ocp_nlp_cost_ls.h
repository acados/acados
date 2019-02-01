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
#include "blasfeo/include/blasfeo_common.h"

// acados
#include "acados/ocp_nlp/ocp_nlp_cost_common.h"
#include "acados/utils/external_function_generic.h"
#include "acados/utils/types.h"

/************************************************
 * dims
 ************************************************/

typedef struct
{
    int nx;  // number of states
    int nu;  // number of inputs
    int ny;  // number of outputs
    int ns;  // number of slacks
} ocp_nlp_cost_ls_dims;

//
int ocp_nlp_cost_ls_dims_calculate_size(void *config);
//
void *ocp_nlp_cost_ls_dims_assign(void *config, void *raw_memory);
//
void ocp_nlp_cost_ls_dims_initialize(void *config, void *dims, int nx, int nu, int ny, int ns);
//
void ocp_nlp_cost_ls_dims_set(void *config_, void *dims_, const char *field, int* value);

/************************************************
 * model
 ************************************************/

typedef struct
{
    // slack penalty has the form z^T * s + .5 * s^T * Z * s
    struct blasfeo_dmat Cyt;            // output matrix: Cy * [x, u] = y; in transposed form
    struct blasfeo_dmat W;              // ls norm corresponding to this matrix
    struct blasfeo_dvec y_ref;          // yref
    struct blasfeo_dvec Z;              // diagonal Hessian of slacks as vector (lower and upper)
    struct blasfeo_dvec z;              // gradient of slacks as vector (lower and upper)
	double scaling;
} ocp_nlp_cost_ls_model;

//
int ocp_nlp_cost_ls_model_calculate_size(void *config, void *dims);
//
void *ocp_nlp_cost_ls_model_assign(void *config, void *dims, void *raw_memory);
//
int ocp_nlp_cost_ls_model_set(void *config_, void *dims_, void *model_,
                              const char *field, void *value_);

/************************************************
 * options
 ************************************************/

typedef struct
{
    int dummy; // struct can't be void
} ocp_nlp_cost_ls_opts;

//
int ocp_nlp_cost_ls_opts_calculate_size(void *config, void *dims);
//
void *ocp_nlp_cost_ls_opts_assign(void *config, void *dims, void *raw_memory);
//
void ocp_nlp_cost_ls_opts_initialize_default(void *config, void *dims, void *opts);
//
void ocp_nlp_cost_ls_opts_update(void *config, void *dims, void *opts);

/************************************************
 * memory
 ************************************************/

typedef struct
{
    struct blasfeo_dmat hess;    // hessian of cost function
    struct blasfeo_dmat W_chol;  // cholesky factor of weight matrix
    struct blasfeo_dvec res;     // ls residual r(x)
    struct blasfeo_dvec grad;    // gradient of cost function
    struct blasfeo_dvec *ux;     // pointer to ux in nlp_out
    struct blasfeo_dmat *RSQrq;  // pointer to RSQrq in qp_in
    struct blasfeo_dvec *Z;      // pointer to Z in qp_in
} ocp_nlp_cost_ls_memory;

//
int ocp_nlp_cost_ls_memory_calculate_size(void *config, void *dims, void *opts);
//
void *ocp_nlp_cost_ls_memory_assign(void *config, void *dims, void *opts, void *raw_memory);
//
struct blasfeo_dvec *ocp_nlp_cost_ls_memory_get_grad_ptr(void *memory_);
//
void ocp_nlp_cost_ls_memory_set_RSQrq_ptr(struct blasfeo_dmat *RSQrq, void *memory);
//
void ocp_nlp_cost_ls_memory_set_Z_ptr(struct blasfeo_dvec *Z, void *memory);
//
void ocp_nlp_cost_ls_memory_set_ux_ptr(struct blasfeo_dvec *ux, void *memory_);

/************************************************
 * workspace
 ************************************************/

typedef struct
{
    struct blasfeo_dmat tmp_nv_ny;
    struct blasfeo_dvec tmp_ny;
} ocp_nlp_cost_ls_workspace;

//
int ocp_nlp_cost_ls_workspace_calculate_size(void *config, void *dims, void *opts);

/************************************************
 * functions
 ************************************************/

//
void ocp_nlp_cost_ls_config_initialize_default(void *config);
//
void ocp_nlp_cost_ls_initialize(void *config_, void *dims, void *model_, void *opts_, void *mem_,
                                void *work_);
//
void ocp_nlp_cost_ls_update_qp_matrices(void *config_, void *dims, void *model_, void *opts_,
                                        void *memory_, void *work_);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_OCP_NLP_OCP_NLP_COST_LS_H_
