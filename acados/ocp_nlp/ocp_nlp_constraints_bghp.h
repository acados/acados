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

#ifndef ACADOS_OCP_NLP_OCP_NLP_CONSTRAINTS_BGHP_H_
#define ACADOS_OCP_NLP_OCP_NLP_CONSTRAINTS_BGHP_H_

#ifdef __cplusplus
extern "C" {
#endif

// acados
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/ocp_nlp/ocp_nlp_constraints_common.h"
#include "acados/utils/external_function_generic.h"
#include "acados/utils/types.h"



/* dims */

typedef struct
{
    int nx;
    int nu;
    int nb;  // nbx + nbu
    int nbu;
    int nbx;
    int ng;  // number of general linear constraints
    int nh;  // number of nonlinear path constraints
    int ns;  // nsbu + nsbx + nsg + nsh
    int nsbu;  // number of softed input bounds
    int nsbx;  // number of softed state bounds
    int nsg;  // number of softed general linear constraints
    int nsh;  // number of softed nonlinear constraints
    int np;  // dimension of nonlinear function in quadratic_over_nonlinear constraint
} ocp_nlp_constraints_bghp_dims;

//
int ocp_nlp_constraints_bghp_dims_calculate_size(void *config);
//
void *ocp_nlp_constraints_bghp_dims_assign(void *config, void *raw_memory);
//
void ocp_nlp_constraints_bghp_dims_initialize(void *config, void *dims, int nx, int nu, int nbx,
                                         int nbu, int ng, int nh, int nq, int ns);
//
void ocp_nlp_constraints_bghp_dims_get(void *config_, void *dims_, const char *field, int* value);


/* model */

typedef struct
{
    //  ocp_nlp_constraints_bghp_dims *dims;
    int *idxb;
    int *idxs;
    struct blasfeo_dvec d;
    struct blasfeo_dmat DCt;
    external_function_generic *nl_constr_h_fun_jac;
    external_function_generic *p;
} ocp_nlp_constraints_bghp_model;

//
int ocp_nlp_constraints_bghp_calculate_size(void *config, void *dims);
//
void *ocp_nlp_constraints_bghp_assign(void *config, void *dims, void *raw_memory);
//
int ocp_nlp_constraints_bghp_model_set(void *config_, void *dims_,
                         void *model_, const char *field, void *value);

/* options */

typedef struct
{
    int compute_adj;
} ocp_nlp_constraints_bghp_opts;

//
int ocp_nlp_constraints_bghp_opts_calculate_size(void *config, void *dims);
//
void *ocp_nlp_constraints_bghp_opts_assign(void *config, void *dims, void *raw_memory);
//
void ocp_nlp_constraints_bghp_opts_initialize_default(void *config, void *dims, void *opts);
//
void ocp_nlp_constraints_bghp_opts_update(void *config, void *dims, void *opts);
//
void ocp_nlp_constraints_bghp_opts_set(void *config, void *dims, void *opts, enum acados_opts name,
    void *ptr_value);

/* memory */

typedef struct
{
    struct blasfeo_dvec fun;
    struct blasfeo_dvec adj;
    struct blasfeo_dvec *ux;     // pointer to ux in nlp_out
    struct blasfeo_dvec *lam;    // pointer to lam in nlp_out
    struct blasfeo_dmat *DCt;    // pointer to DCt in qp_in
    struct blasfeo_dmat *RSQrq;  // pointer to RSQrq in qp_in
    int *idxb;                   // pointer to idxb[ii] in qp_in
    int *idxs;                   // pointer to idxs[ii] in qp_in
} ocp_nlp_constraints_bghp_memory;

//
int ocp_nlp_constraints_bghp_memory_calculate_size(void *config, void *dims, void *opts);
//
void *ocp_nlp_constraints_bghp_memory_assign(void *config, void *dims, void *opts,
    void *raw_memory);
//
struct blasfeo_dvec *ocp_nlp_constraints_bghp_memory_get_fun_ptr(void *memory_);
//
struct blasfeo_dvec *ocp_nlp_constraints_bghp_memory_get_adj_ptr(void *memory_);
//
void ocp_nlp_constraints_bghp_memory_set_ux_ptr(struct blasfeo_dvec *ux, void *memory_);
void ocp_nlp_constraints_bghp_memory_set_lam_ptr(struct blasfeo_dvec *lam, void *memory_);
//
void ocp_nlp_constraints_bghp_memory_set_DCt_ptr(struct blasfeo_dmat *DCt, void *memory);
//
void ocp_nlp_constraints_bghp_memory_set_idxb_ptr(int *idxb, void *memory_);
//
void ocp_nlp_constraints_bghp_memory_set_idxs_ptr(int *idxs, void *memory_);

/* workspace */

typedef struct
{
    struct blasfeo_dvec tmp_ni;
    struct blasfeo_dmat jacobian_quadratic;
} ocp_nlp_constraints_bghp_workspace;

//
int ocp_nlp_constraints_bghp_workspace_calculate_size(void *config, void *dims, void *opts);

/* functions */

//
void ocp_nlp_constraints_bghp_config_initialize_default(void *config);
//
void ocp_nlp_constraints_bghp_initialize(void *config, void *dims, void *model, void *opts,
                                    void *mem, void *work);
//
void ocp_nlp_constraints_bghp_update_qp_matrices(void *config_, void *dims, void *model_,
                                            void *opts_, void *memory_, void *work_);



#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_OCP_NLP_OCP_NLP_CONSTRAINTS_BGHP_H_
