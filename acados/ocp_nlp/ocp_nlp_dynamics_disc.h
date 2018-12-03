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

#ifndef ACADOS_OCP_NLP_OCP_NLP_DYNAMICS_DISC_H_
#define ACADOS_OCP_NLP_OCP_NLP_DYNAMICS_DISC_H_

#ifdef __cplusplus
extern "C" {
#endif

// blasfeo
#include "blasfeo/include/blasfeo_common.h"

// acados
#include "acados/ocp_nlp/ocp_nlp_dynamics_common.h"
#include "acados/utils/external_function_generic.h"
#include "acados/utils/types.h"

/************************************************
 * dims
 ************************************************/

typedef struct
{
    int nx;   // number of states at the current stage
    int nu;   // number of inputs at the current stage
    int nx1;  // number of states at the next stage
    int nu1;  // number of inputes at the next stage
} ocp_nlp_dynamics_disc_dims;

//
int ocp_nlp_dynamics_disc_dims_calculate_size(void *config);
//
void *ocp_nlp_dynamics_disc_dims_assign(void *config, void *raw_memory);
//
void ocp_nlp_dynamics_disc_dims_initialize(void *config, void *dims, int nx, int nu, int nx1,
                                           int nu1, int nz);

//
void ocp_nlp_dynamics_disc_dims_set(void *config_, void *dims_, const char *dim, int* value);


/************************************************
 * options
 ************************************************/

typedef struct
{
    int compute_adj;
} ocp_nlp_dynamics_disc_opts;

//
int ocp_nlp_dynamics_disc_opts_calculate_size(void *config, void *dims);
//
void *ocp_nlp_dynamics_disc_opts_assign(void *config, void *dims, void *raw_memory);
//
void ocp_nlp_dynamics_disc_opts_initialize_default(void *config, void *dims, void *opts);
//
void ocp_nlp_dynamics_disc_opts_update(void *config, void *dims, void *opts);
//
int ocp_nlp_dynamics_disc_precompute(void *config_, void *dims, void *model_, void *opts_,
                                        void *mem_, void *work_);


/************************************************
 * memory
 ************************************************/

typedef struct
{
    struct blasfeo_dvec fun;
    struct blasfeo_dvec adj;
    struct blasfeo_dvec *ux;     // pointer to ux in nlp_out at current stage
    struct blasfeo_dvec *ux1;    // pointer to ux in nlp_out at next stage
    struct blasfeo_dvec *pi;     // pointer to pi in nlp_out at current stage
    struct blasfeo_dmat *BAbt;   // pointer to BAbt in qp_in
    struct blasfeo_dmat *RSQrq;  // pointer to RSQrq in qp_in
} ocp_nlp_dynamics_disc_memory;

//
int ocp_nlp_dynamics_disc_memory_calculate_size(void *config, void *dims, void *opts);
//
void *ocp_nlp_dynamics_disc_memory_assign(void *config, void *dims, void *opts, void *raw_memory);
//
struct blasfeo_dvec *ocp_nlp_dynamics_disc_memory_get_fun_ptr(void *memory);
//
struct blasfeo_dvec *ocp_nlp_dynamics_disc_memory_get_adj_ptr(void *memory);
//
void ocp_nlp_dynamics_disc_memory_set_ux_ptr(struct blasfeo_dvec *ux, void *memory);
//
void ocp_nlp_dynamics_disc_memory_set_ux1_ptr(struct blasfeo_dvec *ux1, void *memory);
//
void ocp_nlp_dynamics_disc_memory_set_pi_ptr(struct blasfeo_dvec *pi, void *memory);
//
void ocp_nlp_dynamics_disc_memory_set_BAbt_ptr(struct blasfeo_dmat *BAbt, void *memory);



/************************************************
 * workspace
 ************************************************/

typedef struct
{
    double *discrete_model_in;
    double *discrete_model_out;
} ocp_nlp_dynamics_disc_workspace;

int ocp_nlp_dynamics_disc_workspace_calculate_size(void *config, void *dims, void *opts);



/************************************************
 * model
 ************************************************/

typedef struct
{
    external_function_generic *discrete_model;
} ocp_nlp_dynamics_disc_model;

//
int ocp_nlp_dynamics_disc_model_calculate_size(void *config, void *dims);
//
void *ocp_nlp_dynamics_disc_model_assign(void *config, void *dims, void *raw_memory);
//
void ocp_nlp_dynamics_disc_model_set_T(double T, void *model);



/************************************************
 * functions
 ************************************************/

//
void ocp_nlp_dynamics_disc_config_initialize_default(void *config);
//
void ocp_nlp_dynamics_disc_initialize(void *config_, void *dims, void *model_, void *opts,
                                      void *mem, void *work_);
//
void ocp_nlp_dynamics_disc_update_qp_matrices(void *config_, void *dims, void *model_, void *opts,
                                              void *mem, void *work_);



#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_OCP_NLP_OCP_NLP_DYNAMICS_DISC_H_
