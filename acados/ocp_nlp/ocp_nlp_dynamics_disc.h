/*
 * Copyright (c) The acados authors.
 *
 * This file is part of acados.
 *
 * Licensed under the 2-Clause BSD License.
 */


/// \addtogroup ocp_nlp
/// @{
/// \addtogroup ocp_nlp_dynamics
/// @{

#ifndef ACADOS_OCP_NLP_OCP_NLP_DYNAMICS_DISC_H_
#define ACADOS_OCP_NLP_OCP_NLP_DYNAMICS_DISC_H_

#ifdef __cplusplus
extern "C" {
#endif

// blasfeo
#include "blasfeo_common.h"

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
    int np;   // number of parameters
    int np_global;   // number of global parameters

} ocp_nlp_dynamics_disc_dims;

//
acados_size_t ocp_nlp_dynamics_disc_dims_calculate_size(void *config);
//
void *ocp_nlp_dynamics_disc_dims_assign(void *config, void *raw_memory);
//
void ocp_nlp_dynamics_disc_dims_set(void *config_, void *dims_, const char *dim, int* value);


/************************************************
 * options
 ************************************************/

typedef struct
{
    int compute_adj;
    int compute_hess;
    int cost_computation;
    int with_solution_sens_wrt_params_forw;
    int with_solution_sens_wrt_params_adj;
} ocp_nlp_dynamics_disc_opts;

//
acados_size_t ocp_nlp_dynamics_disc_opts_calculate_size(void *config, void *dims);
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
    struct blasfeo_dmat *dyn_jac_p_global;  // pointer to jacobian of the dynamics wrt the parameters
    struct blasfeo_dmat *jac_lag_stat_p_global;    // pointer to jacobian of stationarity condition wrt parameters
    struct blasfeo_dvec fun;
    struct blasfeo_dvec adj;
    struct blasfeo_dvec *ux;     // pointer to ux in nlp_out at current stage
    struct blasfeo_dvec *ux1;    // pointer to ux in nlp_out at next stage
    struct blasfeo_dvec *pi;     // pointer to pi in nlp_out at current stage
    struct blasfeo_dmat *BAbt;   // pointer to BAbt in qp_in
    struct blasfeo_dmat *RSQrq;  // pointer to RSQrq in qp_in

    struct blasfeo_dvec *seed_ux;
    struct blasfeo_dvec *seed_pi;
    struct blasfeo_dvec *adj_lag_p_global;
} ocp_nlp_dynamics_disc_memory;

//
acados_size_t ocp_nlp_dynamics_disc_memory_calculate_size(void *config, void *dims, void *opts);
//
void *ocp_nlp_dynamics_disc_memory_assign(void *config, void *dims, void *opts, void *raw_memory);
//
struct blasfeo_dvec *ocp_nlp_dynamics_disc_memory_get_fun_ptr(void *memory);
//
struct blasfeo_dvec *ocp_nlp_dynamics_disc_memory_get_adj_ptr(void *memory);


/************************************************
 * workspace
 ************************************************/

typedef struct
{
    struct blasfeo_dmat tmp_nv_nv;
    struct blasfeo_dvec adj_dyn_ux_pdiff;
} ocp_nlp_dynamics_disc_workspace;

acados_size_t ocp_nlp_dynamics_disc_workspace_calculate_size(void *config, void *dims, void *opts);



/************************************************
 * model
 ************************************************/

typedef struct
{
    external_function_generic *disc_dyn_fun;
    external_function_generic *disc_dyn_fun_jac;
    external_function_generic *disc_dyn_fun_jac_hess;
    external_function_generic *disc_dyn_phi_jac_p_hess_xu_p;
    external_function_generic *disc_dyn_phi_hess_ux_pdiff_adj_pdiff;
    external_function_generic *disc_dyn_adj_p;
} ocp_nlp_dynamics_disc_model;

//
acados_size_t ocp_nlp_dynamics_disc_model_calculate_size(void *config, void *dims);
//
void *ocp_nlp_dynamics_disc_model_assign(void *config, void *dims, void *raw_memory);
//
void ocp_nlp_dynamics_disc_model_set(void *config_, void *dims_, void *model_, const char *field, void *value);



/************************************************
 * functions
 ************************************************/

//
void ocp_nlp_dynamics_disc_config_initialize_default(void *config, int stage);
//
void ocp_nlp_dynamics_disc_initialize(void *config_, void *dims, void *model_, void *opts, void *mem, void *work_);
//
void ocp_nlp_dynamics_disc_update_qp_matrices(void *config_, void *dims, void *model_, void *opts, void *mem, void *work_);
//
void ocp_nlp_dynamics_disc_compute_fun(void *config_, void *dims, void *model_, void *opts, void *mem, void *work_);
//
void ocp_nlp_dynamics_disc_compute_jac_hess_p(void *config_, void *dims, void *model_, void *opts, void *mem, void *work_);
//
void ocp_nlp_dynamics_disc_compute_adj_p(void* config_, void *dims_, void *model_, void *opts_, void *mem_, struct blasfeo_dvec *out);
//
void ocp_nlp_dynamics_disc_reset(void *config_, void *dims_, void *model_, void *opts_, void *mem_, void *work_);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_OCP_NLP_OCP_NLP_DYNAMICS_DISC_H_
/// @}
/// @}
