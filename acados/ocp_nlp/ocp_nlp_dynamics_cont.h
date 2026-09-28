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

#ifndef ACADOS_OCP_NLP_OCP_NLP_DYNAMICS_CONT_H_
#define ACADOS_OCP_NLP_OCP_NLP_DYNAMICS_CONT_H_

#ifdef __cplusplus
extern "C" {
#endif



// blasfeo
#include "blasfeo_common.h"

// acados
#include "acados/ocp_nlp/ocp_nlp_dynamics_common.h"
#include "acados/utils/external_function_generic.h"
#include "acados/utils/types.h"
#include "acados_c/sim_interface.h"



/************************************************
 * dims
 ************************************************/

typedef struct
{
    void *sim;
    int nx;   // number of states at the current stage
    int nz;   // number of algebraic states at the current stage
    int nu;   // number of inputs at the current stage
    int nx1;  // number of states at the next stage
    int nu1;  // number of inputes at the next stage
    int np;   // number of parameters at the current stage
} ocp_nlp_dynamics_cont_dims;

//
acados_size_t ocp_nlp_dynamics_cont_dims_calculate_size(void *config);
//
void *ocp_nlp_dynamics_cont_dims_assign(void *config, void *raw_memory);
//
void ocp_nlp_dynamics_cont_dims_set(void *config_, void *dims_, const char *field, int* value);

/************************************************
 * options
 ************************************************/

typedef struct
{
    void *sim_solver;
    int compute_adj;
    int compute_hess;
} ocp_nlp_dynamics_cont_opts;

//
acados_size_t ocp_nlp_dynamics_cont_opts_calculate_size(void *config, void *dims);
//
void *ocp_nlp_dynamics_cont_opts_assign(void *config, void *dims, void *raw_memory);
//
void ocp_nlp_dynamics_cont_opts_initialize_default(void *config, void *dims, void *opts);
//
void ocp_nlp_dynamics_cont_opts_update(void *config, void *dims, void *opts);
//
void ocp_nlp_dynamics_cont_opts_set(void *config, void *opts, const char *field, void* value);



/************************************************
 * memory
 ************************************************/

typedef struct
{
    struct blasfeo_dvec fun;
    struct blasfeo_dvec adj;
    struct blasfeo_dvec *ux;            // pointer to ux in nlp_out at current stage
    struct blasfeo_dvec *ux1;           // pointer to ux in nlp_out at next stage
    struct blasfeo_dvec *pi;            // pointer to pi in nlp_out at current stage
    struct blasfeo_dmat *BAbt;          // pointer to BAbt in qp_in
    struct blasfeo_dmat *RSQrq;         // pointer to RSQrq in qp_in
    struct blasfeo_dvec *z_alg;         // pointer to output z at t = 0
    bool *set_sim_guess;                 // indicate if initialization for integrator is set from outside
    void *cost_capsule;  // pointer to ocp_nlp_cost_capsule of the cost module
    struct blasfeo_dvec *sim_guess;     // initializations for integrator
    // struct blasfeo_dvec *z;             // pointer to (input) z in nlp_out at current stage
    struct blasfeo_dmat *dzduxt;        // pointer to dzdux transposed
    void *sim_solver;                   // sim solver memory
    acados_size_t workspace_size;
    acados_size_t sim_workspace_size;

} ocp_nlp_dynamics_cont_memory;

//
acados_size_t ocp_nlp_dynamics_cont_memory_calculate_size(void *config, void *dims, void *opts);
//
void *ocp_nlp_dynamics_cont_memory_assign(void *config, void *dims, void *opts, void *raw_memory);
//
struct blasfeo_dvec *ocp_nlp_dynamics_cont_memory_get_fun_ptr(void *memory);
//
struct blasfeo_dvec *ocp_nlp_dynamics_cont_memory_get_adj_ptr(void *memory);



/************************************************
 * workspace
 ************************************************/

typedef struct
{
    struct blasfeo_dmat hess;
    sim_in *sim_in;
    sim_out *sim_out;
    void *sim_solver;  // sim solver workspace
} ocp_nlp_dynamics_cont_workspace;

acados_size_t ocp_nlp_dynamics_cont_workspace_calculate_size(void *config, void *dims, void *opts);



/************************************************
 * model
 ************************************************/

typedef struct
{
    void *sim_model;
    // double *state_transition; // TODO
    double T;  // simulation time
} ocp_nlp_dynamics_cont_model;

//
acados_size_t ocp_nlp_dynamics_cont_model_calculate_size(void *config, void *dims);
//
void *ocp_nlp_dynamics_cont_model_assign(void *config, void *dims, void *raw_memory);
//
void ocp_nlp_dynamics_cont_model_set(void *config_, void *dims_, void *model_, const char *field, void *value);



/************************************************
 * functions
 ************************************************/

//
void ocp_nlp_dynamics_cont_config_initialize_default(void *config, int stage);
//
void ocp_nlp_dynamics_cont_initialize(void *config_, void *dims, void *model_, void *opts, void *mem, void *work_);
//
void ocp_nlp_dynamics_cont_update_qp_matrices(void *config_, void *dims, void *model_, void *opts, void *mem, void *work_);
//
void ocp_nlp_dynamics_cont_compute_fun(void *config_, void *dims, void *model_, void *opts, void *mem, void *work_);
//
void ocp_nlp_dynamics_cont_compute_fun_and_adj(void *config_, void *dims, void *model_, void *opts, void *mem, void *work_);
//
int ocp_nlp_dynamics_cont_precompute(void *config_, void *dims, void *model_, void *opts_, void *mem_, void *work_);
//
void ocp_nlp_dynamics_cont_compute_jac_hess_p(void *config_, void *dims, void *model_, void *opts, void *mem, void *work_);
//
void ocp_nlp_dynamics_cont_compute_adj_p(void* config_, void *dims_, void *model_, void *opts_, void *mem_, struct blasfeo_dvec *out);
//
void ocp_nlp_dynamics_cont_reset(void *config_, void *dims_, void *model_, void *opts_, void *mem_, void *work_);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_OCP_NLP_OCP_NLP_DYNAMICS_CONT_H_
/// @}
/// @}
