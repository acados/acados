/*
 * Copyright (c) The acados authors.
 *
 * This file is part of acados.
 *
 * Licensed under the 2-Clause BSD License.
 */


/// \addtogroup ocp_nlp
/// @{
/// \addtogroup ocp_nlp_reg
/// @{

#ifndef ACADOS_OCP_NLP_OCP_NLP_REG_PROJECT_REDUC_HESS_H_
#define ACADOS_OCP_NLP_OCP_NLP_REG_PROJECT_REDUC_HESS_H_

#ifdef __cplusplus
extern "C" {
#endif



// blasfeo
#include "blasfeo_common.h"

// acados
#include "acados/ocp_nlp/ocp_nlp_reg_common.h"



/************************************************
 * dims
 ************************************************/

// use the functions in ocp_nlp_reg_common

/************************************************
 * options
 ************************************************/

typedef struct
{
    double thr_eig;
    double min_eig;
    double min_pivot;
    int pivoting;
} ocp_nlp_reg_project_reduc_hess_opts;

//
acados_size_t ocp_nlp_reg_project_reduc_hess_opts_calculate_size(void);
//
void *ocp_nlp_reg_project_reduc_hess_opts_assign(void *raw_memory);
//
void ocp_nlp_reg_project_reduc_hess_opts_initialize_default(void *config_, ocp_nlp_reg_dims *dims, void *opts_);
//
void ocp_nlp_reg_project_reduc_hess_opts_set(void *config_, void *opts_, const char *field, void* value);



/************************************************
 * memory
 ************************************************/

typedef struct
{
    double *reg_hess; // TODO move to workspace
    double *V; // TODO move to workspace
    double *d; // TODO move to workspace
    double *e; // TODO move to workspace

    // giaf's
    struct blasfeo_dmat L; // TODO move to workspace
    struct blasfeo_dmat L2; // TODO move to workspace
    struct blasfeo_dmat L3; // TODO move to workspace
    struct blasfeo_dmat Ls; // TODO move to workspace
    struct blasfeo_dmat P; // TODO move to workspace
    struct blasfeo_dmat AL; // TODO move to workspace

    struct blasfeo_dmat **RSQrq;  // pointer to RSQrq in qp_in
    struct blasfeo_dmat **BAbt;  // pointer to RSQrq in qp_in
} ocp_nlp_reg_project_reduc_hess_memory;

//
acados_size_t ocp_nlp_reg_project_reduc_hess_memory_calculate_size(void *config, ocp_nlp_reg_dims *dims, void *opts);
//
void *ocp_nlp_reg_project_reduc_hess_memory_assign(void *config, ocp_nlp_reg_dims *dims, void *opts, void *raw_memory);

/************************************************
 * workspace
 ************************************************/

 // TODO

/************************************************
 * functions
 ************************************************/

//
void ocp_nlp_reg_project_reduc_hess_config_initialize_default(ocp_nlp_reg_config *config);



#ifdef __cplusplus
}
#endif

#endif  // ACADOS_OCP_NLP_OCP_NLP_REG_PROJECT_REDUC_HESS_H_
/// @}
/// @}
