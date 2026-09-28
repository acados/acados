/*
 * Copyright (c) The acados authors.
 *
 * This file is part of acados.
 *
 * Licensed under the 2-Clause BSD License.
 */


/// \ingroup ocp_nlp
/// @{

/// \defgroup ocp_nlp_reg ocp_nlp_reg
/// @{

#ifndef ACADOS_OCP_NLP_OCP_NLP_REG_COMMON_H_
#define ACADOS_OCP_NLP_OCP_NLP_REG_COMMON_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "acados/ocp_qp/ocp_qp_common.h"



/* dims */

// same as qp_dims
typedef struct
{
    int *nx;
    int *nu;
    int *nbu;
    int *nbx;
    int *ng;
    int N;
} ocp_nlp_reg_dims;

//
acados_size_t ocp_nlp_reg_dims_calculate_size(int N);
//
ocp_nlp_reg_dims *ocp_nlp_reg_dims_assign(int N, void *raw_memory);
//
void ocp_nlp_reg_dims_set(void *config_, ocp_nlp_reg_dims *dims, int stage, char *field, int* value);



/* config */

typedef struct
{
    /* dims */
    acados_size_t (*dims_calculate_size)(int N);
    ocp_nlp_reg_dims *(*dims_assign)(int N, void *raw_memory);
    void (*dims_set)(void *config, ocp_nlp_reg_dims *dims, int stage, char *field, int *value);
    /* opts */
    acados_size_t (*opts_calculate_size)(void);
    void *(*opts_assign)(void *raw_memory);
    void (*opts_initialize_default)(void *config, ocp_nlp_reg_dims *dims, void *opts);
    void (*opts_set)(void *config, void *opts, const char *field, void* value);
    /* memory */
    acados_size_t (*memory_calculate_size)(void *config, ocp_nlp_reg_dims *dims, void *opts);
    void *(*memory_assign)(void *config, ocp_nlp_reg_dims *dims, void *opts, void *raw_memory);
    void (*memory_set)(void *config, ocp_nlp_reg_dims *dims, void *memory, char *field, void* value);
    /* functions */
    void (*regularize)(void *config, ocp_nlp_reg_dims *dims, void *opts, void *memory);
    void (*regularize_lhs)(void *config, ocp_nlp_reg_dims *dims, void *opts, void *memory);
    void (*regularize_rhs)(void *config, ocp_nlp_reg_dims *dims, void *opts, void *memory);
    void (*correct_dual_sol)(void *config, ocp_nlp_reg_dims *dims, void *opts, void *memory);
} ocp_nlp_reg_config;

//
acados_size_t ocp_nlp_reg_config_calculate_size(void);
//
void *ocp_nlp_reg_config_assign(void *raw_memory);



/* regularization help functions */
void acados_reconstruct_A(int dim, double *A, double *V, double *d);
void acados_mirror(int dim, double *A, double *V, double *d, double *e, double epsilon);
void acados_mirror_adaptive_eps(int dim, double *A, double *V, double *d, double *e, double max_cond_block, double min_eps);
void acados_project(int dim, double *A, double *V, double *d, double *e, double epsilon);
void acados_project_adaptive_eps(int dim, double *A, double *V, double *d, double *e, double max_cond_block, double min_eps);


#ifdef __cplusplus
}
#endif

#endif  // ACADOS_OCP_NLP_OCP_NLP_REG_COMMON_H_
/// @}
/// @}
