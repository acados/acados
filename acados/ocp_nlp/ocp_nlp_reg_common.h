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

//typedef ocp_qp_dims ocp_nlp_reg_dims;
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
int ocp_nlp_reg_dims_calculate_size(int N);
//
ocp_nlp_reg_dims *ocp_nlp_reg_dims_assign(int N, void *raw_memory);
//
void ocp_nlp_reg_dims_set(void *config_, ocp_nlp_reg_dims *dims, int stage, char *field, int* value);



/* config */

typedef struct
{
    /* dims */
    int (*dims_calculate_size)(int N);
    ocp_nlp_reg_dims *(*dims_assign)(int N, void *raw_memory);
    void (*dims_set)(void *config, ocp_nlp_reg_dims *dims, int stage, char *field, int *value);
    /* opts */
    int (*opts_calculate_size)(void);
    void *(*opts_assign)(void *raw_memory);
    void (*opts_initialize_default)(void *config, ocp_nlp_reg_dims *dims, void *opts);
    void (*opts_set)(void *config, ocp_nlp_reg_dims *dims, void *opts, char *field, void* value);
    /* memory */
    int (*memory_calculate_size)(void *config, ocp_nlp_reg_dims *dims, void *opts);
    void *(*memory_assign)(void *config, ocp_nlp_reg_dims *dims, void *opts, void *raw_memory);
    void (*memory_set)(void *config, ocp_nlp_reg_dims *dims, void *memory, char *field, void* value);
    void (*memory_set_RSQrq_ptr)(ocp_nlp_reg_dims *dims, struct blasfeo_dmat *mat, void *memory);
    void (*memory_set_rq_ptr)(ocp_nlp_reg_dims *dims, struct blasfeo_dvec *vec, void *memory);
    void (*memory_set_BAbt_ptr)(ocp_nlp_reg_dims *dims, struct blasfeo_dmat *mat, void *memory);
    void (*memory_set_b_ptr)(ocp_nlp_reg_dims *dims, struct blasfeo_dvec *vec, void *memory);
    void (*memory_set_idxb_ptr)(ocp_nlp_reg_dims *dims, int **idxb, void *memory);
    void (*memory_set_DCt_ptr)(ocp_nlp_reg_dims *dims, struct blasfeo_dmat *mat, void *memory);
    void (*memory_set_ux_ptr)(ocp_nlp_reg_dims *dims, struct blasfeo_dvec *vec, void *memory);
    void (*memory_set_pi_ptr)(ocp_nlp_reg_dims *dims, struct blasfeo_dvec *vec, void *memory);
    void (*memory_set_lam_ptr)(ocp_nlp_reg_dims *dims, struct blasfeo_dvec *vec, void *memory);
    /* functions */
    void (*regularize_hessian)(void *config, ocp_nlp_reg_dims *dims, void *opts, void *memory);
    void (*correct_dual_sol)(void *config, ocp_nlp_reg_dims *dims, void *opts, void *memory);
} ocp_nlp_reg_config;

//
int ocp_nlp_reg_config_calculate_size(void);
//
void *ocp_nlp_reg_config_assign(void *raw_memory);



/* regularization help functions */
void acados_reconstruct_A(int dim, double *A, double *V, double *d);
void acados_mirror(int dim, double *A, double *V, double *d, double *e, double epsilon);
void acados_project(int dim, double *A, double *V, double *d, double *e, double epsilon);



#ifdef __cplusplus
}
#endif

#endif  // ACADOS_OCP_NLP_OCP_NLP_REG_COMMON_H_
/// @}
/// @}
