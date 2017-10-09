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

#ifndef ACADOS_OCP_NLP_OCP_NLP_SM_GN_H_
#define ACADOS_OCP_NLP_OCP_NLP_SM_GN_H_

#include "acados/ocp_nlp/ocp_nlp_common.h"
#include "acados/sim/sim_common.h"
#include "acados/utils/types.h"

typedef struct {
    void (*fun)(const real_t **, real_t **, int *, real_t *, int);
    void (*jac_fun)(const real_t **, real_t **, int *, real_t *, int);
    const int_t *ny;
    real_t *W;
    real_t *y_ref;
} ocp_nlp_ls_cost;

typedef struct {
    int_t dummy;
} ocp_nlp_sm_gn_args;

typedef struct {
    sim_solver *sim;
    ocp_nlp_ls_cost *ls_cost;
} ocp_nlp_sm_gn_memory;

typedef struct {
    real_t *w; // should be max_i (nx[i]+nu[i])
    real_t *F; // should be max_i ny[i]
    real_t *DF; // should be max_i ny[i]*(nx[i]+nu[i])
    real_t *DFT;  // should be max_i (nx[i]*nu[i])*ny[i]
    real_t *DFTW; // should be max_i (nx[i]+nu[i])*ny[i]
} ocp_nlp_sm_gn_workspace;

// ocp_nlp_sm_gn_args *ocp_nlp_sm_gn_create_arguments();

// int_t ocp_nlp_sm_gn_calculate_memory_size(const ocp_nlp_in *nlp_in,
//                                                     void *args_);

// char *ocp_nlp_sm_gn_assign_memory(const ocp_nlp_in *nlp_in, void *args_,
//                                             void **mem_, void *raw_memory);

// ocp_nlp_sm_gn_memory *ocp_nlp_sm_gn_create_memory(
//     const ocp_nlp_in *nlp_in, void *args_);

// int_t ocp_nlp_sm_gn_calculate_workspace_size(const ocp_nlp_in *nlp_in,
//                                                        void *args_);

int_t ocp_nlp_sm_gn(const ocp_nlp_in *nlp_in, ocp_nlp_memory *nlp_mem,
                              void *args_, void *memory_, void *workspace_);

// void ocp_nlp_sm_gn_initialize(const ocp_nlp_in *nlp_in, void *args_,
//                                         void **mem, void **work);

// void ocp_nlp_sm_gn_destroy(void *mem, void *work);

#endif  // ACADOS_OCP_NLP_OCP_NLP_SM_GN_H_