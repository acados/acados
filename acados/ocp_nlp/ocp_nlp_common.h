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

#ifndef ACADOS_OCP_NLP_OCP_NLP_COMMON_H_
#define ACADOS_OCP_NLP_OCP_NLP_COMMON_H_

#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/sim/sim_common.h"
#include "acados/utils/types.h"

#ifdef __cplusplus
extern "C" {
#endif


typedef struct {
    int *nx;
    int *nu;
    int *nb;  // nbx + nbu
    int *nbx;
    int *nbu;
    int *ng;  // number of general linear constraints
    int *ns;  // number of soft constraints
    int *nh;  // number of path constraints - ONLY difference with ocp_qp_dims atm
    int N;
} ocp_nlp_dims;



typedef struct {
    // const int_t *sparsity;
    // const int_t *idx_in;
    int_t dim_in;
    int_t dim_out;
    void (*fun)(const real_t *, real_t *);
    void (*jac_fun)(const real_t *, real_t *);
    // TODO(rien): other directional and second order derivatives
    // TODO(rien): define the overlapping 'sets' of functions, jacobians, hessians etc..
} ocp_nlp_function;



typedef struct {
    // TODO(rien): only for least squares cost with state and control reference atm
    // void *fun;
    // const int_t *ny;
    real_t **W;
    real_t **y_ref;
} ocp_nlp_ls_cost;



typedef struct {
    ocp_nlp_function *fun;
} ocp_nlp_stage_cost;



typedef struct {
    ocp_nlp_dims *dims;
    // TODO(dimitris): decide on the blasfeo format for those fields
    const int_t **idxb;
    const real_t **lb;
    const real_t **ub;
    const real_t **Cx;
    const real_t **Cu;
    const real_t **lc;
    const real_t **uc;
    const real_t **lg;
    const real_t **ug;

    void *cost;
    sim_solver *sim;
    ocp_nlp_function *g;  // nonlinear constraints
    // TODO(rien): what about invariants, e.g., algebraic constraints?

    bool freezeSens;
} ocp_nlp_in;



typedef struct {
    int_t maxIter;
} ocp_nlp_args;



typedef struct {
    int_t num_vars;
    real_t **x;
    real_t **u;
    real_t **pi;
    real_t **lam;
} ocp_nlp_memory;



typedef struct {
    real_t *w;
} ocp_nlp_work;



typedef struct {
    real_t **x;
    real_t **u;
    real_t **pi;
    real_t **lam;
} ocp_nlp_out;

void cast_nlp_dims_to_qp_dims(ocp_qp_dims *qp_dims, ocp_nlp_dims *nlp_dims);
//
int ocp_nlp_calculate_args_size(ocp_nlp_dims *dims);
//
ocp_nlp_args *ocp_nlp_assign_args(ocp_nlp_dims *dims, void *raw_memory);
//
// TODO(dimitris): write those after deciding on nlp_in format
// int ocp_nlp_in_calculate_size(ocp_nlp_dims *dims, ocp_nlp_args *args);
//
// ocp_nlp_in *ocp_nlp_in_assign(ocp_nlp_dims *dims, ocp_nlp_args *args, void *raw_memory);
//
int ocp_nlp_calculate_memory_size(ocp_nlp_dims *dims, ocp_nlp_args *args);
//
ocp_nlp_memory *ocp_nlp_assign_memory(ocp_nlp_dims *dims, ocp_nlp_args *args, void *raw_memory);
//
int ocp_nlp_calculate_workspace_size(ocp_nlp_dims *dims, ocp_nlp_args *args);
//
void ocp_nlp_cast_workspace(ocp_nlp_work *work, ocp_nlp_memory *mem);
//
int ocp_nlp_out_calculate_size(ocp_nlp_dims *dims, ocp_nlp_args *args);
//
ocp_nlp_out *ocp_nlp_out_assign(ocp_nlp_dims *dims, ocp_nlp_args *args, void *raw_memory);


#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_OCP_NLP_OCP_NLP_COMMON_H_
