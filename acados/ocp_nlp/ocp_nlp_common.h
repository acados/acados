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
    int *nh;  // number of path constraints - ONLY difference with ocp_qp_dims atm
    int *ns;  // number of soft constraints
    int N;
} ocp_nlp_dims;



typedef struct {
    // const int *sparsity;
    // const int *idx_in;
    int dim_in;
    int dim_out;
    void (*fun)(const double *, double *);
    void (*jac_fun)(const double *, double *);
    // TODO(rien): other directional and second order derivatives
    // TODO(rien): define the overlapping 'sets' of functions, jacobians, hessians etc..
} ocp_nlp_function;



typedef struct {
    // TODO(rien): only for least squares cost with state and control reference atm
    // void *fun;
    // const int *ny;
    double **W;
    double **y_ref;
} ocp_nlp_ls_cost;



typedef struct {
    ocp_nlp_function *fun;
} ocp_nlp_stage_cost;



typedef struct {
    ocp_nlp_dims *dims;

    // TODO(dimitris): decide on the blasfeo format for those fields
    int **idxb;
    double **lb;
    double **ub;

    double **Cx;
    double **Cu;
    double **lg;
    double **ug;

    double **lh;
    double **uh;
    ocp_nlp_function *h;  // nonlinear path constraints

    void *cost;
    sim_solver *sim;
    // TODO(rien): what about invariants, e.g., algebraic constraints?

    bool freezeSens;
} ocp_nlp_in;



typedef struct {
    int maxIter;
} ocp_nlp_args;



typedef struct {
    int num_vars;
    double **x;
    double **u;
    double **pi;
    double **lam;
} ocp_nlp_memory;



typedef struct {
    double *w;
} ocp_nlp_work;



typedef struct {
    double **x;
    double **u;
    double **pi;
    double **lam;
} ocp_nlp_out;

void cast_nlp_dims_to_qp_dims(ocp_qp_dims *qp_dims, ocp_nlp_dims *nlp_dims);
//
int ocp_nlp_calculate_args_size(ocp_nlp_dims *dims);
//
ocp_nlp_args *ocp_nlp_assign_args(ocp_nlp_dims *dims, void *raw_memory);
//
int ocp_nlp_in_calculate_size(ocp_nlp_dims *dims, ocp_nlp_args *args);
//
ocp_nlp_in *ocp_nlp_in_assign(ocp_nlp_dims *dims, ocp_nlp_args *args, int num_stages, void *raw_memory);
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
