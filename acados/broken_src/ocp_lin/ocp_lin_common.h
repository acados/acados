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

#ifndef ACADOS_OCP_LIN_OCP_LIN_COMMON_H_
#define ACADOS_OCP_LIN_OCP_LIN_COMMON_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "acados/sim/sim_common.h"
#include "acados/utils/casadi_wrapper.h"
#include "acados/utils/types.h"

typedef struct {
    int_t nx;
    int_t nu;
    int_t np;

    int_t ny;

    casadi_wrapper_in *in;
    casadi_wrapper_out *out;
    casadi_wrapper_args *args;
    casadi_wrapper_workspace *work;
} ocp_nlp_function;

typedef struct {
    int *nx;
    int *nu;
    int *nb;  // nbx + nbu
    int *nbx;
    int *nbu;
    int *ng;  // number of general linear constraints
    int *nh;  // number of path constraints - ONLY difference with ocp_qp_dims
              // atm
    int *ns;  // number of soft constraints
    int *num_stages;
    int N;
} ocp_lin_dims;

typedef struct {
    void *cost;
    sim_solver_fcn_ptrs **sim;
    ocp_nlp_function **path_constraints;

    const real_t **x;
    const real_t **u;
    const real_t **pi;
    const real_t **lam;

    // TODO(nielsvd): should go, old interface
    bool freezeSens;
} ocp_lin_in;

typedef struct {
    const real_t **hess_l;  // TODO(nielsvd): Hessians of stage-wise
                            // Lagrangians, document precise definition.
    const real_t **grad_f;  // Gradients of stage-wise cost terms.
    const real_t **jac_h;   // TODO(niels): rename (maybe phi?). Jacobians of
                            // stage-wise integration operator.
    const real_t **jac_g;   // Jacobians of stage-wise path constraints.
    const real_t **h;       // TODO(nielsvd): rename. Evaluation of stage-wise
                            // integration operator.
    const real_t **g;       // Evaluation of stage-wise path constraints.
} ocp_lin_out;

typedef struct {
    int (*fun)(ocp_lin_in *qp_in, ocp_lin_out *qp_out, void *args, void *mem, void *work);
    int (*calculate_args_size)(ocp_lin_dims *dims);
    void *(*assign_args)(ocp_lin_dims *dims, void *raw_memory);
    void *(*copy_args)(ocp_lin_dims *dims, void *raw_memory, void *source_);
    void (*initialize_default_args)(void *args);
    int (*calculate_memory_size)(ocp_lin_dims *dims, void *args);
    void *(*assign_memory)(ocp_lin_dims *dims, void *args, void *raw_memory);
    int (*calculate_workspace_size)(ocp_lin_dims *dims, void *args);
} ocp_lin_method_fcn_ptrs;

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_OCP_LIN_OCP_LIN_COMMON_H_
