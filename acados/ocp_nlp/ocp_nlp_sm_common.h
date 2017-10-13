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

#ifndef ACADOS_OCP_NLP_OCP_NLP_SM_COMMON_H_
#define ACADOS_OCP_NLP_OCP_NLP_SM_COMMON_H_

#include "acados/sim/sim_common.h"
#include "acados/utils/types.h"
#include "acados/utils/casadi_wrapper.h"

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
    int_t N;
    const int_t *nx;
    const int_t *nu;
    const int_t *nb;
    const int_t *ng;
    
    void *cost;
    sim_solver *sim;
    ocp_nlp_function *path_constraints;
    
    const real_t **x;
    const real_t **u;
    const real_t **pi;
    const real_t **lam;
} ocp_nlp_sm_in;

typedef struct {
    const real_t **hess_l;  // TODO(nielsvd): Hessians of stage-wise Lagrangians, document precise definition.
    const real_t **grad_f;  // Gradients of stage-wise cost terms.
    const real_t **jac_h;   // TODO(niels): rename (maybe phi?). Jacobians of stage-wise integration operator.
    const real_t **jac_g;   // Jacobians of stage-wise path constraints.
    const real_t **h;       // TODO(nielsvd): rename. Evaluation of stage-wise integration operator.
    const real_t **g;       // Evaluation of stage-wise path constraints.
} ocp_nlp_sm_out;

typedef struct {
    void (*fun)(const ocp_nlp_sm_in *sm_in, ocp_nlp_sm_out *sm_out,
                void *args, void *mem, void *work);
    void (*initialize)(const ocp_nlp_sm_in *sm_in, void *args, void **mem,
                       void **work);
    void (*destroy)(void *mem, void *work);
    void *args;
    void *mem;
    void *work;
} ocp_nlp_sm;

#endif  // ACADOS_OCP_NLP_OCP_NLP_SM_COMMON_H_