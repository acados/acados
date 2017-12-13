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

#ifdef __cplusplus
extern "C" {
#endif

#include "acados/utils/types.h"

typedef struct {
    real_t **hess_l;  // TODO(nielsvd): Hessians of stage-wise
                            // Lagrangians, document precise definition.
    real_t **grad_f;  // Gradients of stage-wise cost terms.
    real_t **jac_h;   // TODO(niels): rename (maybe phi?). Jacobians of
                            // stage-wise integration operator.
    real_t **jac_g;   // Jacobians of stage-wise path constraints.
    real_t **h;       // TODO(nielsvd): rename (maybe phi?). Evaluation of
                            // stage-wise integration operator.
    real_t **g;       // Evaluation of stage-wise path constraints.
    real_t **x;
    real_t **u;
    real_t **pi;
    real_t **lam;
    // TODO(nielsvd): what about parameters and addtionaly variables such as
    // lifted variables?
} ocp_nlp_memory;

typedef struct {
    int_t N;
    const int_t *nx;
    const int_t *nu;
    const int_t *nb;
    const int_t *ng;
    const int_t **idxb;
    const real_t **lb;
    const real_t **ub;
    const real_t **lg;
    const real_t **ug;

    void *cost;
    void **sim;
    void **path_constraints;
} ocp_nlp_in;

typedef struct {
    real_t **x;
    real_t **u;
    real_t **pi;
    real_t **lam;
} ocp_nlp_out;

typedef struct {
    int_t (*fun)(const ocp_nlp_in *, ocp_nlp_out *, void *args, void *mem, void *work);
    void (*initialize)(const ocp_nlp_in *nlp_in, void *args, void **mem, void **work);
    void (*destroy)(void *mem, void *work);
    ocp_nlp_in *nlp_in;
    ocp_nlp_out *nlp_out;
    void *args;
    void *mem;
    void *work;
} ocp_nlp_solver;

int_t ocp_nlp_calculate_memory_size(const ocp_nlp_in *nlp_in);
char *ocp_nlp_assign_memory(const ocp_nlp_in *nlp_in, void **mem_,
                            void *raw_memory);
ocp_nlp_memory *ocp_nlp_create_memory(const ocp_nlp_in *nlp_in);

void ocp_nlp_destroy(void *mem_);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_OCP_NLP_OCP_NLP_COMMON_H_
