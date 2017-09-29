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

typedef struct {
    real_t f;
    real_t *grad_f;
    real_t *hess_lag;
    real_t *g;
    real_t *jac_g;
} ocp_nlp_stage_sensitivities;

typedef struct {
    real_t **x;
    real_t **u;
    real_t **pi;
    real_t **lam;
    ocp_nlp_stage_sensitivities *sens;
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

    ocp_nlp_sensitivities_provider *sensitivities_provider;
    ocp_qp_solver *qp_solver; // TODO(nielsvd): rename to something more general
} ocp_nlp_in;

typedef struct {
    int_t maxIter;
} ocp_nlp_args;

typedef struct {
    real_t *w;
} ocp_nlp_work;

typedef struct {
    real_t **x;
    real_t **u;
    real_t **pi;
    real_t **lam;
} ocp_nlp_out;

typedef struct {
    int_t (*fun)(const ocp_nlp_in *, ocp_nlp_out *, void *args, void *mem,
                 void *work);
    ocp_nlp_in *nlp_in;
    ocp_nlp_out *nlp_out;
    void *args;
    void *mem;
    void *work;
} ocp_nlp_solver;

void ocp_nlp_create_memory(const ocp_nlp_in *in, void *mem_);
void ocp_nlp_free_memory(int_t N, void *mem_);

int_t ocp_nlp_calculate_workspace_size(const ocp_nlp_in *in, void *args_);
void ocp_nlp_cast_workspace(ocp_nlp_work *work, ocp_nlp_memory *mem);

#endif  // ACADOS_OCP_NLP_OCP_NLP_COMMON_H_
