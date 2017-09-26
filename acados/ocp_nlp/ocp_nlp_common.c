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

#include "acados/ocp_nlp/ocp_nlp_common.h"

#include <stdlib.h>
#include <string.h>

static int_t number_of_primal_vars(const ocp_nlp_in *in) {
    int_t num_vars = 0;
    for (int_t i = 0; i <= in->N; i++) {
        num_vars += in->nx[i] + in->nu[i];
    }
    return num_vars;
}

void ocp_nlp_create_memory(const ocp_nlp_in *in, void *mem_) {
    ocp_nlp_memory *mem = (ocp_nlp_memory *)mem_;
    mem->num_vars = number_of_primal_vars(in);
    mem->x = (real_t **)calloc(in->N + 1, sizeof(*mem->x));
    mem->u = (real_t **)calloc(in->N + 1, sizeof(*mem->u));
    mem->pi = (real_t **)calloc(in->N, sizeof(*mem->pi));
    mem->lam = (real_t **)calloc(in->N + 1, sizeof(*mem->lam));
    for (int_t i = 0; i < in->N; i++) {
        mem->x[i] = (real_t *)calloc(in->nx[i], sizeof(*mem->x[i]));
        mem->u[i] = (real_t *)calloc(in->nu[i], sizeof(*mem->u[i]));
        mem->pi[i] = (real_t *)calloc(in->nx[i+1], sizeof(*mem->pi[i]));
        mem->lam[i] = (real_t *)calloc(2 * (in->nb[i] + in->nc[i] + in->ng[i]),
            sizeof(*mem->lam[i]));
    }
    mem->x[in->N] = (real_t *)calloc(in->nx[in->N], sizeof(*mem->x[in->N]));
    mem->u[in->N] = (real_t *)calloc(in->nu[in->N], sizeof(*mem->u[in->N]));
    mem->lam[in->N] = (real_t *)calloc(2 * (in->nb[in->N] + in->nc[in->N] + in->ng[in->N]),
        sizeof(*mem->lam[in->N]));
}

void ocp_nlp_free_memory(int_t N, void *mem_) {
    ocp_nlp_memory *mem = (ocp_nlp_memory *)mem_;

    for (int_t i = 0; i < N; i++) {
        free(mem->x[i]);
        free(mem->u[i]);
        free(mem->pi[i]);
        free(mem->lam[i]);
    }
    free(mem->x[N]);
    free(mem->u[N]);
    free(mem->lam[N]);
    free(mem->x);
    free(mem->u);
    free(mem->pi);
    free(mem->lam);
}

int_t ocp_nlp_calculate_workspace_size(const ocp_nlp_in *in, void *args_) {
    ocp_nlp_args *args = (ocp_nlp_args *)args_;

    int_t size;
    int_t num_vars = number_of_primal_vars(in);

    args->dummy = 1;

    size = sizeof(ocp_nlp_work);
    size += num_vars * sizeof(real_t);

    // allocate mem for least-squares cost
    int_t *nr = ocp_nlp_in->ls_cost->nr;
    int_t nr_ = 0;
    for (i = 0; i < ocp_nlp_in->N; i++ ) nr_+=nr[i];

    return size;
}

void ocp_nlp_cast_workspace(ocp_nlp_work *work, ocp_nlp_memory *mem) {
    char *ptr = (char *)work;

    ptr += sizeof(ocp_nlp_work);
    work->w = (real_t *)ptr;
    ptr += mem->num_vars * sizeof(real_t);  // TODO(robin): ptr never used again?
}
