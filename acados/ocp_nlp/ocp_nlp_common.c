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

    for (int_t i = 0; i < in->N; i++) {
        num_vars += in->nx[i] + in->nu[i];
    }
    num_vars += in->nx[in->N];

    return num_vars;
}


void ocp_nlp_create_memory(const ocp_nlp_in *in, ocp_nlp_mem *mem) {
    mem->num_vars = number_of_primal_vars(in);
    mem->x = (real_t **) malloc(sizeof(*mem->x) * (in->N+1));
    mem->u = (real_t **) malloc(sizeof(*mem->u) * in->N);
    mem->lam = (real_t **) malloc(sizeof(*mem->lam) * in->N);
    for (int_t i = 0; i < in->N; i++) {
        mem->x[i] = (real_t *) malloc(sizeof(*mem->x[i]) * (in->nx[i]));
        mem->u[i] = (real_t *) malloc(sizeof(*mem->u[i]) * (in->nu[i]));
        mem->lam[i] = (real_t *) malloc(sizeof(*mem->lam[i]) * (in->nx[i]));
    }
    mem->x[in->N] = (real_t *) malloc(sizeof(*mem->x[in->N]) * (in->nx[in->N]));
}


int_t ocp_nlp_calculate_workspace_size(const ocp_nlp_in *in, void *args_) {
    ocp_nlp_args *args = (ocp_nlp_args*) args_;

    int_t size;
    int_t num_vars = number_of_primal_vars(in);

    args->dummy = 1;

    size = sizeof(ocp_nlp_work);
    size += num_vars*sizeof(real_t);

    return size;
}


void ocp_nlp_cast_workspace(ocp_nlp_work *work, ocp_nlp_mem *mem) {
    char *ptr = (char *)work;

    ptr += sizeof(ocp_nlp_work);
    work->w = (real_t*)ptr;
    ptr += (mem->num_vars)*sizeof(real_t);
}
