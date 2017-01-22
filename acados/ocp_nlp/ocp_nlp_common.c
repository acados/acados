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

#include <stdlib.h>
#include <string.h>

#include "acados/ocp_nlp/ocp_nlp_common.h"

void ocp_nlp_create_memory(const ocp_nlp_in *in, ocp_nlp_mem *mem) {
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
