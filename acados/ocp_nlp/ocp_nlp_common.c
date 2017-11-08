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


static int number_of_primal_vars(ocp_nlp_dims *dims)
{
    int num_vars = 0;
    for (int ii = 0; ii <= dims->N; ii++) {
        num_vars += dims->nx[ii] + dims->nu[ii];
    }
    return num_vars;
}



// TODO(dimitris): update to new convention!
void ocp_nlp_create_memory(ocp_nlp_dims *dims, void *mem_) {

    // TODO(dimitris): FIX THIS!!!
    int ng_old[dims->N+1];
    for (int i = 0; i < dims->N+1; i++) ng_old[i] = 0;

    ocp_nlp_memory *mem = (ocp_nlp_memory *)mem_;
    mem->num_vars = number_of_primal_vars(dims);
    mem->x = (double **)calloc(dims->N + 1, sizeof(*mem->x));
    mem->u = (double **)calloc(dims->N + 1, sizeof(*mem->u));
    mem->pi = (double **)calloc(dims->N, sizeof(*mem->pi));
    mem->lam = (double **)calloc(dims->N + 1, sizeof(*mem->lam));
    for (int i = 0; i < dims->N; i++) {
        mem->x[i] = (double *)calloc(dims->nx[i], sizeof(*mem->x[i]));
        mem->u[i] = (double *)calloc(dims->nu[i], sizeof(*mem->u[i]));
        mem->pi[i] = (double *)calloc(dims->nx[i+1], sizeof(*mem->pi[i]));
        mem->lam[i] = (double *)calloc(2 * (dims->nb[i] + dims->ng[i] + ng_old[i]),
            sizeof(*mem->lam[i]));
    }
    mem->x[dims->N] = (double *)calloc(dims->nx[dims->N], sizeof(*mem->x[dims->N]));
    mem->u[dims->N] = (double *)calloc(dims->nu[dims->N], sizeof(*mem->u[dims->N]));
    mem->lam[dims->N] = (double *)calloc(2 * (dims->nb[dims->N] + dims->ng[dims->N] + ng_old[dims->N]),
        sizeof(*mem->lam[dims->N]));
}



// TODO(dimitris): remove!
void ocp_nlp_free_memory(int N, void *mem_) {
    ocp_nlp_memory *mem = (ocp_nlp_memory *)mem_;

    for (int i = 0; i < N; i++) {
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



int ocp_nlp_calculate_workspace_size(ocp_nlp_dims *dims, ocp_nlp_args *args)
{
    int size = 0;
    int num_vars = number_of_primal_vars(dims);

    size += sizeof(ocp_nlp_work);
    size += num_vars * sizeof(double);  // w

    return size;
}



void ocp_nlp_cast_workspace(ocp_nlp_work *work, ocp_nlp_memory *mem)
{
    char *ptr = (char *)work;

    ptr += sizeof(ocp_nlp_work);
    work->w = (double *)ptr;
    ptr += mem->num_vars * sizeof(double);
    // TODO(dimitris): check that this pointer does not go outside allocated memory
}
