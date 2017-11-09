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
#include "acados/utils/mem.h"

#include <stdlib.h>
#include <string.h>
#include <assert.h>


static int number_of_primal_vars(ocp_nlp_dims *dims)
{
    int num_vars = 0;
    for (int ii = 0; ii <= dims->N; ii++) {
        num_vars += dims->nx[ii] + dims->nu[ii];
    }
    return num_vars;
}



// TODO(dimitris): THAT'S NOT REALLY MEMORY, THAT'S OUTPUT
int ocp_nlp_calculate_memory_size(ocp_nlp_dims *dims, ocp_nlp_args *args)
{
    int N = dims->N;

    // TODO(dimitris): ADD DIMENSION OF PATH CONSTRAINTS TO nlp dims !!!
    int ng_old[dims->N+1];
    for (int i = 0; i < dims->N+1; i++) ng_old[i] = 0;

    int size = sizeof(ocp_nlp_memory);

    size += sizeof(double *) * (N + 1);  // x
    size += sizeof(double *) * (N + 1);  // u
    size += sizeof(double *) * (N + 1);  // lam
    size += sizeof(double *) * N;  // pi

    for (int ii = 0; ii <= N; ii++)
    {
        size += sizeof(double)*dims->nx[ii];  // x
        size += sizeof(double)*dims->nu[ii];  // u
        size += sizeof(double)*(2*dims->nb[ii] + 2*dims->ng[ii] + 2*ng_old[ii]);  // lam
        if (ii < N)
        {
            size += sizeof(double)*dims->nx[ii+1];  // pi
        }
    }

    make_int_multiple_of(64, &size);
    size += 1 * 64;

    return size + 100*(N+1);  // TODO(dimitris): FIND BUG IN CODE WHERE WE WRITE OUTSIDE MEMORY!
}



ocp_nlp_memory *ocp_nlp_assign_memory(ocp_nlp_dims *dims, ocp_nlp_args *args, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    int N = dims->N;

    // TODO(dimitris): ADD DIMENSION OF PATH CONSTRAINTS TO nlp dims !!!
    int ng_old[N+1];
    for (int i = 0; i < N+1; i++) ng_old[i] = 0;

    ocp_nlp_memory *mem = (ocp_nlp_memory *)c_ptr;
    c_ptr += sizeof(ocp_nlp_memory);

    mem->num_vars = number_of_primal_vars(dims);

    // double pointers
    mem->x = (double **)c_ptr;
    c_ptr += sizeof(double *) * (N + 1);

    mem->u = (double **)c_ptr;
    c_ptr += sizeof(double *) * (N + 1);

    mem->pi = (double **)c_ptr;
    c_ptr += sizeof(double *) * N;

    mem->lam = (double **)c_ptr;
    c_ptr += sizeof(double *) * (N + 1);

    // doubles
    align_char_to(64, &c_ptr);

    for (int ii = 0; ii <= N; ii++)
    {
        mem->x[ii] = (double *)c_ptr;
        c_ptr += dims->nx[ii]+100;  // TODO(dimitris): FIND BUG IN CODE WHERE WE WRITE OUTSIDE MEMORY!
        mem->u[ii] = (double *)c_ptr;
        c_ptr += dims->nu[ii];
        if (ii < N)
        {
            mem->pi[ii] = (double *)c_ptr;
            c_ptr += dims->nx[ii+1];
        }
        mem->lam[ii] = (double *)c_ptr;
        c_ptr += 2*(dims->nb[ii] + dims->ng[ii] + ng_old[ii])+100;  // TODO(dimitris): FIND BUG IN CODE WHERE WE WRITE OUTSIDE MEMORY!
    }

    assert((char *)raw_memory + ocp_nlp_calculate_memory_size(dims, args) >= c_ptr);

    return mem;
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
