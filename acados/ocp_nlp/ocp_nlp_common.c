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
#include <assert.h>

#include "acados/utils/mem.h"
#include "hpipm/include/hpipm_d_ocp_qp_size.h"


static int number_of_primal_vars(ocp_nlp_dims *dims)
{
    int num_vars = 0;
    for (int ii = 0; ii <= dims->N; ii++) {
        num_vars += dims->nx[ii] + dims->nu[ii];
    }
    return num_vars;
}



void cast_nlp_dims_to_qp_dims(ocp_qp_dims *qp_dims, ocp_nlp_dims *nlp_dims)
{
    qp_dims->N = nlp_dims->N;
    qp_dims->nx = nlp_dims->nx;
    qp_dims->nu = nlp_dims->nu;
    qp_dims->nb = nlp_dims->nb;
    qp_dims->nbx = nlp_dims->nbx;
    qp_dims->nbu = nlp_dims->nbu;
    qp_dims->ng = nlp_dims->ng;
    qp_dims->ns = nlp_dims->ns;

    // TODO(dimitris): probably redundant (can also remove hpipm header)
    qp_dims->memsize = d_memsize_ocp_qp_size(qp_dims->N);
}



int ocp_nlp_calculate_memory_size(ocp_nlp_dims *dims, ocp_nlp_args *args)
{
    int N = dims->N;

    int size = sizeof(ocp_nlp_memory);

    size += sizeof(double *) * (N + 1);  // x
    size += sizeof(double *) * (N + 1);  // u
    size += sizeof(double *) * (N + 1);  // lam
    size += sizeof(double *) * N;  // pi

    for (int ii = 0; ii <= N; ii++)
    {
        size += sizeof(double)*dims->nx[ii];  // x
        size += sizeof(double)*dims->nu[ii];  // u
        size += sizeof(double)*2*(dims->nb[ii] + dims->ng[ii] + dims->nh[ii]);  // lam
        if (ii < N)
        {
            size += sizeof(double)*dims->nx[ii+1];  // pi
        }
    }

    make_int_multiple_of(64, &size);
    size += 1 * 64;

    return size;
}



ocp_nlp_memory *ocp_nlp_assign_memory(ocp_nlp_dims *dims, ocp_nlp_args *args, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    // TODO(dimitris): check if aligning pointers is really necessary and use such asserts
    // assert((size_t)c_ptr % 8 == 0);

    int N = dims->N;

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
        // mem->x[ii] = (double *)c_ptr;
        // c_ptr += dims->nx[ii]*sizeof(double)+100;
        // mem->u[ii] = (double *)c_ptr;
        // c_ptr += dims->nu[ii]*sizeof(double);
        assign_double(dims->nx[ii], &mem->x[ii], &c_ptr);
        assign_double(dims->nu[ii], &mem->u[ii], &c_ptr);
        if (ii < N)
        {
            // mem->pi[ii] = (double *)c_ptr;
            // c_ptr += dims->nx[ii+1]*sizeof(double)+100;
            assign_double(dims->nx[ii+1], &mem->pi[ii], &c_ptr);
        }
        // mem->lam[ii] = (double *)c_ptr;
        // c_ptr += 2*(dims->nb[ii] + dims->ng[ii] + dims->nh[ii])*sizeof(double);
        assign_double(2*(dims->nb[ii] + dims->ng[ii] + dims->nh[ii]), &mem->lam[ii], &c_ptr);
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



int ocp_nlp_out_calculate_size(ocp_nlp_dims *dims, ocp_nlp_args *args)
{
    int N = dims->N;

    int size = sizeof(ocp_nlp_out);

    size += sizeof(double *) * (N + 1);  // x
    size += sizeof(double *) * (N + 1);  // u
    size += sizeof(double *) * (N + 1);  // lam
    size += sizeof(double *) * N;  // pi

    for (int ii = 0; ii <= N; ii++)
    {
        size += sizeof(double)*dims->nx[ii];  // x
        size += sizeof(double)*dims->nu[ii];  // u
        size += sizeof(double)*2*(dims->nb[ii] + dims->ng[ii] + dims->nh[ii]);  // lam
        if (ii < N)
        {
            size += sizeof(double)*dims->nx[ii+1];  // pi
        }
    }

    make_int_multiple_of(64, &size);
    size += 1 * 64;

    return size;
}



ocp_nlp_out *ocp_nlp_out_assign(ocp_nlp_dims *dims, ocp_nlp_args *args, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    int N = dims->N;

    ocp_nlp_out *out = (ocp_nlp_out *)c_ptr;
    c_ptr += sizeof(ocp_nlp_out);

    // double pointers
    out->x = (double **)c_ptr;
    c_ptr += sizeof(double *) * (N + 1);

    out->u = (double **)c_ptr;
    c_ptr += sizeof(double *) * (N + 1);

    out->pi = (double **)c_ptr;
    c_ptr += sizeof(double *) * N;

    out->lam = (double **)c_ptr;
    c_ptr += sizeof(double *) * (N + 1);

    // doubles
    align_char_to(64, &c_ptr);

    for (int ii = 0; ii <= N; ii++)
    {
        assign_double(dims->nx[ii], &out->x[ii], &c_ptr);
        assign_double(dims->nu[ii], &out->u[ii], &c_ptr);
        if (ii < N)
        {
            assign_double(dims->nx[ii+1], &out->pi[ii], &c_ptr);
        }
        assign_double(2*(dims->nb[ii] + dims->ng[ii] + dims->nh[ii]), &out->lam[ii], &c_ptr);
    }

    assert((char *)raw_memory + ocp_nlp_calculate_memory_size(dims, args) >= c_ptr);

    return out;
}
