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

// external
#include <assert.h>
// hpipm
#include "hpipm/include/hpipm_d_ocp_qp.h"
#include "hpipm/include/hpipm_d_ocp_qp_sol.h"
#include "hpipm/include/hpipm_d_ocp_qp_dim.h"
// acados
#include "acados/utils/types.h"
#include "acados/ocp_qp/ocp_qp_common.h"



int ocp_qp_dims_calculate_size(int N)
{
    int size = sizeof(ocp_qp_dims);
    size += d_memsize_ocp_qp_dim(dims->N);
    return size;
}



ocp_qp_dims *assign_ocp_qp_dims(int N, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    ocp_qp_dims *dims = (ocp_qp_dims *) c_ptr;
    c_ptr += sizeof(ocp_qp_dims);

    d_create_ocp_qp_dim(N, dims, c_ptr);
    c_ptr += d_memsize_ocp_qp_dim(N);

    assert((char *) raw_memory + ocp_qp_dims_calculate_size(N) == c_ptr);

    return dims;
}



int ocp_qp_in_calculate_size(int N)
{
    int size = sizeof(ocp_qp_in);
    size += d_memsize_ocp_qp(dims);
    size += ocp_qp_dims_calculate_size(N);
    return size;
}



ocp_qp_in *assign_ocp_qp_in(ocp_qp_dims *dims, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    ocp_qp_in *qp_in = (ocp_qp_in *) c_ptr;
    c_ptr += sizeof(ocp_qp_in);

    d_create_ocp_qp(dims, qp_in, c_ptr);
    c_ptr += d_memsize_ocp_qp(dims);

    ocp_qp_dims *dims_copy = assign_ocp_qp_dims(dims->N, c_ptr);
    c_ptr += ocp_qp_in_calculate_size(dims->N);

    dims_copy->N = dims->N;

    for (int ii = 0; ii < dims->N+1; ii++)
    {
        dims_copy->nx[ii] = dims->nx[ii];
        dims_copy->nu[ii] = dims->nu[ii];
        dims_copy->nb[ii] = dims->nb[ii];
        dims_copy->ng[ii] = dims->ng[ii];
        dims_copy->ns[ii] = dims->ns[ii];
        dims_copy->nbu[ii] = dims->nbu[ii];
        dims_copy->nbx[ii] = dims->nbx[ii];
    }

    qp_in->dim = dims_copy;

    assert((char*) raw_memory + ocp_qp_in_calculate_size(dims) == c_ptr);

    return qp_in;
}



int ocp_qp_out_calculate_size(ocp_qp_dims *dims)
{
    int size = sizeof(ocp_qp_out);
    size += d_memsize_ocp_qp_sol(dims);
    return size;
}



ocp_qp_out *assign_ocp_qp_out(ocp_qp_dims *dims, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    ocp_qp_out *qp_out = (ocp_qp_out *) c_ptr;
    c_ptr += sizeof(ocp_qp_out);

    d_create_ocp_qp_sol(dims, qp_out, c_ptr);
    c_ptr += d_memsize_ocp_qp_sol(dims);

    assert((char*) raw_memory + ocp_qp_out_calculate_size(dims) == c_ptr);

    return qp_out;
}



// void form_nbu_nbx_rev(int N, int *nbu, int *nbx, int *nb, int* nx, int *nu, int **idxb_rev)
// {
//     for (int ii = 0; ii < N+1; ii++)
//     {
//         nbu[ii] = 0;
//         nbx[ii] = 0;
//         for (int jj = 0; jj < nb[ii]; jj++)
//         {
//             if (idxb_rev[ii][jj] < nx[ii])
//             {
//                 nbx[ii]++;
//             }
//             else
//             {
//                 nbu[ii]++;
//             }
//         }
//     }
// }