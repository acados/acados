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
// blasfeo
#include "blasfeo/include/blasfeo_target.h"
#include "blasfeo/include/blasfeo_common.h"
#include "blasfeo/include/blasfeo_d_aux.h"
#include "blasfeo/include/blasfeo_d_aux_ext_dep.h"
// hpipm
#include "hpipm_d_dense_qp.h"
#include "hpipm_d_dense_qp_sol.h"
#include "hpipm_d_dense_qp_dim.h"
// acados
#include "acados/utils/types.h"
#include "acados/dense_qp/dense_qp_common.h"



int dense_qp_dims_calculate_size()
{
    int size = sizeof(dense_qp_dims);

    size += d_memsize_dense_qp_dim();

    return size;
}



dense_qp_dims *assign_dense_qp_dims(void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    dense_qp_dims *dims = (dense_qp_dims *) c_ptr;
    c_ptr += sizeof(dense_qp_dims);

    d_create_dense_qp_dim(dims, c_ptr);
    c_ptr += d_memsize_dense_qp_dim();

    assert((char *) raw_memory + dense_qp_dims_calculate_size() == c_ptr);

    return dims;
}



int dense_qp_in_calculate_size(dense_qp_dims *dims)
{
    int size = sizeof(dense_qp_in);
    size += sizeof(dense_qp_dims);
    size += d_memsize_dense_qp(dims);
    return size;
}



dense_qp_in *assign_dense_qp_in(dense_qp_dims *dims, void *raw_memory)
{
    char *c_ptr = (char *)raw_memory;

    dense_qp_in *qp_in = (dense_qp_in *) c_ptr;
    c_ptr += sizeof(dense_qp_in);

    assert((size_t)c_ptr % 8 == 0 && "memory not 8-byte aligned!");

    d_create_dense_qp(dims, qp_in, c_ptr);
    c_ptr += d_memsize_dense_qp(dims);

    qp_in->dim = (dense_qp_dims *) c_ptr;
    c_ptr += sizeof(dense_qp_dims);

    qp_in->dim->nv = dims->nv;
    qp_in->dim->ne = dims->ne;
    qp_in->dim->nb = dims->nb;
    qp_in->dim->ng = dims->ng;
    qp_in->dim->ns = dims->ns;

    assert((char*)raw_memory + dense_qp_in_calculate_size(dims) == c_ptr);

    return qp_in;
}



int dense_qp_out_calculate_size(dense_qp_dims *dims)
{
    int size = sizeof(dense_qp_out);
    size += d_memsize_dense_qp_sol(dims);

    return size;
}



dense_qp_out *assign_dense_qp_out(dense_qp_dims *dims, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    dense_qp_out *qp_out = (dense_qp_out *) c_ptr;
    c_ptr += sizeof(dense_qp_out);

    assert((size_t)c_ptr % 8 == 0 && "memory not 8-byte aligned!");

    d_create_dense_qp_sol(dims, qp_out, c_ptr);
    c_ptr += d_memsize_dense_qp_sol(dims);

    assert((char*)raw_memory + dense_qp_out_calculate_size(dims) == c_ptr);

    return qp_out;
}
