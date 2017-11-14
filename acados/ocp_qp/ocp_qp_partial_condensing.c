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
// acados
#include "acados/ocp_qp/ocp_qp_partial_condensing.h"
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/utils/mem.h"
// hpipm
#include "hpipm/include/hpipm_d_dense_qp.h"
#include "hpipm/include/hpipm_d_dense_qp_sol.h"
#include "hpipm/include/hpipm_d_ocp_qp.h"
#include "hpipm/include/hpipm_d_ocp_qp_sol.h"
#include "hpipm/include/hpipm_d_ocp_qp_size.h"
#include "hpipm/include/hpipm_d_cond.h"
#include "hpipm/include/hpipm_d_part_cond.h"

int ocp_qp_partial_condensing_calculate_args_size(ocp_qp_dims *dims)
{
    int size = 0;
    size += sizeof(ocp_qp_partial_condensing_args);
    size += d_memsize_ocp_qp_size(dims->N);  // worst-case size of new QP
    return size;
}



ocp_qp_partial_condensing_args *ocp_qp_partial_condensing_assign_args(ocp_qp_dims *dims, void *mem)
{
    char *c_ptr = (char *) mem;

    ocp_qp_partial_condensing_args *args = (ocp_qp_partial_condensing_args *)c_ptr;
    c_ptr += sizeof(ocp_qp_partial_condensing_args);

    args->new_dims = (ocp_qp_dims *)c_ptr;
    c_ptr += d_memsize_ocp_qp_size(dims->N);

    assert((char*)mem + ocp_qp_partial_condensing_calculate_args_size(dims) == c_ptr);

    return args;
}

// TODO(dimitris): DEFAULT ARGS!!!

int ocp_qp_partial_condensing_calculate_memory_size(ocp_qp_dims *dims, ocp_qp_partial_condensing_args *args)
{
    int size = 0;

    // populate dimensions of new ocp_qp based on N2
    args->new_dims->N = args->N2;
    d_compute_qp_size_ocp2ocp(dims, args->new_dims);

    size += sizeof(ocp_qp_partial_condensing_memory);
    size += sizeof(struct d_cond_qp_ocp2dense_workspace);
    size += d_memsize_cond_qp_ocp2ocp(dims, args->new_dims);

    make_int_multiple_of(8, &size);
    size += 1 * 8;

    return size;
}
