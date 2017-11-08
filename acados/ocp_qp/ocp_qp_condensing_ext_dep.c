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
#include <stdlib.h>
#include <assert.h>
// hpipm
#include "hpipm/include/hpipm_d_dense_qp.h"
#include "hpipm/include/hpipm_d_dense_qp_sol.h"
#include "hpipm/include/hpipm_d_dense_qp_ipm.h"
#include "hpipm/include/hpipm_d_ocp_qp.h"
#include "hpipm/include/hpipm_d_ocp_qp_sol.h"
#include "hpipm/include/hpipm_d_cond.h"
// acados
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/dense_qp/dense_qp_common.h"
#include "acados/ocp_qp/ocp_qp_condensing_ext_dep.h"
#include "acados/ocp_qp/ocp_qp_condensing.h"
#include "acados/utils/types.h"


ocp_qp_condensing_args *ocp_qp_condensing_create_arguments(ocp_qp_dims *dims) {

    int memory_size = ocp_qp_condensing_calculate_args_size(dims);
    void *ptr = malloc(memory_size);
    ocp_qp_condensing_args *args = ocp_qp_condensing_assign_args(dims, ptr);
    return args;
}



ocp_qp_condensing_memory *ocp_qp_condensing_create_memory(ocp_qp_dims *dims, dense_qp_in *out) {

    ocp_qp_condensing_memory *memory;

    int size = ocp_qp_condensing_calculate_memory_size(dims, NULL);
    void *ptr = malloc(size);
    char *ptr_end = assign_ocp_qp_condensing_memory(dims, &memory, ptr);
    assert((char*)ptr + size >= ptr_end); (void) ptr_end;

    return memory;
}
