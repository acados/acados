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
#if defined(RUNTIME_CHECKS)
#include <assert.h>
#endif
// hpipm
#include "hpipm/include/hpipm_d_ocp_qp.h"
#include "hpipm/include/hpipm_d_ocp_qp_sol.h"
#include "hpipm/include/hpipm_d_ocp_qp_ipm.h"
// acados
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/ocp_qp/ocp_qp_hpipm.h"
#include "acados/utils/types.h"
#include "acados/utils/mem.h"


int ocp_qp_hpipm_calculate_args_size(ocp_qp_dims *dims) {

    int size = 0;
    size += sizeof(ocp_qp_hpipm_args);
    size += sizeof(struct d_ocp_qp_ipm_arg);
    size += d_memsize_ocp_qp_ipm_arg(dims);
    size += 1 * 8; // alignment to double word
    make_int_multiple_of(8, &size);
    return size;
}



ocp_qp_hpipm_args *ocp_qp_hpipm_assign_args(ocp_qp_dims *dims, void *mem) {

    ocp_qp_hpipm_args *args;

    char *c_ptr = (char *) mem;

    args = (ocp_qp_hpipm_args *) c_ptr;
    c_ptr += sizeof(ocp_qp_hpipm_args);

    args->hpipm_args = (struct d_ocp_qp_ipm_arg *) c_ptr;
    c_ptr += sizeof(struct d_ocp_qp_ipm_arg);

    align_char_to(8, &c_ptr);
    d_create_ocp_qp_ipm_arg(dims, args->hpipm_args, c_ptr);
    c_ptr += d_memsize_ocp_qp_ipm_arg(dims);

#if defined(RUNTIME_CHECKS)
    assert((char*)mem + ocp_qp_hpipm_calculate_args_size(dims) >= c_ptr);
#endif

    return args;
}



void ocp_qp_hpipm_initialize_default_args(ocp_qp_hpipm_args *args) {

    d_set_default_ocp_qp_ipm_arg(args->hpipm_args);
	// overwrite some default options
    args->hpipm_args->res_g_max = 1e-6;
    args->hpipm_args->res_b_max = 1e-8;
    args->hpipm_args->res_d_max = 1e-8;
    args->hpipm_args->res_m_max = 1e-8;
    args->hpipm_args->iter_max = 50;
    args->hpipm_args->stat_max = 50;
    args->hpipm_args->alpha_min = 1e-8;
    args->hpipm_args->mu0 = 1;
}



int ocp_qp_hpipm_calculate_memory_size(ocp_qp_dims *dims, ocp_qp_hpipm_args *args) {

    int size = 0;
    size += sizeof(ocp_qp_hpipm_memory);

    size += sizeof(struct d_ocp_qp_ipm_workspace);  // ipm_workspace

    size += d_memsize_ocp_qp_ipm(dims, args->hpipm_args);

    make_int_multiple_of(8, &size);
    size += 1 * 8;

    return size;
}



char *ocp_qp_hpipm_assign_memory(ocp_qp_dims *dims, ocp_qp_hpipm_args *args,
        void **mem_, void *raw_memory) {

    ocp_qp_hpipm_memory **hpipm_memory = (ocp_qp_hpipm_memory **) mem_;

    // char pointer
    char *c_ptr = (char *)raw_memory;

    *hpipm_memory = (ocp_qp_hpipm_memory *) c_ptr;
    c_ptr += sizeof(ocp_qp_hpipm_memory);

    //
    (*hpipm_memory)->hpipm_workspace = (struct d_ocp_qp_ipm_workspace *)c_ptr;
    c_ptr += sizeof(struct d_ocp_qp_ipm_workspace);

    struct d_ocp_qp_ipm_workspace *ipm_workspace = (*hpipm_memory)->hpipm_workspace;

    // ipm workspace structure
    align_char_to(8, &c_ptr);
    d_create_ocp_qp_ipm(dims, args->hpipm_args, ipm_workspace, c_ptr);
    c_ptr += ipm_workspace->memsize;

    return c_ptr;
}



int ocp_qp_hpipm(ocp_qp_in *qp_in, ocp_qp_out *qp_out, void *args_, void *mem_) {

    ocp_qp_hpipm_args *args = (ocp_qp_hpipm_args *) args_;
    ocp_qp_hpipm_memory *memory = (ocp_qp_hpipm_memory *) mem_;

    // initialize return code
    int acados_status = ACADOS_SUCCESS;
    // solve ipm
    int hpipm_status = d_solve_ocp_qp_ipm(qp_in, qp_out, args->hpipm_args, memory->hpipm_workspace);

    // printf("HPIPM iter = %d\n", memory->hpipm_workspace->iter);

    // check max number of iterations
    // TODO(dimitris): check ACADOS_MIN_STEP (not implemented in HPIPM yet)
    if (hpipm_status == 1) acados_status = ACADOS_MAXITER;

    // return
    return acados_status;
}
