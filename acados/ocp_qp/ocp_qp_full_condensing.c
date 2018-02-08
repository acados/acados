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
#include "hpipm/include/hpipm_d_dense_qp.h"
#include "hpipm/include/hpipm_d_dense_qp_sol.h"
#include "hpipm/include/hpipm_d_dense_qp_ipm.h"
#include "hpipm/include/hpipm_d_ocp_qp.h"
#include "hpipm/include/hpipm_d_ocp_qp_sol.h"
#include "hpipm/include/hpipm_d_cond.h"
// acados
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/dense_qp/dense_qp_common.h"
#include "acados/ocp_qp/ocp_qp_full_condensing.h"
#include "acados/utils/types.h"
#include "acados/utils/mem.h"


void compute_dense_qp_dims(ocp_qp_dims *dims, dense_qp_dims *ddims)
{
    d_compute_qp_dim_ocp2dense(dims, ddims);
}



int ocp_qp_full_condensing_calculate_args_size(ocp_qp_dims *dims, void *submodules_)
{
    int size = 0;
    size += sizeof(ocp_qp_full_condensing_args);
    return size;
}



ocp_qp_full_condensing_args *ocp_qp_full_condensing_assign_args(ocp_qp_dims *dims, void **submodules_, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    ocp_qp_full_condensing_args *args = (ocp_qp_full_condensing_args *) c_ptr;
    c_ptr += sizeof(ocp_qp_full_condensing_args);

    assert((char*)raw_memory + ocp_qp_full_condensing_calculate_args_size(dims, *submodules_) == c_ptr);

    // Update submodules pointer
    *submodules_ = NULL;

    return args;
}



ocp_qp_full_condensing_args *ocp_qp_full_condensing_copy_args(ocp_qp_dims *dims, void *raw_memory, ocp_qp_full_condensing_args *source_)
{
    ocp_qp_full_condensing_args *dest;

    void *submodules;

    dest = ocp_qp_full_condensing_assign_args(dims, &submodules, raw_memory);

    return dest;
}



void ocp_qp_full_condensing_initialize_default_args(ocp_qp_full_condensing_args *args) {
	
	// condense both Hessian and gradient by default
	args->condense_rhs_only = 0;
	// expand only primal solution (linear MPC, Gauss-Newton)
	args->expand_primal_sol_only = 0;
}

int ocp_qp_full_condensing_calculate_memory_size(ocp_qp_dims *dims, ocp_qp_full_condensing_args *args)
{
    int size = 0;

    size += sizeof(ocp_qp_full_condensing_memory);
    size += sizeof(struct d_cond_qp_ocp2dense_workspace);
    size += d_memsize_cond_qp_ocp2dense(dims);

    return size;
}



ocp_qp_full_condensing_memory *ocp_qp_full_condensing_assign_memory(ocp_qp_dims *dims,
    ocp_qp_full_condensing_args *args, void *raw_memory)
{
    char *c_ptr = (char *)raw_memory;

    ocp_qp_full_condensing_memory *mem = (ocp_qp_full_condensing_memory *) c_ptr;
    c_ptr += sizeof(ocp_qp_full_condensing_memory);

    mem->hpipm_workspace = (struct d_cond_qp_ocp2dense_workspace *)c_ptr;
    c_ptr += sizeof(struct d_cond_qp_ocp2dense_workspace);

    assert((size_t)c_ptr % 8 == 0 && "memory not 8-byte aligned!");

    // hpipm workspace structure
    d_create_cond_qp_ocp2dense(dims, mem->hpipm_workspace, c_ptr);
    c_ptr += mem->hpipm_workspace->memsize;

    assert((char*)raw_memory + ocp_qp_full_condensing_calculate_memory_size(dims, args) == c_ptr);

    return mem;
}



int ocp_qp_full_condensing_calculate_workspace_size(ocp_qp_dims *dims, ocp_qp_full_condensing_args *args)
{
    return 0;
}



void ocp_qp_full_condensing(ocp_qp_in *in, dense_qp_in *out, ocp_qp_full_condensing_args *args,
    ocp_qp_full_condensing_memory *mem, void *work)
{
    // save pointer to ocp_qp_in in memory (needed for expansion)
    mem->qp_in = in;

    // convert to dense qp structure
	if(args->condense_rhs_only == 1) {
		// condense gradient only
		d_cond_rhs_qp_ocp2dense(in, out, mem->hpipm_workspace);
	} else {
		// condense gradient and Hessian
		d_cond_qp_ocp2dense(in, out, mem->hpipm_workspace);

		// ++ for debugging ++
		//
		// printf("gradient with full condensing:\n\n"); 
		// blasfeo_print_dvec(out->g->m, out->g, 0);

		// d_cond_rhs_qp_ocp2dense(in, out, mem->hpipm_workspace);

		// printf("gradient with gradient-only condensing:\n\n");	
		// blasfeo_print_dvec(out->g->m, out->g, 0);
	}
}



void ocp_qp_full_expansion(dense_qp_out *in, ocp_qp_out *out, ocp_qp_full_condensing_args *args,
    ocp_qp_full_condensing_memory *mem, void *work)
{
	if(args->expand_primal_sol_only == 1) {
		d_expand_primal_sol_dense2ocp(mem->qp_in, in, out, mem->hpipm_workspace);
	} else {
		d_expand_sol_dense2ocp(mem->qp_in, in, out, mem->hpipm_workspace);

	}
}
