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
#include <string.h>
// hpipm
#include "hpipm/include/hpipm_d_cond.h"
#include "hpipm/include/hpipm_d_dense_qp.h"
#include "hpipm/include/hpipm_d_dense_qp_ipm.h"
#include "hpipm/include/hpipm_d_dense_qp_sol.h"
#include "hpipm/include/hpipm_d_ocp_qp.h"
#include "hpipm/include/hpipm_d_ocp_qp_sol.h"
// acados
#include "acados/dense_qp/dense_qp_common.h"
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/ocp_qp/ocp_qp_full_condensing.h"
#include "acados/utils/mem.h"
#include "acados/utils/types.h"



/************************************************
 * dims
 ************************************************/

void compute_dense_qp_dims(ocp_qp_dims *dims, dense_qp_dims *ddims)
{
    d_compute_qp_dim_ocp2dense(dims, ddims);
}



/************************************************
 * opts
 ************************************************/

int ocp_qp_full_condensing_opts_calculate_size(ocp_qp_dims *dims)
{
    int size = 0;
    size += sizeof(ocp_qp_full_condensing_opts);
    // hpipm opts
    size += sizeof(struct d_cond_qp_ocp2dense_arg);
    size += d_memsize_cond_qp_ocp2dense_arg();
    //
    size += 1 * 8;
    make_int_multiple_of(8, &size);

    return size;
}



void *ocp_qp_full_condensing_opts_assign(ocp_qp_dims *dims, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    // opts
    ocp_qp_full_condensing_opts *opts = (ocp_qp_full_condensing_opts *) c_ptr;
    c_ptr += sizeof(ocp_qp_full_condensing_opts);

    // hpipm_opts
    opts->hpipm_opts = (struct d_cond_qp_ocp2dense_arg *) c_ptr;
    c_ptr += sizeof(struct d_cond_qp_ocp2dense_arg);

    align_char_to(8, &c_ptr);

    // hpipm_opts
    d_create_cond_qp_ocp2dense_arg(opts->hpipm_opts, c_ptr);
    c_ptr += opts->hpipm_opts->memsize;

    assert((char *) raw_memory + ocp_qp_full_condensing_opts_calculate_size(dims) >= c_ptr);

    return opts;
}



void ocp_qp_full_condensing_opts_initialize_default(ocp_qp_dims *dims, void *opts_)
{
    ocp_qp_full_condensing_opts *opts = opts_;

    // condense both Hessian and gradient by default
    opts->cond_hess = 1;
    // expand only primal solution (linear MPC, Gauss-Newton)
    opts->expand_dual_sol = 1;
    // hpipm_opts
    d_set_default_cond_qp_ocp2dense_arg(opts->hpipm_opts);
}



void ocp_qp_full_condensing_opts_update(ocp_qp_dims *dims, void *opts_)
{
    ocp_qp_full_condensing_opts *opts = opts_;

    // hpipm_opts
	d_set_cond_qp_ocp2dense_arg_ric_alg(opts->ric_alg, opts->hpipm_opts);

	return;
}



void ocp_qp_full_condensing_opts_set(void *opts_, const char *field, void* value)
{

    ocp_qp_full_condensing_opts *opts = opts_;

	if(!strcmp(field, "ric_alg"))
	{
		int *tmp_ptr = value;
		opts->ric_alg = *tmp_ptr;
	}
	else if(!strcmp(field, "hess"))
	{
		int *tmp_ptr = value;
		opts->cond_hess = *tmp_ptr;
	}
	else if(!strcmp(field, "dual_sol"))
	{
		int *tmp_ptr = value;
		opts->expand_dual_sol = *tmp_ptr;
	}
	else
	{
		printf("\nerror: field %s not available in ocp_qp_full_condensing_opts_set\n", field);
		exit(1);
	}

	return;

}



/************************************************
 * memory
 ************************************************/

int ocp_qp_full_condensing_memory_calculate_size(ocp_qp_dims *dims, void *opts_)
{
    ocp_qp_full_condensing_opts *opts = opts_;
    int size = 0;

    size += sizeof(ocp_qp_full_condensing_memory);
    size += sizeof(struct d_cond_qp_ocp2dense_workspace);
    size += d_memsize_cond_qp_ocp2dense(dims, opts->hpipm_opts);

    return size;
}



void *ocp_qp_full_condensing_memory_assign(ocp_qp_dims *dims, void *opts_, void *raw_memory)
{
    ocp_qp_full_condensing_opts *opts = opts_;

    char *c_ptr = (char *) raw_memory;

    ocp_qp_full_condensing_memory *mem = (ocp_qp_full_condensing_memory *) c_ptr;
    c_ptr += sizeof(ocp_qp_full_condensing_memory);

    mem->hpipm_workspace = (struct d_cond_qp_ocp2dense_workspace *) c_ptr;
    c_ptr += sizeof(struct d_cond_qp_ocp2dense_workspace);

    assert((size_t) c_ptr % 8 == 0 && "memory not 8-byte aligned!");

    // hpipm workspace
    d_create_cond_qp_ocp2dense(dims, opts->hpipm_opts, mem->hpipm_workspace, c_ptr);
    c_ptr += mem->hpipm_workspace->memsize;

    assert((char *) raw_memory + ocp_qp_full_condensing_memory_calculate_size(dims, opts) >= c_ptr);

    return mem;
}



int ocp_qp_full_condensing_workspace_calculate_size(ocp_qp_dims *dims, void *opts_)
{
return 0;
}



void ocp_qp_full_condensing(ocp_qp_in *in, dense_qp_in *out, ocp_qp_full_condensing_opts *opts,
                            ocp_qp_full_condensing_memory *mem, void *work)
{
    // save pointer to ocp_qp_in in memory (needed for expansion)
    mem->qp_in = in;

    // convert to dense qp structure
    if (opts->cond_hess == 0)
    {
        // condense gradient only
        d_cond_rhs_qp_ocp2dense(in, out, opts->hpipm_opts, mem->hpipm_workspace);
    }
    else
    {
        // condense gradient and Hessian
        d_cond_qp_ocp2dense(in, out, opts->hpipm_opts, mem->hpipm_workspace);

        // ++ for debugging ++
        //
        // printf("gradient with full condensing:\n\n");
        // blasfeo_print_dvec(out->g->m, out->g, 0);

        // d_cond_rhs_qp_ocp2dense(in, out, mem->hpipm_workspace);

        // printf("gradient with gradient-only condensing:\n\n");
        // blasfeo_print_dvec(out->g->m, out->g, 0);
    }
}



void ocp_qp_full_expansion(dense_qp_out *in, ocp_qp_out *out, ocp_qp_full_condensing_opts *opts,
                           ocp_qp_full_condensing_memory *mem, void *work)
{
    if (opts->expand_dual_sol == 0)
    {
        d_expand_primal_sol_dense2ocp(mem->qp_in, in, out, opts->hpipm_opts, mem->hpipm_workspace);
    }
    else
    {
        d_expand_sol_dense2ocp(mem->qp_in, in, out, opts->hpipm_opts, mem->hpipm_workspace);
    }
}



void ocp_qp_full_condensing_config_initialize_default(void *config_)
{
    ocp_qp_condensing_config *config = config_;

    config->opts_calculate_size = &ocp_qp_full_condensing_opts_calculate_size;
    config->opts_assign = &ocp_qp_full_condensing_opts_assign;
    config->opts_initialize_default = &ocp_qp_full_condensing_opts_initialize_default;
    config->opts_update = &ocp_qp_full_condensing_opts_update;
    config->opts_set = &ocp_qp_full_condensing_opts_set;
    config->memory_calculate_size = &ocp_qp_full_condensing_memory_calculate_size;
    config->memory_assign = &ocp_qp_full_condensing_memory_assign;
    config->workspace_calculate_size = &ocp_qp_full_condensing_workspace_calculate_size;
    // TODO(dimitris): either do casting as below or void in defs, as above (also pass config or
    // not?)
    config->condensing = (int (*)(void *, void *, void *, void *, void *)) & ocp_qp_full_condensing;
    config->expansion = (int (*)(void *, void *, void *, void *, void *)) & ocp_qp_full_expansion;

    return;
}
