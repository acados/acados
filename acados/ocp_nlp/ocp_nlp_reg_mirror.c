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

#include "acados/ocp_nlp/ocp_nlp_reg_mirror.h"

#include <assert.h>

#include "acados/ocp_nlp/ocp_nlp_reg_common.h"
#include "acados/utils/math.h"

#include "blasfeo/include/blasfeo_d_aux.h"
#include "blasfeo/include/blasfeo_d_blas.h"



/************************************************
 * memory
 ************************************************/

int ocp_nlp_reg_mirror_memory_calculate_size(ocp_nlp_reg_dims *dims)
{
    int *nx = dims->nx;
	int *nu = dims->nu;
	int N = dims->N;

	int ii;

	int nuxM = nu[0]+nx[0];
	for(ii=1; ii<=N; ii++)
	{
		nuxM = nu[ii]+nx[ii]>nuxM ? nu[ii]+nx[ii] : nuxM;
	}

    int size = 0;

    size += sizeof(ocp_nlp_reg_mirror_memory);

    size += nuxM*nuxM*sizeof(double);  // reg_hess
    size += nuxM*nuxM*sizeof(double);  // V
    size += 2*nuxM*sizeof(double);     // d e
	size += (N+1)*sizeof(struct blasfeo_dmat *); // RSQrq

    return size;
}



void *ocp_nlp_reg_mirror_memory_assign(ocp_nlp_reg_dims *dims, void *raw_memory)
{
    int *nx = dims->nx;
	int *nu = dims->nu;
	int N = dims->N;

	int ii;

	int nuxM = nu[0]+nx[0];
	for(ii=1; ii<=N; ii++)
	{
		nuxM = nu[ii]+nx[ii]>nuxM ? nu[ii]+nx[ii] : nuxM;
	}

    char *c_ptr = (char *) raw_memory;

    ocp_nlp_reg_mirror_memory *mem = (ocp_nlp_reg_mirror_memory *) c_ptr;
    c_ptr += sizeof(ocp_nlp_reg_mirror_memory);

    mem->reg_hess = (double *) c_ptr;
    c_ptr += nuxM*nuxM*sizeof(double);  // reg_hess

    mem->V = (double *) c_ptr;
    c_ptr += nuxM*nuxM*sizeof(double);  // V

    mem->d = (double *) c_ptr;
    c_ptr += nuxM*sizeof(double); // d

    mem->e = (double *) c_ptr;
    c_ptr += nuxM*sizeof(double); // e

	mem->RSQrq = (struct blasfeo_dmat **) c_ptr;
	c_ptr += (N+1)*sizeof(struct blasfeo_dmat *); // RSQrq

    assert((char *) mem + ocp_nlp_reg_mirror_memory_calculate_size(dims) >= c_ptr);

    return mem;
}



void ocp_nlp_reg_mirror_memory_set_RSQrq_ptr(int N, struct blasfeo_dmat *RSQrq, void *memory_)
{
    ocp_nlp_reg_mirror_memory *memory = memory_;

	int ii;

	for(ii=0; ii<=N; ii++)
	{
		memory->RSQrq[ii] = RSQrq+ii;
	}

    return;
}



/************************************************
 * functions
 ************************************************/

void ocp_nlp_reg_mirror(void *config, ocp_nlp_reg_dims *dims, ocp_nlp_reg_in *in,
                        ocp_nlp_reg_out *out, ocp_nlp_reg_opts *opts, void *mem_)
{
    ocp_nlp_reg_mirror_memory *mem = (ocp_nlp_reg_mirror_memory *) mem_;

	int ii;

    int *nx = dims->nx;
	int *nu = dims->nu;
	int N = dims->N;

    for(ii=0; ii<=dims->N; ii++)
    {
        // make symmetric
        blasfeo_dtrtr_l(nu[ii]+nx[ii], &in->RSQrq[ii], 0, 0, &in->RSQrq[ii], 0, 0);

        // regularize
        blasfeo_unpack_dmat(nu[ii]+nx[ii], nu[ii]+nx[ii], &in->RSQrq[ii], 0, 0, mem->reg_hess, nu[ii]+nx[ii]);
//        acados_mirror(nx+nu, mem->reg_hess, mem->V, mem->d, mem->e, 1e-4);
        acados_project(nu[ii]+nx[ii], mem->reg_hess, mem->V, mem->d, mem->e, 1e-4);
        blasfeo_pack_dmat(nu[ii]+nx[ii], nu[ii]+nx[ii], mem->reg_hess, nu[ii]+nx[ii], &in->RSQrq[ii], 0, 0);
    }
}



void ocp_nlp_reg_mirror_config_initialize_default(ocp_nlp_reg_config *config)
{
	// dims
    config->dims_calculate_size = &ocp_nlp_reg_dims_calculate_size;
    config->dims_assign = &ocp_nlp_reg_dims_assign;
    config->dims_set = &ocp_nlp_reg_dims_set;
	// opts
    config->opts_calculate_size = &ocp_nlp_reg_opts_calculate_size;
    config->opts_assign = &ocp_nlp_reg_opts_assign;
    config->opts_initialize_default = &ocp_nlp_reg_opts_initialize_default;
    config->opts_set = &ocp_nlp_reg_opts_set;
	// memory
    config->memory_calculate_size = &ocp_nlp_reg_mirror_memory_calculate_size;
    config->memory_assign = &ocp_nlp_reg_mirror_memory_assign;
    config->memory_set_RSQrq_ptr = &ocp_nlp_reg_mirror_memory_set_RSQrq_ptr;
	// functions
    config->evaluate = &ocp_nlp_reg_mirror;
}
