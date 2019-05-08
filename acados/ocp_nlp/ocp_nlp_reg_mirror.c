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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "acados/ocp_nlp/ocp_nlp_reg_common.h"
#include "acados/utils/math.h"

#include "blasfeo/include/blasfeo_d_aux.h"
#include "blasfeo/include/blasfeo_d_blas.h"



/************************************************
 * opts
 ************************************************/

int ocp_nlp_reg_mirror_opts_calculate_size(void)
{
    return sizeof(ocp_nlp_reg_mirror_opts);
}



void *ocp_nlp_reg_mirror_opts_assign(void *raw_memory)
{
    return raw_memory;
}



void ocp_nlp_reg_mirror_opts_initialize_default(void *config_, ocp_nlp_reg_dims *dims, void *opts_)
{
    ocp_nlp_reg_mirror_opts *opts = opts_;

    opts->epsilon = 1e-4;

    return;
}



void ocp_nlp_reg_mirror_opts_set(void *config_, ocp_nlp_reg_dims *dims, void *opts_, char *field, void* value)
{

    ocp_nlp_reg_mirror_opts *opts = opts_;

    if (!strcmp(field, "epsilon"))
    {
        double *d_ptr = value;
        opts->epsilon = *d_ptr;
    }
    else
    {
        printf("\nerror: field %s not available in ocp_nlp_reg_mirror_opts_set\n", field);
        exit(1);
    }

    return;
}



/************************************************
 * memory
 ************************************************/

int ocp_nlp_reg_mirror_memory_calculate_size(void *config_, ocp_nlp_reg_dims *dims, void *opts_)
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



void *ocp_nlp_reg_mirror_memory_assign(void *config_, ocp_nlp_reg_dims *dims, void *opts_, void *raw_memory)
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

    assert((char *) mem + ocp_nlp_reg_mirror_memory_calculate_size(config_, dims, opts_) >= c_ptr);

    return mem;
}



void ocp_nlp_reg_mirror_memory_set_RSQrq_ptr(ocp_nlp_reg_dims *dims, struct blasfeo_dmat *RSQrq, void *memory_)
{
    ocp_nlp_reg_mirror_memory *memory = memory_;

    int ii;

    int N = dims->N;
    int *nx = dims->nx;
    int *nu = dims->nu;

    for(ii=0; ii<=N; ii++)
    {
        memory->RSQrq[ii] = RSQrq+ii;
//        blasfeo_print_dmat(nu[ii]+nx[ii]+1, nu[ii]+nx[ii], memory->RSQrq[ii], 0, 0);
    }

    return;
}



void ocp_nlp_reg_mirror_memory_set_rq_ptr(ocp_nlp_reg_dims *dims, struct blasfeo_dvec *rq, void *memory_)
{
    return;
}



void ocp_nlp_reg_mirror_memory_set_BAbt_ptr(ocp_nlp_reg_dims *dims, struct blasfeo_dmat *BAbt, void *memory_)
{
    return;
}



void ocp_nlp_reg_mirror_memory_set_b_ptr(ocp_nlp_reg_dims *dims, struct blasfeo_dvec *b, void *memory_)
{
    return;
}



void ocp_nlp_reg_mirror_memory_set_ux_ptr(ocp_nlp_reg_dims *dims, struct blasfeo_dvec *ux, void *memory_)
{
    return;
}



void ocp_nlp_reg_mirror_memory_set_pi_ptr(ocp_nlp_reg_dims *dims, struct blasfeo_dvec *pi, void *memory_)
{
    return;
}



void ocp_nlp_reg_mirror_memory_set(void *config_, ocp_nlp_reg_dims *dims, void *memory_, char *field, void *value)
{

    if(!strcmp(field, "RSQrq_ptr"))
    {
        struct blasfeo_dmat *RSQrq = value;
        ocp_nlp_reg_mirror_memory_set_RSQrq_ptr(dims, RSQrq, memory_);
    }
    else
    {
        printf("\nerror: field %s not available in ocp_nlp_reg_mirror_set\n", field);
        exit(1);
    }

    return;
}



/************************************************
 * functions
 ************************************************/

void ocp_nlp_reg_mirror_regularize_hessian(void *config, ocp_nlp_reg_dims *dims, void *opts_, void *mem_)
{
    ocp_nlp_reg_mirror_memory *mem = (ocp_nlp_reg_mirror_memory *) mem_;
    ocp_nlp_reg_mirror_opts *opts = opts_;

    int ii;

    int *nx = dims->nx;
    int *nu = dims->nu;
    int N = dims->N;

    for(ii=0; ii<=dims->N; ii++)
    {
        // make symmetric
        blasfeo_dtrtr_l(nu[ii]+nx[ii], mem->RSQrq[ii], 0, 0, mem->RSQrq[ii], 0, 0);

        // regularize
        blasfeo_unpack_dmat(nu[ii]+nx[ii], nu[ii]+nx[ii], mem->RSQrq[ii], 0, 0, mem->reg_hess, nu[ii]+nx[ii]);
        acados_mirror(nu[ii]+nx[ii], mem->reg_hess, mem->V, mem->d, mem->e, opts->epsilon);
        blasfeo_pack_dmat(nu[ii]+nx[ii], nu[ii]+nx[ii], mem->reg_hess, nu[ii]+nx[ii], mem->RSQrq[ii], 0, 0);
    }
}



void ocp_nlp_reg_mirror_correct_dual_sol(void *config, ocp_nlp_reg_dims *dims, void *opts_, void *mem_)
{
    return;
}



void ocp_nlp_reg_mirror_config_initialize_default(ocp_nlp_reg_config *config)
{
    // dims
    config->dims_calculate_size = &ocp_nlp_reg_dims_calculate_size;
    config->dims_assign = &ocp_nlp_reg_dims_assign;
    config->dims_set = &ocp_nlp_reg_dims_set;
    // opts
    config->opts_calculate_size = &ocp_nlp_reg_mirror_opts_calculate_size;
    config->opts_assign = &ocp_nlp_reg_mirror_opts_assign;
    config->opts_initialize_default = &ocp_nlp_reg_mirror_opts_initialize_default;
    config->opts_set = &ocp_nlp_reg_mirror_opts_set;
    // memory
    config->memory_calculate_size = &ocp_nlp_reg_mirror_memory_calculate_size;
    config->memory_assign = &ocp_nlp_reg_mirror_memory_assign;
    config->memory_set = &ocp_nlp_reg_mirror_memory_set;
    config->memory_set_RSQrq_ptr = &ocp_nlp_reg_mirror_memory_set_RSQrq_ptr;
    config->memory_set_rq_ptr = &ocp_nlp_reg_mirror_memory_set_rq_ptr;
    config->memory_set_BAbt_ptr = &ocp_nlp_reg_mirror_memory_set_BAbt_ptr;
    config->memory_set_b_ptr = &ocp_nlp_reg_mirror_memory_set_b_ptr;
    config->memory_set_ux_ptr = &ocp_nlp_reg_mirror_memory_set_ux_ptr;
    config->memory_set_pi_ptr = &ocp_nlp_reg_mirror_memory_set_pi_ptr;
    // functions
    config->regularize_hessian = &ocp_nlp_reg_mirror_regularize_hessian;
    config->correct_dual_sol = &ocp_nlp_reg_mirror_correct_dual_sol;
}
