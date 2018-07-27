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

int ocp_nlp_reg_mirror_memory_calculate_size(ocp_nlp_reg_dims *dims)
{
    int size = 0;

    size += sizeof(ocp_nlp_reg_mirror_memory);

    size += (dims->nx[0]+dims->nu[0])*(dims->nx[0]+dims->nu[0])*sizeof(double);

    return size;
}

void *ocp_nlp_reg_mirror_memory_assign(ocp_nlp_reg_dims *dims, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    ocp_nlp_reg_mirror_memory *mem = (ocp_nlp_reg_mirror_memory *) c_ptr;
    c_ptr += sizeof(ocp_nlp_reg_mirror_memory);

    mem->reg_hess = (double *) c_ptr;
    c_ptr += (dims->nx[0]+dims->nu[0])*(dims->nx[0]+dims->nu[0])*sizeof(double);

    assert((char *) mem + ocp_nlp_reg_mirror_memory_calculate_size(dims) >= c_ptr);

    return mem;
}

void ocp_nlp_reg_mirror(void *config, ocp_nlp_reg_dims *dims, ocp_nlp_reg_in *in,
                        ocp_nlp_reg_out *out, ocp_nlp_reg_opts *opts, void *mem_)
{
    ocp_nlp_reg_mirror_memory *mem = (ocp_nlp_reg_mirror_memory *) mem_;

    int nx = dims->nx[0], nu = dims->nu[0];
    for (int i = 0; i <= dims->N; ++i)
    {
        // make symmetric
        blasfeo_dtrtr_l(nx+nu, &in->RSQrq[i], 0, 0, &in->RSQrq[i], 0, 0);

        // regularize
        blasfeo_unpack_dmat(nx+nu, nx+nu, &in->RSQrq[i], 0, 0, mem->reg_hess, nx+nu);
        project(nx+nu, mem->reg_hess, 1e-4);
        blasfeo_pack_dmat(nx+nu, nx+nu, mem->reg_hess, nx+nu, &in->RSQrq[i], 0, 0);
    }
}

void ocp_nlp_reg_mirror_config_initialize_default(ocp_nlp_reg_config *config)
{
    config->opts_calculate_size = &ocp_nlp_reg_opts_calculate_size;
    config->opts_assign = &ocp_nlp_reg_opts_assign;
    config->memory_calculate_size = &ocp_nlp_reg_mirror_memory_calculate_size;
    config->memory_assign = &ocp_nlp_reg_mirror_memory_assign;

    config->evaluate = &ocp_nlp_reg_mirror;
}
