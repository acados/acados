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

#include "acados/ocp_nlp/ocp_nlp_reg_common.h"

#include <assert.h>

int ocp_nlp_reg_in_calculate_size(void)
{
    int size = 0;

    size += sizeof(ocp_nlp_reg_in);

    return size;
}

void *ocp_nlp_reg_in_assign(void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    ocp_nlp_reg_in *in = (ocp_nlp_reg_in *) c_ptr;
    c_ptr += sizeof(ocp_nlp_reg_in);

    assert((char *) in + ocp_nlp_reg_in_calculate_size() >= c_ptr);

    return in;
}

int ocp_nlp_reg_out_calculate_size(void)
{
    int size = 0;

    size += sizeof(ocp_nlp_reg_out);

    return size;
}

void *ocp_nlp_reg_out_assign(void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    ocp_nlp_reg_out *out = (ocp_nlp_reg_out *) c_ptr;

    assert((char *) out + ocp_nlp_reg_out_calculate_size() >= c_ptr);

    return out;
}

int ocp_nlp_reg_opts_calculate_size(void)
{
    return sizeof(ocp_nlp_reg_opts);
}

void *ocp_nlp_reg_opts_assign(void *raw_memory)
{
    return raw_memory;
}

int ocp_nlp_reg_config_calculate_size(void)
{
    return sizeof(ocp_nlp_reg_config);
}

void *ocp_nlp_reg_config_assign(void *raw_memory)
{
    return raw_memory;
}
