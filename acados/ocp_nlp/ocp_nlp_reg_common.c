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
