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

#include "acados/ocp_nlp/ocp_nlp_cost_common.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>

// blasfeo
#include "blasfeo/include/blasfeo_d_aux.h"
#include "blasfeo/include/blasfeo_d_blas.h"
// acados
#include "acados/utils/mem.h"



/************************************************
 * config
 ************************************************/

int ocp_nlp_cost_config_calculate_size()
{
    int size = 0;

    size += sizeof(ocp_nlp_cost_config);

    return size;
}

ocp_nlp_cost_config *ocp_nlp_cost_config_assign(void *raw_memory)
{
    char *c_ptr = raw_memory;

    ocp_nlp_cost_config *config = (ocp_nlp_cost_config *) c_ptr;
    c_ptr += sizeof(ocp_nlp_cost_config);

    return config;
}
