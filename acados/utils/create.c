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

#if defined (EXT_DEPS)

#include "acados/ocp_nlp/ocp_nlp_common.h"
#include "acados/utils/mem.h"

ocp_nlp_in *create_ocp_nlp_in(ocp_nlp_dims *dims, ocp_nlp_args *args, int num_stages)
{
    int size = ocp_nlp_in_calculate_size(dims, args);
    void *ptr = acados_malloc(size, 1);
    ocp_nlp_in *nlp_in = ocp_nlp_in_assign(dims, args, num_stages, ptr);
    return nlp_in;
}

#endif  // EXT_DEPS