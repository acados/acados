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

#include "acados_c/ocp_qp/ocp_qp_full_condensing.h"



// ocp_qp_full_condensing_args *ocp_qp_full_condensing_copy_args(ocp_qp_solver_config *config, ocp_qp_dims *dims, void *raw_memory, ocp_qp_full_condensing_args *source_)
// {
//     ocp_qp_full_condensing_args *dest;

//     dest = ocp_qp_assign_args(config, dims, raw_memory);

//     return dest;
// }



int ocp_qp_full_condensing_calculate_submodules_size(ocp_qp_solver_config *config, ocp_qp_dims *dims)
{
    return 0;
}



void *ocp_qp_full_condensing_assign_submodules(ocp_qp_solver_config *config, ocp_qp_dims *dims, void *raw_memory)
{
    return NULL;
}
