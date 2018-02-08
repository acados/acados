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

#include "acados_c/dense_qp/dense_qp_hpipm.h"



// void *dense_qp_hpipm_copy_args(dense_qp_solver_config *config, dense_qp_dims *dims, void *raw_memory, void *source_)
// {
//     dense_qp_hpipm_args *source = (dense_qp_hpipm_args *)source_;
//     dense_qp_hpipm_args *dest;

//     dest = (dense_qp_hpipm_args *) dense_qp_assign_args(config, dims, raw_memory);

//     *dest->hpipm_args = *source->hpipm_args;

//     return (void *)dest;
// }



int dense_qp_hpipm_calculate_submodules_size(dense_qp_solver_config *config, dense_qp_dims *dims)
{
    return 0;
}



void *dense_qp_hpipm_assign_submodules(dense_qp_solver_config *config, dense_qp_dims *dims, void *raw_memory)
{
    return NULL;
}
