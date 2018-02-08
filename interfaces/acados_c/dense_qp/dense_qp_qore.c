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

#include "acados_c/dense_qp/dense_qp_qore.h"



// void *dense_qp_qore_copy_args(dense_qp_solver_config *config, dense_qp_dims *dims, void *raw_memory, void *source_)
// {
//     dense_qp_qore_args *source = (dense_qp_qore_args *)source_;
//     dense_qp_qore_args *dest;

//     dest = (dense_qp_qore_args *) dense_qp_qore_assign_args(dims, raw_memory);

//     dest->print_freq = source->print_freq;
//     dest->warm_start = source->warm_start;
//     dest->warm_strategy = source->warm_strategy;
//     dest->nsmax = source->nsmax;
//     dest->hot_start = source->hot_start;
//     dest->max_iter = source->max_iter;

//     return (void *)dest;
// }



int dense_qp_qore_calculate_submodules_size(dense_qp_solver_config *config, dense_qp_dims *dims)
{
    return 0;
}



void *dense_qp_qore_assign_submodules(dense_qp_solver_config *config, dense_qp_dims *dims, void *raw_memory)
{
    return NULL;
}
