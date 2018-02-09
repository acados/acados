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

#include "acados_c/dense_qp/dense_qp_qpoases.h"



// void *dense_qp_qpoases_copy_args(dense_qp_solver_config *config, dense_qp_dims *dims, void *raw_memory, void *source_)
// {
//     dense_qp_qpoases_args *source = (dense_qp_qpoases_args *) source_;
//     dense_qp_qpoases_args *dest;

//     dest = (dense_qp_qpoases_args *) dense_qp_qpoases_assign_args(dims, raw_memory);

//     dest->max_cputime = source->max_cputime;
//     dest->warm_start = source->warm_start;
//     dest->max_nwsr = source->max_nwsr;
//     dest->use_precomputed_cholesky = source->use_precomputed_cholesky;
//     dest->hotstart = source->hotstart;

//     return (void *)dest;
// }



int dense_qp_qpoases_calculate_submodules_size(dense_qp_solver_config *config, dense_qp_dims *dims)
{
    return 0;
}



void *dense_qp_qpoases_assign_submodules(dense_qp_solver_config *config, dense_qp_dims *dims, void *raw_memory)
{
    return NULL;
}
