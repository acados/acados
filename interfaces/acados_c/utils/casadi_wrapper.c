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

#include "acados_c/utils/casadi_wrapper.h"

#include <stdlib.h>



// void *casadi_wrapper_copy_args(external_function_config *config, external_function_dims *dims, void *raw_memory, void *source_)
// {
//     casadi_wrapper_args *source = (casadi_wrapper_args *)source_;
//     casadi_wrapper_args *dest;

//     dest = casadi_wrapper_assign_args(dims, NULL, raw_memory);

//     dest->fun = source->fun;

//     dest->dims = source->dims;

//     dest->sparsity = source->sparsity;
    
//     return (void *)dest;
// }



int casadi_wrapper_calculate_submodules_size(external_function_config *config, external_function_dims *dims)
{
    return 0;
}



void *casadi_wrapper_assign_submodules(external_function_config *config, external_function_dims *dims, void *raw_memory)
{
    return NULL;
}
