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

#include "acados_c/external_function.h"

//external
#include <stdlib.h>



external_function_in *create_external_function_in(external_function_dims *dims)
{
    return NULL;
}



external_function_out *create_external_function_out(external_function_dims *dims)
{
    return NULL;
}



int external_function_calculate_args_size(external_function_config *config, external_function_dims *dims)
{
    return 0;
}



void *external_function_assign_args(external_function_config *config, external_function_dims *dims, void *raw_memory)
{
    return NULL;
}



void *external_function_create_args(external_function_config *config, external_function_dims *dims)
{
    return NULL;
}



void *external_function_copy_args(external_function_config  *config, external_function_dims *dims, void *raw_memory, void *source)
{
    return NULL;
}



int external_function_calculate_size(external_function_config *config, external_function_dims *dims, void *args_)
{
    return 0;
}



external_function *external_function_assign(external_function_config *config, external_function_dims *dims, void *args_, void *raw_memory)
{
    return NULL;
}



external_function *external_function_create(external_function_config *config, external_function_dims *dims, void *args_)
{
    return NULL;
}



int external_function_eval(external_function *solver, external_function_in *qp_in, external_function_out *qp_out)
{
    return 0;
}



int set_external_function_fcn_ptrs(external_function_config *config, external_function_fcn_ptrs *fcn_ptrs)
{
    int return_value = ACADOS_SUCCESS;

    return return_value;
}