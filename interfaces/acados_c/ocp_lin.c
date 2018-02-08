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

#include "acados_c/ocp_lin.h"

//external
#include <stdlib.h>



ocp_lin_dims *create_ocp_lin_dims()
{
    return NULL;
}



ocp_lin_in *create_ocp_lin_in(ocp_lin_dims *dims)
{
    return NULL;
}



ocp_lin_out *create_ocp_lin_out(ocp_lin_dims *dims)
{
    return NULL;
}



int ocp_lin_calculate_args_size(ocp_lin_method_fcn_ptrs *fcn_ptrs, ocp_lin_dims *dims)
{
    return 0;
}



void *ocp_lin_assign_args(ocp_lin_method_fcn_ptrs *fcn_ptrs, ocp_lin_dims *dims, void *raw_memory)
{
    return NULL;
}



void *ocp_lin_create_args(ocp_lin_method_fcn_ptrs *fcn_ptrs, ocp_lin_dims *dims)
{
    return NULL;
}



void *ocp_lin_copy_args(ocp_lin_method_fcn_ptrs *fcn_ptrs, ocp_lin_dims *dims, void *raw_memory, void *source)
{
    return NULL;
}



int ocp_lin_calculate_size(ocp_lin_method_fcn_ptrs *fcn_ptrs, ocp_lin_dims *dims, void *args_)
{
    return 0;
}



ocp_lin_method *ocp_lin_assign(ocp_lin_method_fcn_ptrs *fcn_ptrs, ocp_lin_dims *dims, void *args_, void *raw_memory)
{
    return NULL;
}



ocp_lin_method *ocp_lin_create(ocp_lin_method_fcn_ptrs *fcn_ptrs, ocp_lin_dims *dims, void *args_)
{
    return NULL;
}



int ocp_lin_solve(ocp_lin_method *solver, ocp_lin_in *qp_in, ocp_lin_out *qp_out)
{
    return 0;
}



int ocp_lin_calculate_submodules_size(ocp_lin_method_config *config, ocp_lin_dims *dims)
{
    return 0;
}



void *ocp_lin_assign_submodules(ocp_lin_method_config *config, ocp_lin_dims *dims, void *raw_memory)
{
    return NULL;
}



int calculate_ocp_lin_solver_fcn_ptrs_size(ocp_lin_method_config *config, ocp_lin_dims *dims)
{
    return 0;
}



void *assign_ocp_lin_solver_fcn_ptrs(ocp_lin_method_config *config, ocp_lin_dims *dims, void *raw_memory)
{
    return NULL;
}



void *create_ocp_lin_solver_fcn_ptrs(ocp_lin_method_config *config, ocp_lin_dims *dims)
{
    return NULL;
}



int set_ocp_lin_method_fcn_ptrs(ocp_lin_method_config *config, ocp_lin_method_fcn_ptrs *fcn_ptrs)
{
    int return_value = ACADOS_SUCCESS;

    return return_value;
}
