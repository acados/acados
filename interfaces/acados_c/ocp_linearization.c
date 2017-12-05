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

#include "acados_c/ocp_linearization.h"

//external
#include <stdlib.h>



ocp_linearization_in *create_ocp_linearization_in(ocp_linearization_dims *dims)
{
    return NULL;
}



ocp_linearization_out *create_ocp_linearization_out(ocp_linearization_dims *dims)
{
    return NULL;
}



int ocp_linearization_calculate_size(ocp_linearization_dims *dims, void *args_)
{
    return 0;
}



ocp_linearization_method *ocp_linearization_assign(ocp_linearization_dims *dims, void *args_, void *raw_memory)
{
    return NULL;
}



ocp_linearization_method *ocp_linearization_create(ocp_linearization_dims *dims, void *args_)
{
    return NULL;
}



int ocp_linearization_solve(ocp_linearization_method *solver, ocp_linearization_in *qp_in, ocp_linearization_out *qp_out)
{
    return 0;
}



void ocp_linearization_initialize_default_args(ocp_linearization_method *method)
{

}



int ocp_linearization_calculate_args_size(ocp_linearization_dims *dims, void *args_)
{
    return 0;
}



void *ocp_linearization_assign_args(ocp_linearization_dims *dims, void *args_, void *raw_memory)
{
    return NULL;
}



void *ocp_linearization_create_args(ocp_linearization_dims *dims, void *args_)
{
    return NULL;
}



int ocp_linearization_calculate_memory_size(ocp_linearization_dims *dims, void *args_)
{
    return 0;
}



void *ocp_linearization_assign_memory(ocp_linearization_dims *dims, void *args_, void *raw_memory)
{
    return NULL;
}



void *ocp_linearization_create_memory(ocp_linearization_dims *dims, void *args_)
{
    return NULL;
}



int ocp_linearization_calculate_workspace_size(ocp_linearization_dims *dims, void *args_)
{
    return 0;
}
