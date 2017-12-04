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

#include "acados_c/dense_qp.h"

#include <acados/dense_qp/dense_qp_common.h>
#include <acados/dense_qp/dense_qp_hpipm.h>
#include <acados/dense_qp/dense_qp_qore.h>
#include <acados/dense_qp/dense_qp_qpoases.h>



int dense_qp_calculate_size(dense_qp_config *config, dense_qp_dims *dims)
{
    return 0;
}



dense_qp_solver *dense_qp_assign(dense_qp_config *config, dense_qp_dims *dims, void *raw_memory)
{
    return NULL;
}



dense_qp_solver *dense_qp_create(dense_qp_config *config, dense_qp_dims *dims)
{
    return NULL;
}



int dense_qp_solve(dense_qp_solver *solver, dense_qp_in *qp_in, dense_qp_out *qp_out)
{
    return 0;
}



int dense_qp_calculate_args_size(dense_qp_config *config, dense_qp_dims *dims)
{
    return 0;
}



void *dense_qp_assign_args(dense_qp_config *config, dense_qp_dims *dims, void *raw_memory)
{
    return NULL;
}



void *dense_qp_create_args(dense_qp_config *config, dense_qp_dims *dims)
{
    return NULL;
}



void dense_qp_assign_default_args(dense_qp_config *config, void *args_)
{

}



int dense_qp_calculate_memory_size(dense_qp_dims *dims, void *args_)
{
    return 0;
}



void *dense_qp_assign_memory(dense_qp_dims *dims, void *args_, void *raw_memory)
{
    return NULL;
}



void *dense_qp_create_memory(dense_qp_dims *dims, void *args_)
{
    return NULL;
}



int dense_qp_calculate_workspace_size(dense_qp_dims *dims, void *args_)
{
    return 0;
}
