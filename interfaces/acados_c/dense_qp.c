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

//external
#include <stdlib.h>
#include <assert.h>
#include <string.h>
//acados
#include <acados/dense_qp/dense_qp_common.h>
//acados_c
#include "acados_c/dense_qp/dense_qp_hpipm.h"
#ifdef ACADOS_WITH_QORE
#include "acados_c/dense_qp/dense_qp_qore.h"
#endif
#include "acados_c/dense_qp/dense_qp_qpoases.h"



void dense_qp_copy_dims(dense_qp_dims *dest, dense_qp_dims *src)
{
	dest->nv = src->nv;
	dest->ne = src->ne;
	dest->nb = src->nb;
	dest->ng = src->ng;
	dest->ns = src->ns;
	dest->memsize = src->memsize;
}



dense_qp_dims *create_dense_qp_dims()
{
    int bytes = dense_qp_dims_calculate_size();

    void *ptr = malloc(bytes);

    dense_qp_dims *dims = assign_dense_qp_dims(ptr);

    return dims;
}



dense_qp_in *create_dense_qp_in(dense_qp_dims *dims)
{
    int bytes = dense_qp_in_calculate_size(dims);

    void *ptr = malloc(bytes);

    dense_qp_in *in = assign_dense_qp_in(dims, ptr);

    return in;
}



dense_qp_out *create_dense_qp_out(dense_qp_dims *dims)
{
    int bytes = dense_qp_out_calculate_size(dims);

    void *ptr = malloc(bytes);

    dense_qp_out *out = assign_dense_qp_out(dims, ptr);

    return out;
}



int dense_qp_calculate_args_size(dense_qp_solver_fcn_ptrs *fcn_ptrs, dense_qp_dims *dims)
{
    return fcn_ptrs->calculate_args_size(dims, fcn_ptrs->submodules);
}



void *dense_qp_assign_args(dense_qp_solver_fcn_ptrs *fcn_ptrs, dense_qp_dims *dims, void *raw_memory)
{
    void *submodules = fcn_ptrs->submodules;
    void *args = fcn_ptrs->assign_args(dims, &submodules, raw_memory);

    fcn_ptrs->initialize_default_args(args);

    return args;
}



void *dense_qp_create_args(dense_qp_solver_fcn_ptrs *fcn_ptrs, dense_qp_dims *dims)
{
    int bytes = dense_qp_calculate_args_size(fcn_ptrs, dims);

    void *ptr = malloc(bytes);

    void *args = dense_qp_assign_args(fcn_ptrs, dims, ptr);

    return args;
}



void *dense_qp_copy_args(dense_qp_solver_fcn_ptrs *fcn_ptrs, dense_qp_dims *dims, void *raw_memory, void *source)
{
    return fcn_ptrs->copy_args(dims, raw_memory, source);
}



int dense_qp_calculate_size(dense_qp_solver_fcn_ptrs *fcn_ptrs, dense_qp_dims *dims, void *args_)
{
    int bytes = 0;

    bytes += sizeof(dense_qp_solver);

    bytes += sizeof(dense_qp_solver_fcn_ptrs);

    bytes += dense_qp_dims_calculate_size();

    bytes += dense_qp_calculate_args_size(fcn_ptrs, dims);

    bytes += fcn_ptrs->calculate_memory_size(dims, args_);

    bytes += fcn_ptrs->calculate_workspace_size(dims, args_);

    return bytes;
}



dense_qp_solver *dense_qp_assign(dense_qp_solver_fcn_ptrs *fcn_ptrs, dense_qp_dims *dims, void *args_, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    dense_qp_solver *solver = (dense_qp_solver *) c_ptr;
    c_ptr += sizeof(dense_qp_solver);

    solver->fcn_ptrs = (dense_qp_solver_fcn_ptrs *) c_ptr;
    c_ptr += sizeof(dense_qp_solver_fcn_ptrs);

    solver->dims = assign_dense_qp_dims(c_ptr);
    c_ptr += dense_qp_dims_calculate_size();

    solver->args = dense_qp_copy_args(fcn_ptrs, dims, c_ptr, args_);
    c_ptr += dense_qp_calculate_args_size(fcn_ptrs, dims);

    solver->mem = fcn_ptrs->assign_memory(dims, args_, c_ptr);
    c_ptr += fcn_ptrs->calculate_memory_size(dims, args_);

    solver-> work = (void *) c_ptr;
    c_ptr += fcn_ptrs->calculate_workspace_size(dims, args_);

    assert((char*)raw_memory + dense_qp_calculate_size(fcn_ptrs, dims, args_) == c_ptr);

    *solver->fcn_ptrs = *fcn_ptrs;
    solver->fcn_ptrs->submodules = NULL;

    dense_qp_copy_dims(solver->dims, dims);

    return solver;
}



dense_qp_solver *dense_qp_create(dense_qp_solver_fcn_ptrs *fcn_ptrs, dense_qp_dims *dims, void *args_)
{
    int bytes = dense_qp_calculate_size(fcn_ptrs, dims, args_);

    void *ptr = malloc(bytes);

    dense_qp_solver *solver = dense_qp_assign(fcn_ptrs, dims, args_, ptr);

    return solver;
}



int dense_qp_solve(dense_qp_solver *solver, dense_qp_in *qp_in, dense_qp_out *qp_out)
{
    return solver->fcn_ptrs->fun(qp_in, qp_out, solver->args, solver->mem, solver->work);
}



int dense_qp_calculate_submodules_size(dense_qp_solver_config *config, dense_qp_dims *dims)
{
    dense_qp_solver_t solver_name = config->qp_solver;

    int size;

    switch (solver_name)
    {
        case DENSE_QP_HPIPM:
            size = dense_qp_hpipm_calculate_submodules_size(config, dims);
            break;
        case DENSE_QP_QORE:
#ifdef ACADOS_WITH_QORE
            size = dense_qp_qore_calculate_submodules_size(config, dims);
#endif
            break;
        case DENSE_QP_QPOASES:
            size = dense_qp_qpoases_calculate_submodules_size(config, dims);
            break;
        default:
            size = 0;
    }

    return size;
}



void *dense_qp_assign_submodules(dense_qp_solver_config *config, dense_qp_dims *dims, void *raw_memory)
{
    dense_qp_solver_t solver_name = config->qp_solver;

    void *submodules;

    switch (solver_name)
    {
        case DENSE_QP_HPIPM:
            submodules = dense_qp_hpipm_assign_submodules(config, dims, raw_memory);
            break;
        case DENSE_QP_QORE:
#ifdef ACADOS_WITH_QORE
            submodules = dense_qp_qore_assign_submodules(config, dims, raw_memory);
#endif
            break;
        case DENSE_QP_QPOASES:
            submodules = dense_qp_qpoases_assign_submodules(config, dims, raw_memory);
            break;
        default:
            submodules = NULL;
    }

    return submodules;
}



int calculate_dense_qp_solver_fcn_ptrs_size(dense_qp_solver_config *config, dense_qp_dims *dims)
{
    int size = sizeof(dense_qp_solver_fcn_ptrs);

    size += dense_qp_calculate_submodules_size(config, dims);

    return size;
}



void *assign_dense_qp_solver_fcn_ptrs(dense_qp_solver_config *config, dense_qp_dims *dims, void *raw_memory)
{
    char *c_ptr = (char *)raw_memory;

    dense_qp_solver_fcn_ptrs *fcn_ptrs = (dense_qp_solver_fcn_ptrs *)c_ptr;
    c_ptr += sizeof(dense_qp_solver_fcn_ptrs);

    set_dense_qp_solver_fcn_ptrs(config, fcn_ptrs);

    fcn_ptrs->submodules = dense_qp_assign_submodules(config, dims, c_ptr);
    c_ptr += dense_qp_calculate_submodules_size(config, dims);

    assert((char*)raw_memory + calculate_dense_qp_solver_fcn_ptrs_size(config, dims) == c_ptr);

    return (void *)fcn_ptrs;
}



void *create_dense_qp_solver_fcn_ptrs(dense_qp_solver_config *config, dense_qp_dims *dims)
{
    int bytes = calculate_dense_qp_solver_fcn_ptrs_size(config, dims);

    void *ptr = malloc(bytes);

    dense_qp_solver_fcn_ptrs *fcn_ptrs = assign_dense_qp_solver_fcn_ptrs(config, dims, ptr);

    return fcn_ptrs;
}



int set_dense_qp_solver_fcn_ptrs(dense_qp_solver_config *config, dense_qp_solver_fcn_ptrs *fcn_ptrs)
{
    int return_value = ACADOS_SUCCESS;
    dense_qp_solver_t solver_name = config->qp_solver;

    switch (solver_name)
    {
        case DENSE_QP_HPIPM:
            fcn_ptrs->fun = &dense_qp_hpipm;
            fcn_ptrs->calculate_args_size = &dense_qp_hpipm_calculate_args_size;
            fcn_ptrs->assign_args = &dense_qp_hpipm_assign_args;
            fcn_ptrs->copy_args = &dense_qp_hpipm_copy_args;
            fcn_ptrs->initialize_default_args = &dense_qp_hpipm_initialize_default_args;
            fcn_ptrs->calculate_memory_size = &dense_qp_hpipm_calculate_memory_size;
            fcn_ptrs->assign_memory = &dense_qp_hpipm_assign_memory;
            fcn_ptrs->calculate_workspace_size = &dense_qp_hpipm_calculate_workspace_size;
            break;
        case DENSE_QP_QORE:
            #ifdef ACADOS_WITH_QORE
            fcn_ptrs->fun = &dense_qp_qore;
            fcn_ptrs->calculate_args_size = &dense_qp_qore_calculate_args_size;
            fcn_ptrs->assign_args = &dense_qp_qore_assign_args;
            fcn_ptrs->copy_args = &dense_qp_qore_copy_args;
            fcn_ptrs->initialize_default_args = &dense_qp_qore_initialize_default_args;
            fcn_ptrs->calculate_memory_size = &dense_qp_qore_calculate_memory_size;
            fcn_ptrs->assign_memory = &dense_qp_qore_assign_memory;
            fcn_ptrs->calculate_workspace_size = &dense_qp_qore_calculate_workspace_size;
            #else
            return_value = ACADOS_FAILURE;
            #endif
            break;
        case DENSE_QP_QPOASES:
            fcn_ptrs->fun = &dense_qp_qpoases;
            fcn_ptrs->calculate_args_size = &dense_qp_qpoases_calculate_args_size;
            fcn_ptrs->assign_args = &dense_qp_qpoases_assign_args;
            fcn_ptrs->copy_args = &dense_qp_qpoases_copy_args;
            fcn_ptrs->initialize_default_args = &dense_qp_qpoases_initialize_default_args;
            fcn_ptrs->calculate_memory_size = &dense_qp_qpoases_calculate_memory_size;
            fcn_ptrs->assign_memory = &dense_qp_qpoases_assign_memory;
            fcn_ptrs->calculate_workspace_size = &dense_qp_qpoases_calculate_workspace_size;
            break;
        default:
            return_value = ACADOS_FAILURE;
    }

    return return_value;
}
