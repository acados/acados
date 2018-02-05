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

#include "acados_c/ocp_qp.h"

//external
#include <stdlib.h>
#include <assert.h>
#include <string.h>
//acados
#include "acados/ocp_qp/ocp_qp_full_condensing.h"
//acados_c
#include "acados_c/ocp_qp/ocp_qp_full_condensing_solver.h"
#include "acados_c/ocp_qp/ocp_qp_sparse_solver.h"
#include "acados_c/dense_qp/dense_qp_hpipm.h"
#ifdef ACADOS_WITH_QORE
#include "acados_c/dense_qp/dense_qp_qore.h"
#endif
#include "acados_c/dense_qp/dense_qp_qpoases.h"
#include "acados_c/ocp_qp/ocp_qp_hpipm.h"
#ifdef ACADOS_WITH_HPMPC
#include "acados_c/ocp_qp/ocp_qp_hpmpc.h"
#endif
#ifdef ACADOS_WITH_QPDUNES
#include "acados_c/ocp_qp/ocp_qp_qpdunes.h"
#endif

// #include <acados/ocp_qp/ocp_qp_hpmpc.h>
// #include <acados/ocp_qp/ocp_qp_ooqp.h>



void ocp_qp_copy_dims(ocp_qp_dims *dest, ocp_qp_dims *src)
{
    dest->N = src->N;
    dest->memsize = src->memsize;

    for (int ii = 0; ii < src->N+1; ii++)
    {
        dest->nx[ii] = src->nx[ii];
        dest->nu[ii] = src->nu[ii];
        dest->nb[ii] = src->nb[ii];
        dest->ng[ii] = src->ng[ii];
        dest->ns[ii] = src->ns[ii];
        dest->nbu[ii] = src->nbu[ii];
        dest->nbx[ii] = src->nbx[ii];
    }
}



ocp_qp_dims *create_ocp_qp_dims(int N)
{
    int bytes = ocp_qp_dims_calculate_size(N);

    void *ptr = calloc(1, bytes);

    ocp_qp_dims *dims = assign_ocp_qp_dims(N, ptr);
    dims->N = N;

    return dims;
}



ocp_qp_in *create_ocp_qp_in(ocp_qp_dims *dims)
{
    int bytes = ocp_qp_in_calculate_size(dims);

    void *ptr = calloc(1, bytes);

    ocp_qp_in *in = assign_ocp_qp_in(dims, ptr);

    return in;
}



ocp_qp_out *create_ocp_qp_out(ocp_qp_dims *dims)
{
    int bytes = ocp_qp_out_calculate_size(dims);

    void *ptr = calloc(1, bytes);

    ocp_qp_out *out = assign_ocp_qp_out(dims, ptr);

    return out;
}

void *set_submodules_fcn_ptrs(
    ocp_qp_solver_config *config,
    ocp_qp_full_condensing_solver_submodules *full_condensing_solver_submodules,
    ocp_qp_sparse_solver_submodules *sparse_solver_submodules)
{
    void *submodules;

    ocp_qp_solver_t solver_name = config->qp_solver;

    int status = ACADOS_SUCCESS;
    dense_qp_solver_config full_condensing_solver_submodules_config;
    ocp_qp_solver_config sparse_solver_submodules_config;
    switch (solver_name) {
        case PARTIAL_CONDENSING_HPIPM:
            sparse_solver_submodules_config.qp_solver = SPARSE_QP_HPIPM;
            status = set_ocp_qp_solver_fcn_ptrs(&sparse_solver_submodules_config, &(sparse_solver_submodules->qpsol));
            submodules = (void *) sparse_solver_submodules;
            break;
        case PARTIAL_CONDENSING_HPMPC:
            sparse_solver_submodules_config.qp_solver = SPARSE_QP_HPMPC;
            status = set_ocp_qp_solver_fcn_ptrs(&sparse_solver_submodules_config, &(sparse_solver_submodules->qpsol));
            submodules = (void *) sparse_solver_submodules;
            break;
        case PARTIAL_CONDENSING_OOQP:
            sparse_solver_submodules_config.qp_solver = SPARSE_QP_OOQP;
            status = set_ocp_qp_solver_fcn_ptrs(&sparse_solver_submodules_config, &(sparse_solver_submodules->qpsol));
            submodules = (void *) sparse_solver_submodules;
            break;
        case PARTIAL_CONDENSING_QPDUNES:
            sparse_solver_submodules_config.qp_solver = SPARSE_QP_QPDUNES;
            status = set_ocp_qp_solver_fcn_ptrs(&sparse_solver_submodules_config, &(sparse_solver_submodules->qpsol));
            submodules = (void *) sparse_solver_submodules;
            break;
        case FULL_CONDENSING_HPIPM:
            full_condensing_solver_submodules_config.qp_solver = DENSE_QP_HPIPM;
            status = set_dense_qp_solver_fcn_ptrs(&full_condensing_solver_submodules_config, &(full_condensing_solver_submodules->qpsol));
            submodules = (void *) full_condensing_solver_submodules;
            break;
        case FULL_CONDENSING_QPOASES:
            full_condensing_solver_submodules_config.qp_solver = DENSE_QP_QPOASES;
            status = set_dense_qp_solver_fcn_ptrs(&full_condensing_solver_submodules_config, &(full_condensing_solver_submodules->qpsol));
            submodules = (void *) full_condensing_solver_submodules;
            break;
        case FULL_CONDENSING_QORE:
            full_condensing_solver_submodules_config.qp_solver = DENSE_QP_QORE;
            status = set_dense_qp_solver_fcn_ptrs(&full_condensing_solver_submodules_config, &(full_condensing_solver_submodules->qpsol));
            submodules = (void *) full_condensing_solver_submodules;
            break;
        default:
            submodules = NULL;
    }

    if (status != ACADOS_SUCCESS)
    {
        printf("\n\nSpecified solver interface not compiled with acados!\n\n");
        exit(1);
    }

    return submodules;
}



int ocp_qp_calculate_args_size(ocp_qp_solver_config *config, ocp_qp_dims *dims)
{
    ocp_qp_solver_fcn_ptrs fcn_ptrs;

    int status = set_ocp_qp_solver_fcn_ptrs(config, &fcn_ptrs);

    ocp_qp_full_condensing_solver_submodules full_condensing_solver_submodules;
    ocp_qp_sparse_solver_submodules sparse_solver_submodules;
    void *submodules = set_submodules_fcn_ptrs(config, &full_condensing_solver_submodules, &sparse_solver_submodules);

    int size = -1;

    if (status == ACADOS_SUCCESS)
    {
        size = fcn_ptrs.calculate_args_size(dims, submodules);
    }
    else
    {
        printf("\n\nSpecified solver interface not compiled with acados!\n\n");
        exit(1);
    }
    return size;
}



void *ocp_qp_assign_args(ocp_qp_solver_config *config, ocp_qp_dims *dims, void *raw_memory)
{
    ocp_qp_solver_fcn_ptrs fcn_ptrs;

    set_ocp_qp_solver_fcn_ptrs(config, &fcn_ptrs);

    ocp_qp_full_condensing_solver_submodules full_condensing_solver_submodules;
    ocp_qp_sparse_solver_submodules sparse_solver_submodules;
    void *submodules = set_submodules_fcn_ptrs(config, &full_condensing_solver_submodules, &sparse_solver_submodules);

    void *args = fcn_ptrs.assign_args(dims, submodules, raw_memory);

    fcn_ptrs.initialize_default_args(args);

    return args;
}



void *ocp_qp_create_args(ocp_qp_solver_config *config, ocp_qp_dims *dims)
{
    int bytes = ocp_qp_calculate_args_size(config, dims);

    void *ptr = calloc(1, bytes);

    void *args = ocp_qp_assign_args(config, dims, ptr);

    return args;
}



void *ocp_qp_copy_args(ocp_qp_solver_config *config, ocp_qp_dims *dims, void *raw_memory, void *source)
{
    ocp_qp_solver_t solver_name = config->qp_solver;

    void *args;

    switch (solver_name) {
        case SPARSE_QP_HPIPM:
            ocp_qp_hpipm_copy_args(config, dims, raw_memory, source);
            break;
        case SPARSE_QP_HPMPC:
#ifdef ACADOS_WITH_HPMPC
            ocp_qp_hpmpc_copy_args(config, dims, raw_memory, source);
#endif
            break;
        case SPARSE_QP_OOQP:
            // ocp_qp_ooqp_copy_args(config, dims, raw_memory, source);
            break;
        case SPARSE_QP_QPDUNES:
#ifdef ACADOS_WITH_QPDUNES
            ocp_qp_qpdunes_copy_args(config, dims, raw_memory, source);
#endif
            break;
        case PARTIAL_CONDENSING_HPIPM:
            ocp_qp_sparse_solver_copy_args(config, dims, raw_memory, source);
            break;
        case PARTIAL_CONDENSING_HPMPC:
            ocp_qp_sparse_solver_copy_args(config, dims, raw_memory, source);
            break;
        case PARTIAL_CONDENSING_OOQP:
            ocp_qp_sparse_solver_copy_args(config, dims, raw_memory, source);
            break;
        case PARTIAL_CONDENSING_QPDUNES:
            ocp_qp_sparse_solver_copy_args(config, dims, raw_memory, source);
            break;
        case FULL_CONDENSING_HPIPM:
            ocp_qp_full_condensing_solver_copy_args(config, dims, raw_memory, source);
            break;
        case FULL_CONDENSING_QPOASES:
            ocp_qp_full_condensing_solver_copy_args(config, dims, raw_memory, source);
            break;
        case FULL_CONDENSING_QORE:
            ocp_qp_full_condensing_solver_copy_args(config, dims, raw_memory, source);
            break;
    }

    return args;
}



int ocp_qp_calculate_size(ocp_qp_solver_config *config, ocp_qp_dims *dims, void *args_)
{
    ocp_qp_solver_fcn_ptrs fcn_ptrs;

    set_ocp_qp_solver_fcn_ptrs(config, &fcn_ptrs);

    ocp_qp_full_condensing_solver_submodules full_condensing_solver_submodules;
    ocp_qp_sparse_solver_submodules sparse_solver_submodules;
    void *submodules = set_submodules_fcn_ptrs(config, &full_condensing_solver_submodules, &sparse_solver_submodules);

    int bytes;

    bytes = sizeof(ocp_qp_solver);

    bytes += sizeof(ocp_qp_solver_fcn_ptrs);

    bytes += ocp_qp_dims_calculate_size(dims->N);

    bytes += fcn_ptrs.calculate_args_size(dims, submodules);

    bytes += fcn_ptrs.calculate_memory_size(dims, args_);

    bytes += fcn_ptrs.calculate_workspace_size(dims, args_);

    return bytes;
}



ocp_qp_solver *ocp_qp_assign(ocp_qp_solver_config *config, ocp_qp_dims *dims, void *args_, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    ocp_qp_solver *solver = (ocp_qp_solver *) c_ptr;
    c_ptr += sizeof(ocp_qp_solver);

    solver->fcn_ptrs = (ocp_qp_solver_fcn_ptrs *) c_ptr;
    c_ptr += sizeof(ocp_qp_solver_fcn_ptrs);

    solver->dims = assign_ocp_qp_dims(dims->N, c_ptr);
    c_ptr += ocp_qp_dims_calculate_size(dims->N);
    ocp_qp_copy_dims(solver->dims, dims);

    solver->args = ocp_qp_copy_args(config, dims, c_ptr, args_);
    c_ptr += ocp_qp_calculate_args_size(config, dims);

    solver->mem = solver->fcn_ptrs->assign_memory(dims, args_, c_ptr);
    c_ptr += solver->fcn_ptrs->calculate_memory_size(dims, args_);

    solver->work = (void *) c_ptr;
    c_ptr += solver->fcn_ptrs->calculate_workspace_size(dims, args_);

    assert((char*)raw_memory + ocp_qp_calculate_size(config, dims, args_) == c_ptr);

    return solver;
}



ocp_qp_solver *ocp_qp_create(ocp_qp_solver_config *config, ocp_qp_dims *dims, void *args_)
{
    int bytes = ocp_qp_calculate_size(config, dims, args_);

    void *ptr = calloc(1, bytes);

    ocp_qp_solver *solver = ocp_qp_assign(config, dims, args_, ptr);

    return solver;
}



int ocp_qp_solve(ocp_qp_solver *solver, ocp_qp_in *qp_in, ocp_qp_out *qp_out)
{
    return solver->fcn_ptrs->fun(qp_in, qp_out, solver->args, solver->mem, solver->work);
}



void set_ocp_qp_full_condensing_solver_fcn_ptrs(ocp_qp_solver_fcn_ptrs *fcn_ptrs)
{
    fcn_ptrs->calculate_args_size = &ocp_qp_full_condensing_solver_calculate_args_size;
    fcn_ptrs->assign_args = &ocp_qp_full_condensing_solver_assign_args;
    fcn_ptrs->initialize_default_args = &ocp_qp_full_condensing_solver_initialize_default_args;
    fcn_ptrs->calculate_memory_size = &ocp_qp_full_condensing_solver_calculate_memory_size;
    fcn_ptrs->assign_memory = &ocp_qp_full_condensing_solver_assign_memory;
    fcn_ptrs->calculate_workspace_size = &ocp_qp_full_condensing_solver_calculate_workspace_size;
    fcn_ptrs->fun = &ocp_qp_full_condensing_solver;
}



void set_ocp_qp_partial_condensing_solver_fcn_ptrs(ocp_qp_solver_fcn_ptrs *fcn_ptrs)
{
    fcn_ptrs->calculate_args_size = &ocp_qp_sparse_solver_calculate_args_size;
    fcn_ptrs->assign_args = &ocp_qp_sparse_solver_assign_args;
    fcn_ptrs->initialize_default_args = &ocp_qp_sparse_solver_initialize_default_args;
    fcn_ptrs->calculate_memory_size = &ocp_qp_sparse_solver_calculate_memory_size;
    fcn_ptrs->assign_memory = &ocp_qp_sparse_solver_assign_memory;
    fcn_ptrs->calculate_workspace_size = &ocp_qp_sparse_solver_calculate_workspace_size;
    fcn_ptrs->fun = &ocp_qp_sparse_solver;
}



int set_ocp_qp_solver_fcn_ptrs(ocp_qp_solver_config *config, ocp_qp_solver_fcn_ptrs *fcn_ptrs)
{
    int return_value = ACADOS_SUCCESS;
    ocp_qp_solver_t solver_name = config->qp_solver;

    switch (solver_name) {
        case SPARSE_QP_HPIPM:
            fcn_ptrs->calculate_args_size = &ocp_qp_hpipm_calculate_args_size;
            fcn_ptrs->assign_args = &ocp_qp_hpipm_assign_args;
            fcn_ptrs->initialize_default_args = &ocp_qp_hpipm_initialize_default_args;
            fcn_ptrs->calculate_memory_size = &ocp_qp_hpipm_calculate_memory_size;
            fcn_ptrs->assign_memory = &ocp_qp_hpipm_assign_memory;
            fcn_ptrs->calculate_workspace_size = &ocp_qp_hpipm_calculate_workspace_size;
            fcn_ptrs->fun = &ocp_qp_hpipm;
            break;
        case SPARSE_QP_HPMPC:
#ifdef ACADOS_WITH_HPMPC
            fcn_ptrs->calculate_args_size = &ocp_qp_hpmpc_calculate_args_size;
            fcn_ptrs->assign_args = &ocp_qp_hpmpc_assign_args;
            fcn_ptrs->initialize_default_args = &ocp_qp_hpmpc_initialize_default_args;
            fcn_ptrs->calculate_memory_size = &ocp_qp_hpmpc_calculate_memory_size;
            fcn_ptrs->assign_memory = &ocp_qp_hpmpc_assign_memory;
            fcn_ptrs->calculate_workspace_size = &ocp_qp_hpmpc_calculate_workspace_size;
            fcn_ptrs->fun = &ocp_qp_hpmpc;
#else
            return_value = ACADOS_FAILURE;
#endif
            break;
        case SPARSE_QP_OOQP:
            // fcn_ptrs->calculate_args_size = &ocp_qp_ooqp_calculate_args_size;
            // fcn_ptrs->assign_args = &ocp_qp_ooqp_assign_args;
            // fcn_ptrs->initialize_default_args = &ocp_qp_ooqp_initialize_default_args;
            // fcn_ptrs->calculate_memory_size = &ocp_qp_ooqp_calculate_memory_size;
            // fcn_ptrs->assign_memory = &ocp_qp_ooqp_assign_memory;
            // fcn_ptrs->calculate_workspace_size = &ocp_qp_ooqp_calculate_workspace_size;
            // fcn_ptrs->fun = &ocp_qp_ooqp;
            break;
        case SPARSE_QP_QPDUNES:
#ifdef ACADOS_WITH_QPDUNES
            fcn_ptrs->calculate_args_size = &ocp_qp_qpdunes_calculate_args_size;
            fcn_ptrs->assign_args = &ocp_qp_qpdunes_assign_args;
            fcn_ptrs->initialize_default_args = &ocp_qp_qpdunes_initialize_default_args;
            fcn_ptrs->calculate_memory_size = &ocp_qp_qpdunes_calculate_memory_size;
            fcn_ptrs->assign_memory = &ocp_qp_qpdunes_assign_memory;
            fcn_ptrs->calculate_workspace_size = &ocp_qp_qpdunes_calculate_workspace_size;
            fcn_ptrs->fun = &ocp_qp_qpdunes;
#else
            return_value = ACADOS_FAILURE;
#endif
            break;
        case PARTIAL_CONDENSING_HPIPM:
            set_ocp_qp_partial_condensing_solver_fcn_ptrs(fcn_ptrs);
            break;
        case PARTIAL_CONDENSING_HPMPC:
            set_ocp_qp_partial_condensing_solver_fcn_ptrs(fcn_ptrs);
            break;
        case PARTIAL_CONDENSING_OOQP:
            set_ocp_qp_partial_condensing_solver_fcn_ptrs(fcn_ptrs);
            break;
        case PARTIAL_CONDENSING_QPDUNES:
            set_ocp_qp_partial_condensing_solver_fcn_ptrs(fcn_ptrs);
            break;
        case FULL_CONDENSING_HPIPM:
            set_ocp_qp_full_condensing_solver_fcn_ptrs(fcn_ptrs);
            break;
        case FULL_CONDENSING_QPOASES:
            set_ocp_qp_full_condensing_solver_fcn_ptrs(fcn_ptrs);
            break;
        case FULL_CONDENSING_QORE:
            set_ocp_qp_full_condensing_solver_fcn_ptrs(fcn_ptrs);
            break;
        default:
            return_value = ACADOS_FAILURE;
    }

    return return_value;
}
