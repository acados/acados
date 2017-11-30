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

#if defined (EXT_DEPS)

#include "acados/utils/create.h"

#include "acados/dense_qp/dense_qp_common.h"
#include "acados/dense_qp/dense_qp_qpoases.h"
#include "acados/dense_qp/dense_qp_hpipm.h"

#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/ocp_nlp/ocp_nlp_common.h"
#include "acados/ocp_qp/ocp_qp_condensing.h"
#include "acados/ocp_qp/ocp_qp_partial_condensing.h"
#include "acados/ocp_qp/ocp_qp_sparse_solver.h"
#include "acados/ocp_qp/ocp_qp_condensing_solver.h"
#include "acados/ocp_qp/ocp_qp_hpipm.h"
#include "acados/ocp_nlp/ocp_nlp_gn_sqp.h"
#include "acados/utils/mem.h"


ocp_qp_in *create_ocp_qp_in(ocp_qp_dims *dims)
{
    int size = ocp_qp_in_calculate_size(dims);
    void *ptr = acados_malloc(size, 1);
    ocp_qp_in *qp_in = assign_ocp_qp_in(dims, ptr);
    return qp_in;
}



ocp_qp_out *create_ocp_qp_out(ocp_qp_dims *dims)
{
    int size = ocp_qp_out_calculate_size(dims);
    void *ptr = acados_malloc(size, 1);
    ocp_qp_out *qp_out = assign_ocp_qp_out(dims, ptr);
    return qp_out;
}



ocp_qp_hpipm_args *ocp_qp_hpipm_create_arguments(ocp_qp_dims *dims)
{
    int size = ocp_qp_hpipm_calculate_args_size(dims);
    void *ptr = acados_malloc(size, 1);
    void *args = ocp_qp_hpipm_assign_args(dims, ptr);
    ocp_qp_hpipm_initialize_default_args(args);

    return args;
}



ocp_qp_hpipm_memory *ocp_qp_hpipm_create_memory(ocp_qp_dims *dims, void *args_)
{
    ocp_qp_hpipm_args *args = (ocp_qp_hpipm_args *) args_;

    int size = ocp_qp_hpipm_calculate_memory_size(dims, args);
    void *ptr = acados_malloc(size, 1);
    void *mem = ocp_qp_hpipm_assign_memory(dims, args, ptr);

    return mem;
}




ocp_qp_condensing_args *ocp_qp_condensing_create_arguments(ocp_qp_dims *dims)
{
    int size = ocp_qp_condensing_calculate_args_size(dims);
    void *ptr = acados_malloc(size, 1);
    ocp_qp_condensing_args *args = ocp_qp_condensing_assign_args(dims, ptr);
    return args;
}



ocp_qp_condensing_memory *ocp_qp_condensing_create_memory(ocp_qp_dims *dims, ocp_qp_condensing_args *args)
{
    int size = ocp_qp_condensing_calculate_memory_size(dims, args);
    void *ptr = acados_malloc(size, 1);
    ocp_qp_condensing_memory *memory = ocp_qp_condensing_assign_memory(dims, args, ptr);
    return memory;
}



ocp_qp_partial_condensing_args *ocp_qp_partial_condensing_create_arguments(ocp_qp_dims *dims)
{
    int size = ocp_qp_partial_condensing_calculate_args_size(dims);
    void *ptr = acados_malloc(size, 1);
    ocp_qp_partial_condensing_args *args = ocp_qp_partial_condensing_assign_args(dims, ptr);
    ocp_qp_partial_condensing_initialize_default_args(args);
    return args;
}



ocp_qp_partial_condensing_memory *ocp_qp_partial_condensing_create_memory(ocp_qp_dims *dims,
    ocp_qp_partial_condensing_args *args)
{
    int size = ocp_qp_partial_condensing_calculate_memory_size(dims, args);
    void *ptr = acados_malloc(size, 1);
    ocp_qp_partial_condensing_memory *mem = ocp_qp_partial_condensing_assign_memory(dims, args, ptr);
    return mem;
}



ocp_qp_sparse_solver_args *ocp_qp_sparse_solver_create_arguments(ocp_qp_dims *dims, qp_solver_t solver_name)
{
    ocp_qp_solver solver;
    set_qp_solver_fun_ptrs(solver_name, &solver);

    int size = ocp_qp_sparse_solver_calculate_args_size(dims, &solver);
    void *ptr = acados_malloc(size, 1);
    ocp_qp_sparse_solver_args *args = ocp_qp_sparse_solver_assign_args(dims, &solver, ptr);
    ocp_qp_sparse_solver_initialize_default_args(args);
    return args;
}



ocp_qp_sparse_solver_memory *ocp_qp_sparse_solver_create_memory(ocp_qp_dims *dims, void *args_)
{
    ocp_qp_sparse_solver_args *args = (ocp_qp_sparse_solver_args *) args_;

    int size = ocp_qp_sparse_solver_calculate_memory_size(dims, args);
    void *ptr = acados_malloc(size, 1);
    ocp_qp_sparse_solver_memory *mem = ocp_qp_sparse_solver_assign_memory(dims, args, ptr);
    return mem;
}



ocp_qp_condensing_solver_args *ocp_qp_condensing_solver_create_arguments(ocp_qp_dims *dims, qp_solver_t solver_name)
{
    dense_qp_solver solver;
    set_qp_solver_fun_ptrs(solver_name, &solver);

    int size = ocp_qp_condensing_solver_calculate_args_size(dims, &solver);
    void *ptr = acados_malloc(size, 1);
    ocp_qp_condensing_solver_args *args = ocp_qp_condensing_solver_assign_args(dims, &solver, ptr);
    ocp_qp_condensing_solver_initialize_default_args(args);
    return args;
}



ocp_qp_condensing_solver_memory *ocp_qp_condensing_solver_create_memory(ocp_qp_dims *dims, void *args_)
{
    ocp_qp_condensing_solver_args *args = (ocp_qp_condensing_solver_args *) args_;

    int size = ocp_qp_condensing_solver_calculate_memory_size(dims, args);
    void *ptr = acados_malloc(size, 1);
    ocp_qp_condensing_solver_memory *mem = ocp_qp_condensing_solver_assign_memory(dims, args, ptr);
    return mem;
}



dense_qp_in *create_dense_qp_in(dense_qp_dims *dims)
{
    int size = dense_qp_in_calculate_size(dims);
    void *ptr = acados_malloc(size, 1);
    dense_qp_in *qp_in = assign_dense_qp_in(dims, ptr);
    return qp_in;
}



dense_qp_out *create_dense_qp_out(dense_qp_dims *dims)
{
    int size = dense_qp_out_calculate_size(dims);
    void *ptr = acados_malloc(size, 1);
    dense_qp_out *qp_out = assign_dense_qp_out(dims, ptr);
    return qp_out;
}



dense_qp_hpipm_args *dense_qp_hpipm_create_arguments(dense_qp_dims *dims)
{
    int size = dense_qp_hpipm_calculate_args_size(dims);
    void *ptr = acados_malloc(size, 1);
    dense_qp_hpipm_args *args = dense_qp_hpipm_assign_args(dims, ptr);
    dense_qp_hpipm_initialize_default_args(args);

    return args;
}



dense_qp_hpipm_memory *dense_qp_hpipm_create_memory(dense_qp_dims *dims, void *args_)
{
    dense_qp_hpipm_args *args = (dense_qp_hpipm_args *) args_;

    int size = dense_qp_hpipm_calculate_memory_size(dims, args);
    void *ptr = acados_malloc(size, 1);
    dense_qp_hpipm_memory *mem = dense_qp_hpipm_assign_memory(dims, args, ptr);

    return mem;
}



dense_qp_qpoases_args *dense_qp_qpoases_create_arguments(dense_qp_dims *dims)
{
    int size = dense_qp_qpoases_calculate_args_size(dims);
    void *ptr = acados_malloc(size, 1);
    dense_qp_qpoases_args *args = dense_qp_qpoases_assign_args(dims, ptr);
    dense_qp_qpoases_initialize_default_args(args);

    return args;
}



dense_qp_qpoases_memory *dense_qp_qpoases_create_memory(dense_qp_dims *dims, void *args_)
{
    dense_qp_qpoases_args *args = (dense_qp_qpoases_args *) args_;

    int size = dense_qp_qpoases_calculate_memory_size(dims, args);
    void *ptr = acados_malloc(size, 1);
    dense_qp_qpoases_memory *mem = dense_qp_qpoases_assign_memory(dims, args, ptr);
    return mem;
}



ocp_nlp_in *create_ocp_nlp_in(ocp_nlp_dims *dims, int num_stages)
{
    int size = ocp_nlp_in_calculate_size(dims);
    void *ptr = acados_malloc(size, 1);
    ocp_nlp_in *nlp_in = assign_ocp_nlp_in(dims, num_stages, ptr);
    return nlp_in;
}


#ifdef YT
ocp_nlp_gn_sqp_args *ocp_nlp_gn_sqp_create_args(ocp_nlp_dims *dims, qp_solver_t qp_solver_name, sim_solver_t *sim_solver_names)
#else
ocp_nlp_gn_sqp_args *ocp_nlp_gn_sqp_create_args(ocp_nlp_dims *dims, qp_solver_t qp_solver_name)
#endif
{
    ocp_qp_xcond_solver qp_solver;
    dense_qp_solver solver_funs;  // TODO TEMP!
    qp_solver.qp_solver_funs = &solver_funs;

    set_xcond_qp_solver_fun_ptrs(qp_solver_name, &qp_solver);

    #ifdef YT
    int size = ocp_nlp_gn_sqp_calculate_args_size(dims, &qp_solver, sim_solver_names);
    #else
    int size = ocp_nlp_gn_sqp_calculate_args_size(dims, &qp_solver);
    #endif

    void *ptr = acados_malloc(size, 1);

    #ifdef YT
    ocp_nlp_gn_sqp_args *args = ocp_nlp_gn_sqp_assign_args(dims, &qp_solver, sim_solver_names, ptr);
    #else
    ocp_nlp_gn_sqp_args *args = ocp_nlp_gn_sqp_assign_args(dims, &qp_solver, ptr);
    #endif

    // TODO(dimitris): nest in initialize default args of SQP solver!
    args->qp_solver->initialize_default_args(args->qp_solver_args);
    args->maxIter = 30;

    #ifdef YT
    sim_dims sim_dims;
    for (int ii = 0; ii < dims->N; ii++)
    {
        if (sim_solver_names[ii] != PREVIOUS)
        {
            cast_nlp_dims_to_sim_dims(&sim_dims, dims, ii);
            args->sim_solvers[ii]->initialize_default_args(&sim_dims, args->sim_solvers_args[ii]);
        }
        sim_RK_opts *tmp = (sim_RK_opts *)args->sim_solvers_args[ii];
    }
    #endif

    return args;
}



ocp_nlp_gn_sqp_memory *ocp_nlp_gn_sqp_create_memory(ocp_nlp_dims *dims, ocp_nlp_gn_sqp_args *args)
{
    int size = ocp_nlp_gn_sqp_calculate_memory_size(dims, args);
    void *ptr = acados_malloc(size, 1);
    ocp_nlp_gn_sqp_memory *mem = ocp_nlp_gn_sqp_assign_memory(dims, args, ptr);
    return mem;
}


#endif  // EXT_DEPS