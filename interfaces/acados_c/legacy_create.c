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

#include "acados_c/legacy_create.h"

#include <assert.h>
#include <stdlib.h>
#include <acados/utils/mem.h>



ocp_qp_res *create_ocp_qp_res(ocp_qp_dims *dims)
{
    int size = ocp_qp_res_calculate_size(dims);
    void *ptr = acados_malloc(size, 1);
    ocp_qp_res *qp_res = assign_ocp_qp_res(dims, ptr);
    return qp_res;
}



ocp_qp_res_ws *create_ocp_qp_res_ws(ocp_qp_dims *dims)
{
    int size = ocp_qp_res_ws_calculate_size(dims);
    void *ptr = acados_malloc(size, 1);
    ocp_qp_res_ws *res_ws = assign_ocp_qp_res_ws(dims, ptr);
    return res_ws;
}



ocp_qp_full_condensing_args *ocp_qp_full_condensing_create_arguments(ocp_qp_dims *dims, void *submodules_)
{
    int size = ocp_qp_full_condensing_calculate_args_size(dims, submodules_);
    void *ptr = acados_malloc(size, 1);

    void *submodules = submodules_;
    ocp_qp_full_condensing_args *args = ocp_qp_full_condensing_assign_args(dims, &submodules, ptr);

    return args;
}



ocp_qp_full_condensing_memory *ocp_qp_full_condensing_create_memory(ocp_qp_dims *dims, ocp_qp_full_condensing_args *args)
{
    int size = ocp_qp_full_condensing_calculate_memory_size(dims, args);
    void *ptr = acados_malloc(size, 1);
    ocp_qp_full_condensing_memory *memory = ocp_qp_full_condensing_assign_memory(dims, args, ptr);
    return memory;
}



ocp_qp_partial_condensing_args *ocp_qp_partial_condensing_create_arguments(ocp_qp_dims *dims, void *submodules_)
{
    int size = ocp_qp_partial_condensing_calculate_args_size(dims, submodules_);
    void *ptr = acados_malloc(size, 1);

    void *submodules = submodules_;
    ocp_qp_partial_condensing_args *args = ocp_qp_partial_condensing_assign_args(dims, &submodules, ptr);
    
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



// TODO(dimitris): NUM_STAGES NOT NEEDED ANY MORE
ocp_nlp_in *create_ocp_nlp_in(ocp_nlp_dims *dims, int num_stages)
{
    int size = ocp_nlp_in_calculate_size(dims);
    void *ptr = acados_malloc(size, 1);
    ocp_nlp_in *nlp_in = assign_ocp_nlp_in(dims, num_stages, ptr);
    return nlp_in;
}



ocp_nlp_out *create_ocp_nlp_out(ocp_nlp_dims *dims)
{
    int size = ocp_nlp_out_calculate_size(dims);
    void *ptr = acados_malloc(size, 1);
    ocp_nlp_out *nlp_out = assign_ocp_nlp_out(dims, ptr);
    return nlp_out;
}



ocp_nlp_gn_sqp_args *ocp_nlp_gn_sqp_create_args(ocp_nlp_dims *dims, ocp_qp_solver_t qp_solver_name, sim_solver_t *sim_solver_names)
{
    // ocp_qp_xcond_solver_fcn_ptrs qp_solver;
    // module_fcn_ptrs solver_funs;
    // qp_solver.qp_solver = &solver_funs;

    // ocp_qp_solver_config qp_config;
    // qp_config.qp_solver = qp_solver_name;
    // set_ocp_qp_xcond_solver_fcn_ptrs(&qp_config, &qp_solver);

    // int return_value;
    // sim_solver_fcn_ptrs *sim_solvers = acados_malloc(sizeof(sim_solver_fcn_ptrs), dims->N);
    // sim_solver_config sim_config;

    // for (int ii = 0; ii < dims->N; ii++)
    // {
    //     sim_config.sim_solver = sim_solver_names[ii];
    //     return_value = set_sim_solver_fcn_ptrs(&sim_config, &sim_solvers[ii]);
    //     assert(return_value == ACADOS_SUCCESS);
    // }

    // int size = ocp_nlp_gn_sqp_calculate_args_size(dims, &qp_solver, sim_solvers);

    // void *ptr = acados_malloc(size, 1);

    // ocp_nlp_gn_sqp_args *args = ocp_nlp_gn_sqp_assign_args(dims, &qp_solver, sim_solvers, ptr);

    // // TODO(dimitris): nest in initialize default args of SQP solver!
    // args->qp_solver->initialize_default_args(args->qp_solver_args);
    // args->maxIter = 30;

    // sim_dims sim_dims;
    // for (int ii = 0; ii < dims->N; ii++)
    // {
    //     cast_nlp_dims_to_sim_dims(&sim_dims, dims, ii);
    //     args->sim_solvers[ii]->initialize_default_args(&sim_dims, args->sim_solvers_args[ii]);
    // }

    // free(sim_solvers);

    // return args;
    return NULL;
}



ocp_nlp_gn_sqp_memory *ocp_nlp_gn_sqp_create_memory(ocp_nlp_dims *dims, ocp_nlp_gn_sqp_args *args)
{
    // int size = ocp_nlp_gn_sqp_calculate_memory_size(dims, args);
    // void *ptr = acados_malloc(size, 1);
    // ocp_nlp_gn_sqp_memory *mem = ocp_nlp_gn_sqp_assign_memory(dims, args, ptr);
    // return mem;
    return NULL;
}
