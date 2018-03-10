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



// ocp_qp_res *create_ocp_qp_res(ocp_qp_dims *dims)
// {
//     int size = ocp_qp_res_calculate_size(dims);
//     void *ptr = acados_malloc(size, 1);
//     ocp_qp_res *qp_res = ocp_qp_res_assign(dims, ptr);
//     return qp_res;
// }



// ocp_qp_res_ws *create_ocp_qp_res_ws(ocp_qp_dims *dims)
// {
//     int size = ocp_qp_res_workspace_calculate_size(dims);
//     void *ptr = acados_malloc(size, 1);
//     ocp_qp_res_ws *res_ws = ocp_qp_res_workspace_assign(dims, ptr);
//     return res_ws;
// }



ocp_qp_full_condensing_args *ocp_qp_full_condensing_create_arguments(ocp_qp_dims *dims)
{
    int size = ocp_qp_full_condensing_opts_calculate_size(dims);
    void *ptr = acados_malloc(size, 1);
    ocp_qp_full_condensing_args *args = ocp_qp_full_condensing_assign_args(dims, ptr);
    return args;
}



ocp_qp_full_condensing_memory *ocp_qp_full_condensing_create_memory(ocp_qp_dims *dims, ocp_qp_full_condensing_args *args)
{
    int size = ocp_qp_full_condensing_calculate_memory_size(dims, args);
    void *ptr = acados_malloc(size, 1);
    ocp_qp_full_condensing_memory *memory = ocp_qp_full_condensing_assign_memory(dims, args, ptr);
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



ocp_nlp_in *create_ocp_nlp_in(ocp_nlp_solver_config *config, ocp_nlp_dims *dims)
{
    int size = ocp_nlp_in_calculate_size(config, dims);
    void *ptr = acados_malloc(size, 1);
    ocp_nlp_in *nlp_in = ocp_nlp_in_assign(config, dims, ptr);
    return nlp_in;
}



ocp_nlp_out *create_ocp_nlp_out(ocp_nlp_solver_config *config, ocp_nlp_dims *dims)
{
    int size = ocp_nlp_out_calculate_size(config, dims);
    void *ptr = acados_malloc(size, 1);
    ocp_nlp_out *nlp_out = ocp_nlp_out_assign(config, dims, ptr);
    return nlp_out;
}


// ocp_nlp_gn_sqp_opts *ocp_nlp_gn_sqp_create_args(ocp_nlp_dims *dims, ocp_qp_solver_t qp_solver_name, sim_solver_t *sim_solver_names)
// {
//     ocp_qp_xcond_solver_config qp_solver;
//     qp_solver_config solver_funs;
//     qp_solver.qp_solver = &solver_funs;

//     ocp_qp_solver_plan qp_plan;
//     qp_plan.qp_solver = qp_solver_name;
//     set_ocp_qp_xcond_solver_fcn_ptrs(&qp_plan, &qp_solver);

//     int return_value = ACADOS_SUCCESS;
//     sim_solver_config *sim_solvers = acados_malloc(sizeof(sim_solver_config), dims->N);
//     sim_solver_plan sim_plan;

//     for (int ii = 0; ii < dims->N; ii++)
//     {
//         sim_plan.sim_solver = sim_solver_names[ii];
//         return_value = set_sim_solver_fcn_ptrs(&sim_plan, &sim_solvers[ii]);
//         if (return_value != ACADOS_SUCCESS)
//             return NULL;
//     }

//     int size = ocp_nlp_gn_sqp_opts_calculate_size(dims, &qp_solver, sim_solvers);

//     void *ptr = acados_malloc(size, 1);

//     ocp_nlp_gn_sqp_opts *args = ocp_nlp_gn_sqp_assign_args(dims, &qp_solver, sim_solvers, ptr);

//     // TODO(dimitris): nest in initialize default args of SQP solver!
//     args->qp_solver->opts_initialize_default(args->qp_solver_opts);
//     args->maxIter = 30;

//     sim_dims sim_dims;
//     for (int ii = 0; ii < dims->N; ii++)
//     {
//         cast_nlp_dims_to_sim_dims(&sim_dims, dims, ii);
//         args->sim_solvers[ii]->opts_initialize_default(&sim_dims, args->sim_solvers_args[ii]);
//     }

//     free(sim_solvers);
//     return args;
// }



// ocp_nlp_gn_sqp_memory *ocp_nlp_gn_sqp_create_memory(ocp_nlp_dims *dims, ocp_nlp_gn_sqp_opts *args)
// {
//     int size = ocp_nlp_gn_sqp_calculate_memory_size(dims, args);
//     void *ptr = acados_malloc(size, 1);
//     ocp_nlp_gn_sqp_memory *mem = ocp_nlp_gn_sqp_assign_memory(dims, args, ptr);
//     return mem;
// }



// sim_rk_opts *create_sim_irk_opts(sim_dims *dims)
// {
//     int size = sim_irk_opts_calculate_size(dims);

//     void *ptr = acados_malloc(size, 1);

//     sim_rk_opts *opts = sim_irk_opts_assign(dims, ptr);

//     sim_irk_opts_initialize_default(dims, opts);

//     return opts;
// }
