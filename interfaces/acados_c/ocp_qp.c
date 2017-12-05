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
#include <string.h>
//acados
#include <acados/ocp_qp/ocp_qp_full_condensing_solver.h>
#include <acados/ocp_qp/ocp_qp_sparse_solver.h>
#include <acados/dense_qp/dense_qp_hpipm.h>
#include <acados/dense_qp/dense_qp_qore.h>
#include <acados/dense_qp/dense_qp_qpoases.h>
#include <acados/ocp_qp/ocp_qp_hpipm.h>
// #include <acados/ocp_qp/ocp_qp_hpmpc.h>
// #include <acados/ocp_qp/ocp_qp_ooqp.h>
// #include <acados/ocp_qp/ocp_qp_qpdunes.h>

// int set_qp_solver_fun_ptrs(ocp_qp_solver_t qp_solver_name, void *qp_solver) {
//     int return_value = ACADOS_SUCCESS;

//     switch (qp_solver_name) {
//         case HPIPM:
//             ((ocp_qp_solver *)qp_solver)->calculate_args_size =
//                 &ocp_qp_hpipm_calculate_args_size;
//             ((ocp_qp_solver *)qp_solver)->assign_args =
//                 &ocp_qp_hpipm_assign_args;
//             ((ocp_qp_solver *)qp_solver)->initialize_default_args =
//                 &ocp_qp_hpipm_initialize_default_args;
//             ((ocp_qp_solver *)qp_solver)->calculate_memory_size =
//                 &ocp_qp_hpipm_calculate_memory_size;
//             ((ocp_qp_solver *)qp_solver)->assign_memory =
//                 &ocp_qp_hpipm_assign_memory;
//             ((ocp_qp_solver *)qp_solver)->calculate_workspace_size =
//                 &ocp_qp_hpipm_calculate_workspace_size;
//             ((ocp_qp_solver *)qp_solver)->fun = &ocp_qp_hpipm;
//             break;
//         case CONDENSING_HPIPM:
//             ((dense_qp_solver *)qp_solver)->calculate_args_size =
//                 &dense_qp_hpipm_calculate_args_size;
//             ((dense_qp_solver *)qp_solver)->assign_args =
//                 &dense_qp_hpipm_assign_args;
//             ((dense_qp_solver *)qp_solver)->initialize_default_args =
//                 &dense_qp_hpipm_initialize_default_args;
//             ((dense_qp_solver *)qp_solver)->calculate_memory_size =
//                 &dense_qp_hpipm_calculate_memory_size;
//             ((dense_qp_solver *)qp_solver)->assign_memory =
//                 &dense_qp_hpipm_assign_memory;
//             ((dense_qp_solver *)qp_solver)->calculate_workspace_size =
//                 &dense_qp_hpipm_calculate_workspace_size;
//             ((dense_qp_solver *)qp_solver)->fun = &dense_qp_hpipm;
//             break;
//         case CONDENSING_QPOASES:
//             ((dense_qp_solver *)qp_solver)->calculate_args_size =
//                 &dense_qp_qpoases_calculate_args_size;
//             ((dense_qp_solver *)qp_solver)->assign_args =
//                 &dense_qp_qpoases_assign_args;
//             ((dense_qp_solver *)qp_solver)->initialize_default_args =
//                 &dense_qp_qpoases_initialize_default_args;
//             ((dense_qp_solver *)qp_solver)->calculate_memory_size =
//                 &dense_qp_qpoases_calculate_memory_size;
//             ((dense_qp_solver *)qp_solver)->assign_memory =
//                 &dense_qp_qpoases_assign_memory;
//             ((dense_qp_solver *)qp_solver)->calculate_workspace_size =
//                 &dense_qp_qpoases_calculate_workspace_size;
//             ((dense_qp_solver *)qp_solver)->fun = &dense_qp_qpoases;
//             break;
//         case CONDENSING_QORE:
//             ((dense_qp_solver *)qp_solver)->calculate_args_size =
//                 &dense_qp_qore_calculate_args_size;
//             ((dense_qp_solver *)qp_solver)->assign_args =
//                 &dense_qp_qore_assign_args;
//             ((dense_qp_solver *)qp_solver)->initialize_default_args =
//                 &dense_qp_qore_initialize_default_args;
//             ((dense_qp_solver *)qp_solver)->calculate_memory_size =
//                 &dense_qp_qore_calculate_memory_size;
//             ((dense_qp_solver *)qp_solver)->assign_memory =
//                 &dense_qp_qore_assign_memory;
//             ((dense_qp_solver *)qp_solver)->calculate_workspace_size =
//                 &dense_qp_qore_calculate_workspace_size;
//             ((dense_qp_solver *)qp_solver)->fun = &dense_qp_qore;
//         default:
//             return_value = ACADOS_FAILURE;
//     }
//     return return_value;
// }

// void set_xcond_qp_solver_fun_ptrs(ocp_qp_solver_t qp_solver_name,
//                                   ocp_qp_xcond_solver *qp_solver) {
//     if (qp_solver_name < CONDENSING_HPIPM) {
//         qp_solver->calculate_args_size =
//             &ocp_qp_sparse_solver_calculate_args_size;
//         qp_solver->assign_args = &ocp_qp_sparse_solver_assign_args;
//         qp_solver->initialize_default_args =
//             &ocp_qp_sparse_solver_initialize_default_args;
//         qp_solver->calculate_memory_size =
//             &ocp_qp_sparse_solver_calculate_memory_size;
//         qp_solver->assign_memory = &ocp_qp_sparse_solver_assign_memory;
//         qp_solver->calculate_workspace_size =
//             &ocp_qp_sparse_solver_calculate_workspace_size;
//         qp_solver->fun = &ocp_qp_sparse_solver;
//     } else {
//         qp_solver->calculate_args_size =
//             &ocp_qp_condensing_solver_calculate_args_size;
//         qp_solver->assign_args = &ocp_qp_condensing_solver_assign_args;
//         qp_solver->initialize_default_args =
//             &ocp_qp_condensing_solver_initialize_default_args;
//         qp_solver->calculate_memory_size =
//             &ocp_qp_condensing_solver_calculate_memory_size;
//         qp_solver->assign_memory = &ocp_qp_condensing_solver_assign_memory;
//         qp_solver->calculate_workspace_size =
//             &ocp_qp_condensing_solver_calculate_workspace_size;
//         qp_solver->fun = &ocp_qp_condensing_solver;
//     }
//     set_qp_solver_fun_ptrs(qp_solver_name, qp_solver->qp_solver_funs);
// }




int calculate_dims_size(ocp_qp_dims *dims)
{
    return 0;
}



ocp_qp_in *create_ocp_qp_in(ocp_qp_dims *dims)
{
    int bytes = ocp_qp_in_calculate_size(dims);

    void *ptr = malloc(bytes);

    ocp_qp_in *in = assign_ocp_qp_in(dims, ptr);

    return in;
}



ocp_qp_out *create_ocp_qp_out(ocp_qp_dims *dims)
{
    int bytes = ocp_qp_out_calculate_size(dims);

    void *ptr = malloc(bytes);

    ocp_qp_out *out = assign_ocp_qp_out(dims, ptr);

    return out;
}



int ocp_qp_calculate_args_size(ocp_qp_solver_plan *plan, ocp_qp_dims *dims)
{
    ocp_qp_xcond_solver_fcn_ptrs fcn_ptrs;
    module_fcn_ptrs submodule_fcn_ptrs;
    set_ocp_qp_xcond_solver_fcn_ptrs(plan, &fcn_ptrs, &submodule_fcn_ptrs);

    return fcn_ptrs.calculate_args_size(dims, fcn_ptrs.qp_solver);
}



void *ocp_qp_assign_args(ocp_qp_solver_plan *plan, ocp_qp_dims *dims, void *raw_memory)
{
    ocp_qp_xcond_solver_fcn_ptrs fcn_ptrs;
    module_fcn_ptrs submodule_fcn_ptrs;
    set_ocp_qp_xcond_solver_fcn_ptrs(plan, &fcn_ptrs, &submodule_fcn_ptrs);

    return fcn_ptrs.assign_args(dims, fcn_ptrs.qp_solver, raw_memory);
}



void *ocp_qp_create_args(ocp_qp_solver_plan *plan, ocp_qp_dims *dims)
{
    int bytes = ocp_qp_calculate_args_size(plan, dims);

    void *ptr = malloc(bytes);

    void *args = ocp_qp_assign_args(plan, dims, ptr);

    return args;
}



int ocp_qp_calculate_size(ocp_qp_solver_plan *plan, ocp_qp_dims *dims, void *args_)
{
    ocp_qp_xcond_solver_fcn_ptrs fcn_ptrs;
    module_fcn_ptrs submodule_fcn_ptrs;
    set_ocp_qp_xcond_solver_fcn_ptrs(plan, &fcn_ptrs, &submodule_fcn_ptrs);


    int bytes;
    
    bytes = sizeof(ocp_qp_solver);
    
    bytes += sizeof(ocp_qp_solver_fcn_ptrs);
    
    bytes += calculate_dims_size(dims);

    bytes += ocp_qp_calculate_args_size(plan, dims);

    bytes += fcn_ptrs.calculate_memory_size(dims, args_);

    bytes += fcn_ptrs.calculate_workspace_size(dims, args_);

    return bytes;
}



ocp_qp_solver *ocp_qp_assign(ocp_qp_solver_plan *plan, ocp_qp_dims *dims, void *args_, void *raw_memory)
{
    ocp_qp_xcond_solver_fcn_ptrs fcn_ptrs;
    module_fcn_ptrs submodule_fcn_ptrs;
    set_ocp_qp_xcond_solver_fcn_ptrs(plan, &fcn_ptrs, &submodule_fcn_ptrs);

    char *c_ptr = (char *) raw_memory;

    ocp_qp_solver *solver;
    solver = (ocp_qp_solver *) c_ptr;
    c_ptr += sizeof(ocp_qp_solver);

    solver->fcn_ptrs = (ocp_qp_solver_fcn_ptrs *) c_ptr;
    set_ocp_qp_solver_fcn_ptrs(plan, solver->fcn_ptrs);
    c_ptr += sizeof(ocp_qp_solver_fcn_ptrs);

    solver->dims = (ocp_qp_dims *) c_ptr;
    int dims_size = calculate_dims_size(dims);
    c_ptr += dims_size;
    memcpy(solver->dims, dims, dims_size);

    solver->args = ocp_qp_assign_args(plan, dims, c_ptr);
    int args_size = ocp_qp_calculate_args_size(plan, dims);
    c_ptr += args_size;
    memcpy(solver->args, args_, args_size);

    solver->mem = fcn_ptrs.assign_memory(dims, args_, c_ptr);
    c_ptr += fcn_ptrs.calculate_memory_size(dims, args_);

    solver->work = (void *) c_ptr;
    c_ptr += fcn_ptrs.calculate_workspace_size(dims, args_);

    assert((char*)raw_memory + ocp_qp_calculate_size(plan, dims, args_) == c_ptr);

    return solver;
}



ocp_qp_solver *ocp_qp_create(ocp_qp_solver_plan *plan, ocp_qp_dims *dims, void *args_)
{
    int bytes = ocp_qp_calculate_size(plan, dims, args_);

    void *ptr = malloc(bytes);

    ocp_qp_solver *solver = ocp_qp_assign(plan, dims, args_, ptr);

    return solver;
}



int ocp_qp_solve(ocp_qp_solver *solver, ocp_qp_in *qp_in, ocp_qp_out *qp_out)
{
    return solver->fcn_ptrs->fun(qp_in, qp_out, solver->args, solver->mem, solver->work);
}



void ocp_qp_initialize_default_args(ocp_qp_solver *solver)
{
    solver->fcn_ptrs->initialize_default_args(solver->args);
}



int set_ocp_qp_solver_fcn_ptrs(ocp_qp_solver_plan *plan, ocp_qp_solver_fcn_ptrs *fcn_ptrs)
{
    int return_value = ACADOS_SUCCESS;

    return return_value;
}



int set_ocp_qp_xcond_solver_fcn_ptrs(ocp_qp_solver_plan *plan, ocp_qp_xcond_solver_fcn_ptrs *fcn_ptrs, module_fcn_ptrs *submodule_fcn_ptrs)
{
    int return_value = ACADOS_SUCCESS;

    return return_value;
}
