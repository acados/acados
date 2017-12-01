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
 
// X-Condensing QP solvers
#include <acados/ocp_qp/ocp_qp_condensing_solver.h>
#include <acados/ocp_qp/ocp_qp_sparse_solver.h>

// Submodules
#include <acados/dense_qp/dense_qp_hpipm.h>
#include <acados/dense_qp/dense_qp_qore.h>
#include <acados/dense_qp/dense_qp_qpoases.h>
#include <acados/ocp_qp/ocp_qp_hpipm.h>
// #include <acados/ocp_qp/ocp_qp_hpmpc.h>
// #include <acados/ocp_qp/ocp_qp_ooqp.h>
// #include <acados/ocp_qp/ocp_qp_qpdunes.h>

int set_qp_solver_fun_ptrs(qp_solver_t qp_solver_name, void *qp_solver) {
    int return_value = ACADOS_SUCCESS;

    switch (qp_solver_name) {
        case HPIPM:
            ((ocp_qp_solver *)qp_solver)->calculate_args_size =
                &ocp_qp_hpipm_calculate_args_size;
            ((ocp_qp_solver *)qp_solver)->assign_args =
                &ocp_qp_hpipm_assign_args;
            ((ocp_qp_solver *)qp_solver)->initialize_default_args =
                &ocp_qp_hpipm_initialize_default_args;
            ((ocp_qp_solver *)qp_solver)->calculate_memory_size =
                &ocp_qp_hpipm_calculate_memory_size;
            ((ocp_qp_solver *)qp_solver)->assign_memory =
                &ocp_qp_hpipm_assign_memory;
            ((ocp_qp_solver *)qp_solver)->calculate_workspace_size =
                &ocp_qp_hpipm_calculate_workspace_size;
            ((ocp_qp_solver *)qp_solver)->fun = &ocp_qp_hpipm;
            break;
        case CONDENSING_HPIPM:
            ((dense_qp_solver *)qp_solver)->calculate_args_size =
                &dense_qp_hpipm_calculate_args_size;
            ((dense_qp_solver *)qp_solver)->assign_args =
                &dense_qp_hpipm_assign_args;
            ((dense_qp_solver *)qp_solver)->initialize_default_args =
                &dense_qp_hpipm_initialize_default_args;
            ((dense_qp_solver *)qp_solver)->calculate_memory_size =
                &dense_qp_hpipm_calculate_memory_size;
            ((dense_qp_solver *)qp_solver)->assign_memory =
                &dense_qp_hpipm_assign_memory;
            ((dense_qp_solver *)qp_solver)->calculate_workspace_size =
                &dense_qp_hpipm_calculate_workspace_size;
            ((dense_qp_solver *)qp_solver)->fun = &dense_qp_hpipm;
            break;
        case CONDENSING_QPOASES:
            ((dense_qp_solver *)qp_solver)->calculate_args_size =
                &dense_qp_qpoases_calculate_args_size;
            ((dense_qp_solver *)qp_solver)->assign_args =
                &dense_qp_qpoases_assign_args;
            ((dense_qp_solver *)qp_solver)->initialize_default_args =
                &dense_qp_qpoases_initialize_default_args;
            ((dense_qp_solver *)qp_solver)->calculate_memory_size =
                &dense_qp_qpoases_calculate_memory_size;
            ((dense_qp_solver *)qp_solver)->assign_memory =
                &dense_qp_qpoases_assign_memory;
            ((dense_qp_solver *)qp_solver)->calculate_workspace_size =
                &dense_qp_qpoases_calculate_workspace_size;
            ((dense_qp_solver *)qp_solver)->fun = &dense_qp_qpoases;
            break;
        case CONDENSING_QORE:
            ((dense_qp_solver *)qp_solver)->calculate_args_size =
                &dense_qp_qore_calculate_args_size;
            ((dense_qp_solver *)qp_solver)->assign_args =
                &dense_qp_qore_assign_args;
            ((dense_qp_solver *)qp_solver)->initialize_default_args =
                &dense_qp_qore_initialize_default_args;
            ((dense_qp_solver *)qp_solver)->calculate_memory_size =
                &dense_qp_qore_calculate_memory_size;
            ((dense_qp_solver *)qp_solver)->assign_memory =
                &dense_qp_qore_assign_memory;
            ((dense_qp_solver *)qp_solver)->calculate_workspace_size =
                &dense_qp_qore_calculate_workspace_size;
            ((dense_qp_solver *)qp_solver)->fun = &dense_qp_qore;
        default:
            return_value = ACADOS_FAILURE;
    }
    return return_value;
}



void set_xcond_qp_solver_fun_ptrs(qp_solver_t qp_solver_name,
                                  ocp_qp_xcond_solver *qp_solver) {
    if (qp_solver_name < CONDENSING_HPIPM) {
        qp_solver->calculate_args_size =
            &ocp_qp_sparse_solver_calculate_args_size;
        qp_solver->assign_args = &ocp_qp_sparse_solver_assign_args;
        qp_solver->initialize_default_args =
            &ocp_qp_sparse_solver_initialize_default_args;
        qp_solver->calculate_memory_size =
            &ocp_qp_sparse_solver_calculate_memory_size;
        qp_solver->assign_memory = &ocp_qp_sparse_solver_assign_memory;
        qp_solver->calculate_workspace_size =
            &ocp_qp_sparse_solver_calculate_workspace_size;
        qp_solver->fun = &ocp_qp_sparse_solver;
    } else {
        qp_solver->calculate_args_size =
            &ocp_qp_condensing_solver_calculate_args_size;
        qp_solver->assign_args = &ocp_qp_condensing_solver_assign_args;
        qp_solver->initialize_default_args =
            &ocp_qp_condensing_solver_initialize_default_args;
        qp_solver->calculate_memory_size =
            &ocp_qp_condensing_solver_calculate_memory_size;
        qp_solver->assign_memory = &ocp_qp_condensing_solver_assign_memory;
        qp_solver->calculate_workspace_size =
            &ocp_qp_condensing_solver_calculate_workspace_size;
        qp_solver->fun = &ocp_qp_condensing_solver;
    }
    set_qp_solver_fun_ptrs(qp_solver_name, qp_solver->qp_solver_funs);
}