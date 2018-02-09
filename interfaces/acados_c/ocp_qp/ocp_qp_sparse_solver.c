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

#include "acados_c/ocp_qp/ocp_qp_sparse_solver.h"

#include <stdlib.h>

#include "acados_c/ocp_qp/ocp_qp_hpipm.h"
#ifdef ACADOS_WITH_HPMPC
#include "acados_c/ocp_qp/ocp_qp_hpmpc.h"
#endif
#include "acados_c/ocp_qp/ocp_qp_partial_condensing.h"
#ifdef ACADOS_WITH_QPDUNES
#include "acados_c/ocp_qp/ocp_qp_qpdunes.h"
#endif

// void *ocp_qp_sparse_solver_copy_args(ocp_qp_solver_config *config, ocp_qp_dims *dims, void *raw_memory, void *source_)
// {
//     ocp_qp_sparse_solver_args *source = (ocp_qp_sparse_solver_args *) source_;
//     ocp_qp_sparse_solver_args *dest;

//     dest = ocp_qp_assign_args(config, dims, raw_memory);

//     ocp_qp_partial_condensing_copy_args(config, dims, dest->pcond_args, source->pcond_args);

//     ocp_qp_solver_t solver_name = config->qp_solver;

//     ocp_qp_solver_config solver_config;
//     switch(solver_name) {
//         case PARTIAL_CONDENSING_HPIPM:
//             solver_config.qp_solver = SPARSE_QP_HPIPM;
//             ocp_qp_hpipm_copy_args(&solver_config, dims, dest->solver_args, source->solver_args);
//             break;
//         case PARTIAL_CONDENSING_HPMPC:
//             solver_config.qp_solver = SPARSE_QP_HPMPC;
//             ocp_qp_hpmpc_copy_args(&solver_config, dims, dest->solver_args, source->solver_args);
//             break;
//         case PARTIAL_CONDENSING_OOQP:
//             // solver_config.qp_solver = SPARSE_QP_OOQP;
//             // ocp_qp_ooqp_copy_args(&solver_config, dims, dest->solver_args, source->solver_args);
//             break;
//         case PARTIAL_CONDENSING_QPDUNES:
//             solver_config.qp_solver = SPARSE_QP_QPDUNES;
//             ocp_qp_qpdunes_copy_args(&solver_config, dims, dest->solver_args, source->solver_args);
//             break;
//         default:
//             printf(
//                 "\n\nSpecified solver is does not employ partial condensing!\n\n");
//             exit(1);
//     }

//     return (void *) dest;
// }



int ocp_qp_sparse_solver_calculate_submodules_size(ocp_qp_solver_config *config, ocp_qp_dims *dims)
{
    int size = sizeof(ocp_qp_sparse_solver_submodules);

    ocp_qp_solver_t solver_name = config->qp_solver;

    ocp_qp_solver_config subconfig;
    switch (solver_name) {
        case PARTIAL_CONDENSING_HPIPM:
            subconfig.qp_solver = SPARSE_QP_HPIPM;
            break;
        case PARTIAL_CONDENSING_HPMPC:
            subconfig.qp_solver = SPARSE_QP_HPMPC;
            break;
        case PARTIAL_CONDENSING_OOQP:
            subconfig.qp_solver = SPARSE_QP_OOQP;
            break;
        case PARTIAL_CONDENSING_QPDUNES:
            subconfig.qp_solver = SPARSE_QP_QPDUNES;
            break;
    }

    size += calculate_ocp_qp_solver_fcn_ptrs_size(&subconfig, dims);

    return size;
}



void *ocp_qp_sparse_solver_assign_submodules(ocp_qp_solver_config *config, ocp_qp_dims *dims, void *raw_memory)
{
    ocp_qp_sparse_solver_submodules *submodules;

    ocp_qp_solver_t solver_name = config->qp_solver;

    ocp_qp_solver_config subconfig;
    switch (solver_name) {
        case PARTIAL_CONDENSING_HPIPM:
            subconfig.qp_solver = SPARSE_QP_HPIPM;
            break;
        case PARTIAL_CONDENSING_HPMPC:
            subconfig.qp_solver = SPARSE_QP_HPMPC;
            break;
        case PARTIAL_CONDENSING_OOQP:
            subconfig.qp_solver = SPARSE_QP_OOQP;
            break;
        case PARTIAL_CONDENSING_QPDUNES:
            subconfig.qp_solver = SPARSE_QP_QPDUNES;
            break;
    }

    char *c_ptr = (char *) raw_memory;

    submodules = (ocp_qp_sparse_solver_submodules *) c_ptr;
    c_ptr += sizeof(ocp_qp_sparse_solver_submodules);

    submodules->solver = assign_ocp_qp_solver_fcn_ptrs(&subconfig, dims, c_ptr);
    c_ptr += calculate_ocp_qp_solver_fcn_ptrs_size(&subconfig, dims);

    assert((size_t)c_ptr % 8 == 0 && "memory not 8-byte aligned!");

    assert((char*)raw_memory + ocp_qp_sparse_solver_calculate_submodules_size(config, dims) == c_ptr);

    return (void *)submodules;
}
