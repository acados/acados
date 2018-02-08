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

#include "acados_c/ocp_qp/ocp_qp_full_condensing_solver.h"

#include <stdlib.h>

#include "acados_c/dense_qp/dense_qp_hpipm.h"
#include "acados_c/dense_qp/dense_qp_qpoases.h"
#include "acados_c/dense_qp/dense_qp_qore.h"
#include "acados_c/ocp_qp/ocp_qp_full_condensing.h"
#include "acados_c/dense_qp.h"



// void *ocp_qp_full_condensing_solver_copy_args(ocp_qp_solver_config *config, ocp_qp_dims *dims, void *raw_memory, void *source_)
// {
//     ocp_qp_full_condensing_solver_args *source = (ocp_qp_full_condensing_solver_args *) source_;
//     ocp_qp_full_condensing_solver_args *dest;

//     dest = ocp_qp_assign_args(config, dims, raw_memory);

//     ocp_qp_full_condensing_copy_args(config, dims, dest->cond_args, source->cond_args);

//     dense_qp_dims ddims;
//     compute_dense_qp_dims(dims, &ddims);

//     ocp_qp_solver_t solver_name = config->qp_solver;

//     dense_qp_solver_config solver_config;
//     switch(solver_name) {
//         case FULL_CONDENSING_HPIPM:
//             solver_config.qp_solver = DENSE_QP_HPIPM;
//             dense_qp_hpipm_copy_args(&solver_config, &ddims, dest->solver_args, source->solver_args);
//             break;
//         case FULL_CONDENSING_QPOASES:
//             solver_config.qp_solver = DENSE_QP_QPOASES;
//             dense_qp_qpoases_copy_args(&solver_config, &ddims, dest->solver_args, source->solver_args);
//             break;
//         case FULL_CONDENSING_QORE:
//             solver_config.qp_solver = DENSE_QP_QORE;
//             dense_qp_qore_copy_args(&solver_config, &ddims, dest->solver_args, source->solver_args);
//             break;
//         default:
//             printf(
//                 "\n\nSpecified solver is does not employ full condensing!\n\n");
//             exit(1);
//     }

//     return (void*)dest;
// }



int ocp_qp_full_condensing_solver_calculate_submodules_size(ocp_qp_solver_config *config, ocp_qp_dims *dims)
{
    int size = sizeof(ocp_qp_full_condensing_solver_submodules);

    ocp_qp_solver_t solver_name = config->qp_solver;

    dense_qp_solver_config subconfig;
    switch (solver_name) {
        case FULL_CONDENSING_HPIPM:
            subconfig.qp_solver = DENSE_QP_HPIPM;
            break;
        case FULL_CONDENSING_QPOASES:
            subconfig.qp_solver = DENSE_QP_QPOASES;
            break;
        case FULL_CONDENSING_QORE:
            subconfig.qp_solver = DENSE_QP_QORE;
            break;
    }

    dense_qp_dims ddims;
    compute_dense_qp_dims(dims, &ddims);

    size += calculate_dense_qp_solver_fcn_ptrs_size(&subconfig, &ddims);

    return size;
}



void *ocp_qp_full_condensing_solver_assign_submodules(ocp_qp_solver_config *config, ocp_qp_dims *dims, void *raw_memory)
{
    ocp_qp_full_condensing_solver_submodules *submodules;

    ocp_qp_solver_t solver_name = config->qp_solver;

    dense_qp_solver_config subconfig;
    switch (solver_name) {
        case FULL_CONDENSING_HPIPM:
            subconfig.qp_solver = DENSE_QP_HPIPM;
            break;
        case FULL_CONDENSING_QPOASES:
            subconfig.qp_solver = DENSE_QP_QPOASES;
            break;
        case FULL_CONDENSING_QORE:
            subconfig.qp_solver = DENSE_QP_QORE;
            break;
    }

    char *c_ptr = (char *) raw_memory;

    submodules = (ocp_qp_full_condensing_solver_submodules *) c_ptr;
    c_ptr += sizeof(ocp_qp_full_condensing_solver_submodules);

    dense_qp_dims ddims;
    compute_dense_qp_dims(dims, &ddims);

    submodules->solver = assign_dense_qp_solver_fcn_ptrs(&subconfig, &ddims, c_ptr);
    c_ptr += calculate_dense_qp_solver_fcn_ptrs_size(&subconfig, &ddims);

    assert((size_t)c_ptr % 8 == 0 && "memory not 8-byte aligned!");

    assert((char*)raw_memory + ocp_qp_full_condensing_solver_calculate_submodules_size(config, dims) == c_ptr);

    return (void *)submodules;
}
