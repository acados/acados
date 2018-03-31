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

#include "acados_c/condensing_interface.h"

//external
#include <stdlib.h>
#include <assert.h>
#include <string.h>

// acados
#include "acados/ocp_qp/ocp_qp_partial_condensing.h"
#include "acados/ocp_qp/ocp_qp_full_condensing.h"
#include "acados/utils/mem.h"

condensing_config *condensing_config_create(condensing_plan *plan)
{
    // int bytes = condensing_module_config_calculate_size();
    // void *ptr = calloc(1, bytes);
    // condensing_config *solver_config = dense_condensing_config_assign(ptr);

    // dense_qp_solver_t solver_name = plan->qp_solver;

    // // TODO(dimitris): cath error if solver not compiled
    // // printf("\n\nSpecified solver interface not compiled with acados!\n\n");
    // switch (solver_name)
    // {
    //     case DENSE_QP_HPIPM:
    //         dense_qp_hpipm_config_initialize_default(solver_config);
    //         break;
    //     case DENSE_QP_QPOASES:
    //         #ifdef ACADOS_WITH_QPOASES
    //         dense_qp_qpoases_config_initialize_default(solver_config);
    //         #endif
    //         break;
    //     case DENSE_QP_QORE:
    //         #ifdef ACADOS_WITH_QORE
    //         dense_qp_qore_config_initialize_default(solver_config);
    //         #endif
    //         break;
    // }
    // return condensing_config;
}



void *condensing_opts_create(condensing_config *config, void *dims_)
{
    // int bytes = config->opts_calculate_size(config, dims);

    // void *ptr = calloc(1, bytes);

    // void *opts = config->opts_assign(config, dims, ptr);

    // config->opts_initialize_default(config, dims, opts);

    // return opts;
}



int condensing_calculate_size(condensing_config *config, void *dims_, void *opts_)
{
    // int bytes = sizeof(dense_qp_solver);

    // bytes += config->memory_calculate_size(config, dims, opts_);
    // bytes += config->workspace_calculate_size(config, dims, opts_);

    // return bytes;
}



condensing_module *condensing_assign(condensing_config *config, void *dims_, void *opts_, void *raw_memory)
{
    // char *c_ptr = (char *) raw_memory;

    // dense_qp_solver *solver = (dense_qp_solver *) c_ptr;
    // c_ptr += sizeof(dense_qp_solver);

    // solver->config = config;
    // solver->dims = dims;
    // solver->opts = opts_;

    // // TODO(dimitris): CHECK ALIGNMENT!

    // solver->mem = config->memory_assign(config, dims, opts_, c_ptr);
    // c_ptr += config->memory_calculate_size(config, dims, opts_);

    // solver->work = (void *) c_ptr;
    // c_ptr += config->workspace_calculate_size(config, dims, opts_);

    // assert((char*)raw_memory + dense_qp_calculate_size(config, dims, opts_) == c_ptr);

    // return solver;
}



condensing_module *condensing_create(condensing_config *config, void *dims_, void *opts_)
{
    // int bytes = dense_qp_calculate_size(config, dims, opts_);

    // void *ptr = calloc(1, bytes);

    // dense_qp_solver *solver = dense_qp_assign(config, dims, opts_, ptr);

    // return solver;
}



int condense(condensing_module *module, void *qp_in, void *qp_out)
{
    return solver->config->evaluate(solver->config, qp_in, qp_out, solver->opts, solver->mem, solver->work);
}
