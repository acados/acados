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

#include <assert.h>

// osqp
#include "osqp/include/types.h"
#include "osqp/include/osqp.h"

// acados
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/ocp_qp/ocp_qp_osqp.h"
#include "acados/utils/mem.h"
#include "acados/utils/timing.h"
#include "acados/utils/types.h"

/************************************************
 * opts
 ************************************************/

int ocp_qp_osqp_opts_calculate_size(void *config_, void *dims_)
{
    int size = 0;
    size += sizeof(ocp_qp_osqp_opts);
    size += sizeof(OSQPSettings);

    return size;
}



void *ocp_qp_osqp_opts_assign(void *config_, void *dims_, void *raw_memory)
{
    ocp_qp_osqp_opts *opts;

    char *c_ptr = (char *) raw_memory;

    opts = (ocp_qp_osqp_opts *) c_ptr;
    c_ptr += sizeof(ocp_qp_osqp_opts);

    opts->osqp_opts = (OSQPSettings *) c_ptr;
    c_ptr += sizeof(OSQPSettings);

    assert((char *) raw_memory + ocp_qp_osqp_opts_calculate_size(config_, dims) == c_ptr);

    return (void *) opts;
}



void ocp_qp_osqp_opts_initialize_default(void *config_, void *dims_, void *opts_)
{
    ocp_qp_osqp_opts *opts = opts_;

    osqp_set_default_settings(opts->osqp_opts);

    return;
}



void ocp_qp_osqp_opts_update(void *config_, void *dims_, void *opts_)
{
    // ocp_qp_osqp_opts *opts = (ocp_qp_osqp_opts *)opts_;

    return;
}

/************************************************
 * memory
 ************************************************/

int ocp_qp_osqp_memory_calculate_size(void *config_, void *dims_, void *opts_)
{
    int size = 0;
    size += sizeof(ocp_qp_osqp_memory);

    size += sizeof(OSQPData);

    return size;
}



void *ocp_qp_osqp_memory_assign(void *config_, void *dims_, void *opts_, void *raw_memory)
{
    ocp_qp_dims *dims = dims_;
    ocp_qp_osqp_opts *opts = opts_;
    ocp_qp_osqp_memory *mem;

    // char pointer
    char *c_ptr = (char *) raw_memory;

    mem = (ocp_qp_osqp_memory *) c_ptr;
    c_ptr += sizeof(ocp_qp_osqp_memory);

    mem->osqp_data = (OSQPData *) c_ptr;
    c_ptr += sizeof(OSQPData);

    mem->osqp_work = (OSQPWorkspace *) c_ptr;
    c_ptr += sizeof(OSQPWorkspace);

    // TODO(dimitris): implement
    // mem->osqp_data->n = n;
    // mem->osqp_data->m = m;
    // mem->osqp_data->P = csc_matrix(data->n, data->n, P_nnz, P_x, P_i, P_p);
    // mem->osqp_data->q = q;
    // mem->osqp_data->A = csc_matrix(data->m, data->n, A_nnz, A_x, A_i, A_p);
    // mem->osqp_data->l = l;
    // mem->osqp_data->u = u;

    mem->osqp_work = osqp_setup(mem->osqp_data, opts->osqp_opts);

    assert((char *)raw_memory + ocp_qp_osqp_memory_calculate_size(config_, dims, opts_) == c_ptr);

    return mem;
}

/************************************************
 * workspace
 ************************************************/

int ocp_qp_osqp_workspace_calculate_size(void *config_, void *dims_, void *opts_) { return 0; }

/************************************************
 * functions
 ************************************************/

int ocp_qp_osqp(void *config_, void *qp_in_, void *qp_out_, void *opts_, void *mem_, void *work_)
{
    ocp_qp_in *qp_in = qp_in_;
    ocp_qp_out *qp_out = qp_out_;

    ocp_qp_info *info = (ocp_qp_info *) qp_out->misc;
    acados_timer tot_timer, qp_timer, interface_timer;

    acados_tic(&tot_timer);
    // cast data structures
    ocp_qp_osqp_opts *opts = (ocp_qp_osqp_opts *) opts_;
    ocp_qp_osqp_memory *mem = (ocp_qp_osqp_memory *) mem_;

    acados_tic(&interface_timer);
    // TODO(dimitris): convert data to OSQP format
    info->interface_time = acados_toc(&interface_timer);

    acados_tic(&qp_timer);
    // TODO(dimitris): use update functions to insert new data in OSQP

    // solve OSQP
    osqp_solve(mem->osqp_work);

    info->solve_QP_time = acados_toc(&qp_timer);
    info->total_time = acados_toc(&tot_timer);
    info->num_iter = -1;  // TODO(dimitris): extract n_iter
    info->t_computed = 0;

    // // check exit conditions
    // int acados_status = osqp_status;
    // if (osqp_status == 0) acados_status = ACADOS_SUCCESS;
    // if (osqp_status == 1) acados_status = ACADOS_MAXITER;
    // if (osqp_status == 2) acados_status = ACADOS_MINSTEP;
    // return acados_status;
}



void ocp_qp_osqp_config_initialize_default(void *config_)
{
    qp_solver_config *config = config_;

    config->opts_calculate_size = &ocp_qp_osqp_opts_calculate_size;
    config->opts_assign = &ocp_qp_osqp_opts_assign;
    config->opts_initialize_default = &ocp_qp_osqp_opts_initialize_default;
    config->opts_update = &ocp_qp_osqp_opts_update;
    config->memory_calculate_size = &ocp_qp_osqp_memory_calculate_size;
    config->memory_assign = &ocp_qp_osqp_memory_assign;
    config->workspace_calculate_size = &ocp_qp_osqp_workspace_calculate_size;
    config->evaluate = &ocp_qp_osqp;

    return;
}
