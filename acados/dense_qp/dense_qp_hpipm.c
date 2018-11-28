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

// external
#include <assert.h>
// hpipm
#include "hpipm/include/hpipm_d_dense_qp.h"
#include "hpipm/include/hpipm_d_dense_qp_ipm.h"
#include "hpipm/include/hpipm_d_dense_qp_sol.h"
// acados
#include "acados/dense_qp/dense_qp_common.h"
#include "acados/dense_qp/dense_qp_hpipm.h"
#include "acados/utils/mem.h"
#include "acados/utils/timing.h"

/************************************************
 * opts
 ************************************************/

int dense_qp_hpipm_opts_calculate_size(void *config_, void *dims_)
{
    dense_qp_dims *dims = dims_;

    int size = 0;
    size += sizeof(dense_qp_hpipm_opts);
    size += sizeof(struct d_dense_qp_ipm_arg);
    size += d_memsize_dense_qp_ipm_arg(dims);

    return size;
}

void *dense_qp_hpipm_opts_assign(void *config_, void *dims_, void *raw_memory)
{
    dense_qp_dims *dims = dims_;
    dense_qp_hpipm_opts *opts;

    char *c_ptr = (char *) raw_memory;

    opts = (dense_qp_hpipm_opts *) c_ptr;
    c_ptr += sizeof(dense_qp_hpipm_opts);

    opts->hpipm_opts = (struct d_dense_qp_ipm_arg *) c_ptr;
    c_ptr += sizeof(struct d_dense_qp_ipm_arg);

    assert((size_t) c_ptr % 8 == 0 && "memory not 8-byte aligned!");

    d_create_dense_qp_ipm_arg(dims, opts->hpipm_opts, c_ptr);
    c_ptr += d_memsize_dense_qp_ipm_arg(dims);

    assert((char *) raw_memory + dense_qp_hpipm_opts_calculate_size(config_, dims) == c_ptr);

    return (void *) opts;
}

void dense_qp_hpipm_opts_initialize_default(void *config_, void *dims_, void *opts_)
{
    dense_qp_hpipm_opts *opts = opts_;

    d_set_default_dense_qp_ipm_arg(BALANCE, opts->hpipm_opts);
    // overwrite some default options
    opts->hpipm_opts->res_g_max = 1e-6;
    opts->hpipm_opts->res_b_max = 1e-8;
    opts->hpipm_opts->res_d_max = 1e-8;
    opts->hpipm_opts->res_m_max = 1e-8;
    opts->hpipm_opts->iter_max = 50;
    opts->hpipm_opts->stat_max = 50;
    opts->hpipm_opts->alpha_min = 1e-8;
    opts->hpipm_opts->mu0 = 1e0;

    return;
}

void dense_qp_hpipm_opts_update(void *config_, void *dims_, void *opts_)
{
    //    dense_qp_hpipm_opts *opts = (dense_qp_hpipm_opts *)opts_;

    return;
}

/************************************************
 * memory
 ************************************************/

int dense_qp_hpipm_memory_calculate_size(void *config_, void *dims_, void *opts_)
{
    dense_qp_dims *dims = dims_;
    dense_qp_hpipm_opts *opts = opts_;

    int size = 0;
    size += sizeof(dense_qp_hpipm_memory);
    size += sizeof(struct d_dense_qp_ipm_workspace);

    size += d_memsize_dense_qp_ipm(dims, opts->hpipm_opts);

    return size;
}

void *dense_qp_hpipm_memory_assign(void *config_, void *dims_, void *opts_, void *raw_memory)
{
    dense_qp_dims *dims = dims_;
    dense_qp_hpipm_opts *opts = opts_;
    dense_qp_hpipm_memory *mem;

    char *c_ptr = (char *) raw_memory;

    mem = (dense_qp_hpipm_memory *) c_ptr;
    c_ptr += sizeof(dense_qp_hpipm_memory);

    mem->hpipm_workspace = (struct d_dense_qp_ipm_workspace *) c_ptr;
    c_ptr += sizeof(struct d_dense_qp_ipm_workspace);

    struct d_dense_qp_ipm_workspace *ipm_workspace = mem->hpipm_workspace;

    assert((size_t) c_ptr % 8 == 0 && "memory not 8-byte aligned!");

    // ipm workspace structure
    d_create_dense_qp_ipm(dims, opts->hpipm_opts, ipm_workspace, c_ptr);
    c_ptr += ipm_workspace->memsize;

    assert((char *) raw_memory + dense_qp_hpipm_memory_calculate_size(config_, dims, opts) ==
           c_ptr);

    return mem;
}

int dense_qp_hpipm_workspace_calculate_size(void *config_, void *dims_, void *opts_) { return 0; }
int dense_qp_hpipm(void *config, void *qp_in_, void *qp_out_, void *opts_, void *mem_, void *work_)
{
    dense_qp_in *qp_in = qp_in_;
    dense_qp_out *qp_out = qp_out_;

    dense_qp_info *info = (dense_qp_info *) qp_out->misc;
    acados_timer tot_timer, qp_timer;

    acados_tic(&tot_timer);

    // cast structures
    dense_qp_hpipm_opts *opts = opts_;
    dense_qp_hpipm_memory *memory = mem_;

    // solve ipm
    acados_tic(&qp_timer);
    int hpipm_status =
        d_solve_dense_qp_ipm(qp_in, qp_out, opts->hpipm_opts, memory->hpipm_workspace);

    info->solve_QP_time = acados_toc(&qp_timer);
    info->interface_time = 0;  // there are no conversions for hpipm
    info->total_time = acados_toc(&tot_timer);
    info->num_iter = memory->hpipm_workspace->iter;
    info->t_computed = 1;

    // check exit conditions
    int acados_status = hpipm_status;
    if (hpipm_status == 0) acados_status = ACADOS_SUCCESS;
    if (hpipm_status == 1) acados_status = ACADOS_MAXITER;
    if (hpipm_status == 2) acados_status = ACADOS_MINSTEP;
    return acados_status;
}

void dense_qp_hpipm_config_initialize_default(void *config_)
{
    qp_solver_config *config = config_;

    config->dims_set = &dense_qp_dims_set;
    config->opts_calculate_size = &dense_qp_hpipm_opts_calculate_size;
    config->opts_assign = &dense_qp_hpipm_opts_assign;
    config->opts_initialize_default = &dense_qp_hpipm_opts_initialize_default;
    config->opts_update = &dense_qp_hpipm_opts_update;
    config->memory_calculate_size = &dense_qp_hpipm_memory_calculate_size;
    config->memory_assign = &dense_qp_hpipm_memory_assign;
    config->workspace_calculate_size = &dense_qp_hpipm_workspace_calculate_size;
    config->evaluate = &dense_qp_hpipm;

    return;
}
