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
#include <string.h>
#include <stdlib.h>

// acados
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/ocp_qp/ocp_qp_partial_condensing.h"
#include "acados/ocp_qp/ocp_qp_partial_condensing_solver.h"
#include "acados/utils/mem.h"
#include "acados/utils/timing.h"
#include "acados/utils/types.h"

/************************************************
 * opts
 ************************************************/

int ocp_qp_partial_condensing_solver_opts_calculate_size(void *config_, ocp_qp_dims *dims)
{
    ocp_qp_xcond_solver_config *config = config_;
    qp_solver_config *qp_solver = config->qp_solver;

    int size = 0;
    size += sizeof(ocp_qp_partial_condensing_solver_opts);

    size += ocp_qp_partial_condensing_opts_calculate_size(dims);
    size += qp_solver->opts_calculate_size(qp_solver, dims);

    return size;
}

void *ocp_qp_partial_condensing_solver_opts_assign(void *config_, ocp_qp_dims *dims,
                                                   void *raw_memory)
{
    ocp_qp_xcond_solver_config *config = config_;
    qp_solver_config *qp_solver = config->qp_solver;

    char *c_ptr = (char *) raw_memory;

    ocp_qp_partial_condensing_solver_opts *opts = (ocp_qp_partial_condensing_solver_opts *) c_ptr;
    c_ptr += sizeof(ocp_qp_partial_condensing_solver_opts);

    assert((size_t) c_ptr % 8 == 0 && "double not 8-byte aligned!");

    opts->pcond_opts = ocp_qp_partial_condensing_opts_assign(dims, c_ptr);
    c_ptr += ocp_qp_partial_condensing_opts_calculate_size(dims);

    assert((size_t) c_ptr % 8 == 0 && "double not 8-byte aligned!");

    opts->qp_solver_opts = qp_solver->opts_assign(qp_solver, dims, c_ptr);
    c_ptr += qp_solver->opts_calculate_size(qp_solver, dims);

    assert((char *) raw_memory +
               ocp_qp_partial_condensing_solver_opts_calculate_size(config_, dims) ==
           c_ptr);

    return (void *) opts;
}

void ocp_qp_partial_condensing_solver_opts_initialize_default(void *config_, ocp_qp_dims *dims,
                                                              void *opts_)
{
    ocp_qp_xcond_solver_config *config = config_;
    qp_solver_config *qp_solver = config->qp_solver;

    // partial_condensing solver opts
    ocp_qp_partial_condensing_solver_opts *opts = (ocp_qp_partial_condensing_solver_opts *) opts_;
    // partial condensing opts
    ocp_qp_partial_condensing_opts_initialize_default(dims, opts->pcond_opts);
    // ocp_qp_partial_condensing_opts *pcond_opts = opts->pcond_opts;
    // pcond_opts->N2 = config->N2;
    // qp solver opts
    qp_solver->opts_initialize_default(qp_solver, dims, opts->qp_solver_opts);
}

void ocp_qp_partial_condensing_solver_opts_update(void *config_, ocp_qp_dims *dims, void *opts_)
{
    ocp_qp_xcond_solver_config *config = config_;
    qp_solver_config *qp_solver = config->qp_solver;

    // partial_condensing solver opts
    ocp_qp_partial_condensing_solver_opts *opts = (ocp_qp_partial_condensing_solver_opts *) opts_;
    // partial condensing opts
    ocp_qp_partial_condensing_opts_update(dims, opts->pcond_opts);
    // ocp_qp_partial_condensing_opts *pcond_opts = opts->pcond_opts;
    // pcond_opts->N2 = config->N2;
    // qp solver opts
    qp_solver->opts_update(qp_solver, dims, opts->qp_solver_opts);
}


void ocp_qp_partial_condensing_solver_opts_set(void *config_, void *opts_,
                                            const char *field, const void* value)
{
    ocp_qp_partial_condensing_solver_opts *opts = (ocp_qp_partial_condensing_solver_opts *) opts_;
    // ocp_qp_xcond_solver_config *config = config_;

    if (!strcmp(field, "N2") || !strcmp(field, "pcond_N2"))
    {
        int* N2 = (int *) value;
        opts->pcond_opts->N2 = *N2;
    }
    else if (!strcmp(field, "N2_bkp") || !strcmp(field, "pcond_N2_bkp"))
    {
        int* N2_bkp = (int *) value;
        opts->pcond_opts->N2_bkp = *N2_bkp;
    }
    else
    {
        printf("\nerror: option type %s not available in ocp_qp_partial_condense solver module\n",
               field);
        exit(1);
    }
}


/************************************************
 * memory
 ************************************************/

int ocp_qp_partial_condensing_solver_memory_calculate_size(void *config_, ocp_qp_dims *dims,
                                                           void *opts_)
{
    ocp_qp_xcond_solver_config *config = config_;
    qp_solver_config *qp_solver = config->qp_solver;

    ocp_qp_partial_condensing_solver_opts *opts = (ocp_qp_partial_condensing_solver_opts *) opts_;

    int size = 0;
    size += sizeof(ocp_qp_partial_condensing_solver_memory);

    // set up dimesions of partially condensed qp
    ocp_qp_dims *pcond_dims;
    if (opts->pcond_opts->N2 < dims->N)
    {
        pcond_dims = opts->pcond_opts->pcond_dims;
    }
    else
    {
        pcond_dims = dims;
    }

    if (opts->pcond_opts->N2 < dims->N)
    {
        size += ocp_qp_partial_condensing_memory_calculate_size(dims, opts->pcond_opts);
    }

    size += qp_solver->memory_calculate_size(qp_solver, pcond_dims, opts->qp_solver_opts);

    if (opts->pcond_opts->N2 < dims->N)
    {
        size += ocp_qp_in_calculate_size(qp_solver, pcond_dims);
        size += ocp_qp_out_calculate_size(qp_solver, pcond_dims);
    }

    return size;
}

void *ocp_qp_partial_condensing_solver_memory_assign(void *config_, ocp_qp_dims *dims, void *opts_,
                                                     void *raw_memory)
{
    ocp_qp_xcond_solver_config *config = config_;
    qp_solver_config *qp_solver = config->qp_solver;

    ocp_qp_partial_condensing_solver_opts *opts = (ocp_qp_partial_condensing_solver_opts *) opts_;

    assert(opts->pcond_opts->N2 > 0 && "N2 must be positive!");
    assert(opts->pcond_opts->N2 <= dims->N && "N2 cannot be bigger than N!");

    char *c_ptr = (char *) raw_memory;

    // set up dimesions of partially condensed qp
    ocp_qp_dims *pcond_dims;
    if (opts->pcond_opts->N2 < dims->N)
    {
        pcond_dims = opts->pcond_opts->pcond_dims;
    }
    else
    {
        pcond_dims = dims;
    }

    ocp_qp_partial_condensing_solver_memory *mem =
        (ocp_qp_partial_condensing_solver_memory *) c_ptr;
    c_ptr += sizeof(ocp_qp_partial_condensing_solver_memory);

    assert((size_t) c_ptr % 8 == 0 && "double not 8-byte aligned!");

    if (opts->pcond_opts->N2 < dims->N)
    {
        mem->pcond_memory =
            (ocp_qp_partial_condensing_memory *) ocp_qp_partial_condensing_memory_assign(
                dims, opts->pcond_opts, c_ptr);
        c_ptr += ocp_qp_partial_condensing_memory_calculate_size(dims, opts->pcond_opts);
    }
    else
    {
        mem->pcond_memory = NULL;
    }

    assert((size_t) c_ptr % 8 == 0 && "double not 8-byte aligned!");

    mem->solver_memory =
        qp_solver->memory_assign(qp_solver, pcond_dims, opts->qp_solver_opts, c_ptr);
    c_ptr += qp_solver->memory_calculate_size(qp_solver, pcond_dims, opts->qp_solver_opts);

    if (opts->pcond_opts->N2 < dims->N)
    {
        mem->pcond_qp_in = ocp_qp_in_assign(qp_solver, pcond_dims, c_ptr);
        c_ptr += ocp_qp_in_calculate_size(qp_solver, pcond_dims);
    }
    else
    {
        mem->pcond_qp_in = NULL;
    }

    if (opts->pcond_opts->N2 < dims->N)
    {
        mem->pcond_qp_out = ocp_qp_out_assign(qp_solver, pcond_dims, c_ptr);
        c_ptr += ocp_qp_out_calculate_size(qp_solver, pcond_dims);
    }
    else
    {
        mem->pcond_qp_out = NULL;
    }

    assert((char *) raw_memory +
               ocp_qp_partial_condensing_solver_memory_calculate_size(config_, dims, opts_) ==
           c_ptr);

    return mem;
}

/************************************************
 * workspace
 ************************************************/

int ocp_qp_partial_condensing_solver_workspace_calculate_size(void *config_, ocp_qp_dims *dims,
                                                              void *opts_)
{
    ocp_qp_xcond_solver_config *config = config_;
    qp_solver_config *qp_solver = config->qp_solver;

    ocp_qp_partial_condensing_solver_opts *opts = (ocp_qp_partial_condensing_solver_opts *) opts_;

    int size = sizeof(ocp_qp_partial_condensing_solver_workspace);
    size += ocp_qp_partial_condensing_workspace_calculate_size(dims, opts->pcond_opts);

    // set up dimesions of partially condensed qp
    ocp_qp_dims *pcond_dims;
    if (opts->pcond_opts->N2 < dims->N)
    {
        pcond_dims = opts->pcond_opts->pcond_dims;
    }
    else
    {
        pcond_dims = dims;
    }

    size += qp_solver->workspace_calculate_size(qp_solver, pcond_dims, opts->qp_solver_opts);

    return size;
}

static void cast_workspace(void *config_, ocp_qp_dims *dims,
                           ocp_qp_partial_condensing_solver_opts *opts,
                           ocp_qp_partial_condensing_solver_memory *mem,
                           ocp_qp_partial_condensing_solver_workspace *work)
{
    ocp_qp_xcond_solver_config *config = config_;
    qp_solver_config *qp_solver = config->qp_solver;

    // set up dimesions of partially condensed qp
    ocp_qp_dims *pcond_dims;
    if (opts->pcond_opts->N2 < dims->N)
    {
        pcond_dims = opts->pcond_opts->pcond_dims;
    }
    else
    {
        pcond_dims = dims;
    }

    char *c_ptr = (char *) work;

    c_ptr += sizeof(ocp_qp_partial_condensing_solver_workspace);

    work->pcond_work = c_ptr;
    c_ptr += ocp_qp_partial_condensing_workspace_calculate_size(dims, opts->pcond_opts);

    work->solver_work = c_ptr;
    c_ptr += qp_solver->workspace_calculate_size(qp_solver, pcond_dims, opts->qp_solver_opts);
}

/************************************************
 * functions
 ************************************************/

int ocp_qp_partial_condensing_solver(void *config_, ocp_qp_in *qp_in, ocp_qp_out *qp_out,
                                     void *opts_, void *mem_, void *work_)
{
    ocp_qp_xcond_solver_config *config = config_;
    qp_solver_config *qp_solver = config->qp_solver;

    ocp_qp_info *info = (ocp_qp_info *) qp_out->misc;
    acados_timer tot_timer, cond_timer;
    acados_tic(&tot_timer);

    // cast data structures
    ocp_qp_partial_condensing_solver_opts *opts = (ocp_qp_partial_condensing_solver_opts *) opts_;
    ocp_qp_partial_condensing_solver_memory *memory =
        (ocp_qp_partial_condensing_solver_memory *) mem_;
    ocp_qp_partial_condensing_solver_workspace *work =
        (ocp_qp_partial_condensing_solver_workspace *) work_;

    // cast workspace
    cast_workspace(config_, qp_in->dim, opts, memory, work);

    if (opts->pcond_opts->N2 < qp_in->dim->N)
    {  // condensing
        acados_tic(&cond_timer);
        ocp_qp_partial_condensing(qp_in, memory->pcond_qp_in, opts->pcond_opts,
                                  memory->pcond_memory, work->pcond_work);
        info->condensing_time = acados_toc(&cond_timer);
    }
    else
    {
        memory->pcond_qp_in = qp_in;
        memory->pcond_qp_out = qp_out;
        info->condensing_time = 0;
    }

    // solve qp
    int solver_status;
    if (memory->solver_memory != NULL)
    {
        solver_status =
            qp_solver->evaluate(qp_solver, memory->pcond_qp_in, memory->pcond_qp_out,
                                opts->qp_solver_opts, memory->solver_memory, work->solver_work);
    }
    else
    {
        solver_status = ACADOS_FAILURE;
    }

    // expand
    if (opts->pcond_opts->N2 < qp_in->dim->N)
    {
        acados_tic(&cond_timer);
        ocp_qp_partial_expansion(memory->pcond_qp_out, qp_out, opts->pcond_opts,
                                 memory->pcond_memory, work->pcond_work);
        info->condensing_time += acados_toc(&cond_timer);
    }

    info->total_time = acados_toc(&tot_timer);
    info->solve_QP_time = ((ocp_qp_info *) (memory->pcond_qp_out->misc))->solve_QP_time;
    info->interface_time = ((ocp_qp_info *) (memory->pcond_qp_out->misc))->interface_time;
    info->num_iter = ((ocp_qp_info *) (memory->pcond_qp_out->misc))->num_iter;
    info->t_computed = ((ocp_qp_info *) (memory->pcond_qp_out->misc))->t_computed;

    return solver_status;
}

void ocp_qp_partial_condensing_solver_config_initialize_default(void *config_)
{
    ocp_qp_xcond_solver_config *config = config_;

    config->dims_set = &ocp_qp_dims_set;
    config->opts_calculate_size = &ocp_qp_partial_condensing_solver_opts_calculate_size;
    config->opts_assign = &ocp_qp_partial_condensing_solver_opts_assign;
    config->opts_initialize_default = &ocp_qp_partial_condensing_solver_opts_initialize_default;
    config->opts_update = &ocp_qp_partial_condensing_solver_opts_update;
    config->opts_set = &ocp_qp_partial_condensing_solver_opts_set;
    config->memory_calculate_size = &ocp_qp_partial_condensing_solver_memory_calculate_size;
    config->memory_assign = &ocp_qp_partial_condensing_solver_memory_assign;
    config->workspace_calculate_size = &ocp_qp_partial_condensing_solver_workspace_calculate_size;
    config->evaluate = &ocp_qp_partial_condensing_solver;

    return;
}
