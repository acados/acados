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
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

// acados
#include "acados/dense_qp/dense_qp_common.h"
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/ocp_qp/ocp_qp_full_condensing.h"
#include "acados/ocp_qp/ocp_qp_full_condensing_solver.h"
#include "acados/utils/mem.h"
#include "acados/utils/timing.h"
#include "acados/utils/types.h"

/************************************************
 * opts
 ************************************************/

int ocp_qp_full_condensing_solver_opts_calculate_size(void *config_, ocp_qp_dims *dims)
{
    ocp_qp_xcond_solver_config *config = config_;
    qp_solver_config *qp_solver = config->qp_solver;

    int size = 0;
    size += sizeof(ocp_qp_full_condensing_solver_opts);

    dense_qp_dims ddims;
    compute_dense_qp_dims(dims, &ddims);

    size += ocp_qp_full_condensing_opts_calculate_size(dims);

    // TODO(dimitris): shouldn't we pass config->qp_solver instead of config_ below?
    size += qp_solver->opts_calculate_size(config_, &ddims);

    return size;
}

void *ocp_qp_full_condensing_solver_opts_assign(void *config_, ocp_qp_dims *dims, void *raw_memory)
{
    ocp_qp_xcond_solver_config *config = config_;
    qp_solver_config *qp_solver = config->qp_solver;

    char *c_ptr = (char *) raw_memory;

    ocp_qp_full_condensing_solver_opts *opts = (ocp_qp_full_condensing_solver_opts *) c_ptr;
    c_ptr += sizeof(ocp_qp_full_condensing_solver_opts);

    dense_qp_dims ddims;
    compute_dense_qp_dims(dims, &ddims);

    assert((size_t) c_ptr % 8 == 0 && "memory not 8-byte aligned!");

    opts->cond_opts = ocp_qp_full_condensing_opts_assign(dims, c_ptr);
    c_ptr += ocp_qp_full_condensing_opts_calculate_size(dims);

    align_char_to(8, &c_ptr);

    opts->qp_solver_opts = qp_solver->opts_assign(qp_solver, &ddims, c_ptr);
    c_ptr += qp_solver->opts_calculate_size(qp_solver, &ddims);

    assert((char *) raw_memory + ocp_qp_full_condensing_solver_opts_calculate_size(config_, dims) ==
           c_ptr);

    return (void *) opts;
}

void ocp_qp_full_condensing_solver_opts_initialize_default(void *config_, ocp_qp_dims *dims,
                                                           void *opts_)
{
    ocp_qp_xcond_solver_config *config = config_;
    qp_solver_config *qp_solver = config->qp_solver;

    // full cond solver
    ocp_qp_full_condensing_solver_opts *opts = (ocp_qp_full_condensing_solver_opts *) opts_;
    // full condensing
    ocp_qp_full_condensing_opts_initialize_default(dims, opts->cond_opts);
    // qp solver
    qp_solver->opts_initialize_default(qp_solver, NULL,
                                       opts->qp_solver_opts);  // TODO(all): pass dense_qp_dims ???
}

void ocp_qp_full_condensing_solver_opts_update(void *config_, ocp_qp_dims *dims, void *opts_)
{
    ocp_qp_xcond_solver_config *config = config_;
    qp_solver_config *qp_solver = config->qp_solver;

    // full cond solver
    ocp_qp_full_condensing_solver_opts *opts = (ocp_qp_full_condensing_solver_opts *) opts_;
    // full condensing
    ocp_qp_full_condensing_opts_update(dims, opts->cond_opts);
    // qp solver
    qp_solver->opts_update(qp_solver, NULL, opts->qp_solver_opts);  // TODO(all): pass dense_qp_dims
}

void ocp_qp_full_condensing_solver_opts_set(void *config_, void *opts_,
                                            const char *field, const void* value)
{
    ocp_qp_full_condensing_solver_opts *opts = (ocp_qp_full_condensing_solver_opts *) opts_;
    // ocp_qp_xcond_solver_config *config = config_;

    if (!strcmp(field, "condense_rhs_only"))
    {
        int* condense_rhs_only = (int *) value;
        opts->cond_opts->condense_rhs_only = *condense_rhs_only;
    }
    else if (!strcmp(field, "expand_primal_sol_only"))
    {
        int* expand_primal_sol_only = (int *) value;
        opts->cond_opts->expand_primal_sol_only = *expand_primal_sol_only;
    }
    else
    {
        printf("\nerror: option type %s not available in ocp_qp_full_condensing_solver module\n",
               field);
        exit(1);
    }
}



/************************************************
 * memory
 ************************************************/

int ocp_qp_full_condensing_solver_memory_calculate_size(void *config_, ocp_qp_dims *dims,
                                                        void *opts_)
{
    ocp_qp_xcond_solver_config *config = config_;
    qp_solver_config *qp_solver = config->qp_solver;

    ocp_qp_full_condensing_solver_opts *opts = (ocp_qp_full_condensing_solver_opts *) opts_;

    dense_qp_dims ddims;
    compute_dense_qp_dims(dims, &ddims);

    int size = 0;
    size += sizeof(ocp_qp_full_condensing_solver_memory);

    size += ocp_qp_full_condensing_memory_calculate_size(dims, opts->cond_opts);
    size += qp_solver->memory_calculate_size(qp_solver, &ddims, opts->qp_solver_opts);

    size += dense_qp_in_calculate_size(qp_solver, &ddims);
    size += dense_qp_out_calculate_size(qp_solver, &ddims);

    return size;
}

void *ocp_qp_full_condensing_solver_memory_assign(void *config_, ocp_qp_dims *dims, void *opts_,
                                                  void *raw_memory)
{
    ocp_qp_xcond_solver_config *config = config_;
    qp_solver_config *qp_solver = config->qp_solver;

    ocp_qp_full_condensing_solver_opts *opts = (ocp_qp_full_condensing_solver_opts *) opts_;

    dense_qp_dims ddims;
    compute_dense_qp_dims(dims, &ddims);

    char *c_ptr = (char *) raw_memory;

    ocp_qp_full_condensing_solver_memory *mem = (ocp_qp_full_condensing_solver_memory *) c_ptr;
    c_ptr += sizeof(ocp_qp_full_condensing_solver_memory);

    assert((size_t) c_ptr % 8 == 0 && "memory not 8-byte aligned!");

    mem->cond_memory = (ocp_qp_full_condensing_memory *) ocp_qp_full_condensing_memory_assign(
        dims, opts->cond_opts, c_ptr);
    c_ptr += ocp_qp_full_condensing_memory_calculate_size(dims, opts->cond_opts);

    assert((size_t) c_ptr % 8 == 0 && "memory not 8-byte aligned!");

    mem->solver_memory = qp_solver->memory_assign(qp_solver, &ddims, opts->qp_solver_opts, c_ptr);
    c_ptr += qp_solver->memory_calculate_size(qp_solver, &ddims, opts->qp_solver_opts);

    assert((size_t) c_ptr % 8 == 0 && "memory not 8-byte aligned!");

    mem->qpd_in = dense_qp_in_assign(qp_solver, &ddims, c_ptr);
    c_ptr += dense_qp_in_calculate_size(qp_solver, &ddims);

    assert((size_t) c_ptr % 8 == 0 && "memory not 8-byte aligned!");

    mem->qpd_out = dense_qp_out_assign(qp_solver, &ddims, c_ptr);
    c_ptr += dense_qp_out_calculate_size(qp_solver, &ddims);

    assert((char *) raw_memory +
               ocp_qp_full_condensing_solver_memory_calculate_size(config_, dims, opts_) ==
           c_ptr);

    return mem;
}

/************************************************
 * workspace
 ************************************************/

int ocp_qp_full_condensing_solver_workspace_calculate_size(void *config_, ocp_qp_dims *dims,
                                                           void *opts_)
{
    ocp_qp_xcond_solver_config *config = config_;
    qp_solver_config *qp_solver = config->qp_solver;

    ocp_qp_full_condensing_solver_opts *opts = (ocp_qp_full_condensing_solver_opts *) opts_;

    dense_qp_dims ddims;
    compute_dense_qp_dims(dims, &ddims);

    int size = sizeof(ocp_qp_full_condensing_solver_workspace);
    size += ocp_qp_full_condensing_workspace_calculate_size(dims, opts->cond_opts);
    size += qp_solver->workspace_calculate_size(qp_solver, &ddims, opts->qp_solver_opts);

    return size;
}

static void cast_workspace(void *config_, ocp_qp_dims *dims,
                           ocp_qp_full_condensing_solver_opts *opts,
                           ocp_qp_full_condensing_solver_memory *mem,
                           ocp_qp_full_condensing_solver_workspace *work)
{
    ocp_qp_xcond_solver_config *config = config_;
    qp_solver_config *qp_solver = config->qp_solver;

    dense_qp_dims *ddims = mem->qpd_in->dim;

    char *c_ptr = (char *) work;

    c_ptr += sizeof(ocp_qp_full_condensing_solver_workspace);

    work->cond_work = c_ptr;
    c_ptr += ocp_qp_full_condensing_workspace_calculate_size(dims, opts->cond_opts);

    work->solver_workspace = c_ptr;
    c_ptr += qp_solver->workspace_calculate_size(qp_solver, ddims, opts->qp_solver_opts);

    assert((char *) work +
               ocp_qp_full_condensing_solver_workspace_calculate_size(config_, dims, opts) >=
           c_ptr);
}

/************************************************
 * functions
 ************************************************/

int ocp_qp_full_condensing_solver(void *config_, ocp_qp_in *qp_in, ocp_qp_out *qp_out, void *opts_,
                                  void *mem_, void *work_)
{
    ocp_qp_xcond_solver_config *config = config_;
    qp_solver_config *qp_solver = config->qp_solver;

    ocp_qp_info *info = (ocp_qp_info *) qp_out->misc;
    acados_timer tot_timer, cond_timer;
    acados_tic(&tot_timer);

    // cast data structures
    ocp_qp_full_condensing_solver_opts *opts = (ocp_qp_full_condensing_solver_opts *) opts_;
    ocp_qp_full_condensing_solver_memory *memory = (ocp_qp_full_condensing_solver_memory *) mem_;
    ocp_qp_full_condensing_solver_workspace *work =
        (ocp_qp_full_condensing_solver_workspace *) work_;

    // cast workspace
    cast_workspace(config_, qp_in->dim, opts, memory, work);

    // condense
    acados_tic(&cond_timer);
    ocp_qp_full_condensing(qp_in, memory->qpd_in, opts->cond_opts, memory->cond_memory,
                           work->cond_work);
    info->condensing_time = acados_toc(&cond_timer);

    // solve qp
    int solver_status =
        qp_solver->evaluate(qp_solver, memory->qpd_in, memory->qpd_out, opts->qp_solver_opts,
                            memory->solver_memory, work->solver_workspace);

    // expand
    acados_tic(&cond_timer);
    ocp_qp_full_expansion(memory->qpd_out, qp_out, opts->cond_opts, memory->cond_memory,
                          work->cond_work);
    info->condensing_time += acados_toc(&cond_timer);

    info->total_time = acados_toc(&tot_timer);
    info->solve_QP_time = ((dense_qp_info *) (memory->qpd_out->misc))->solve_QP_time;
    info->interface_time = ((dense_qp_info *) (memory->qpd_out->misc))->interface_time;
    info->num_iter = ((dense_qp_info *) (memory->qpd_out->misc))->num_iter;
    info->t_computed = ((dense_qp_info *) (memory->qpd_out->misc))->t_computed;

    return solver_status;
}

void ocp_qp_full_condensing_solver_config_initialize_default(void *config_)
{
    ocp_qp_xcond_solver_config *config = config_;

    config->dims_set = &ocp_qp_dims_set;
    config->opts_calculate_size = &ocp_qp_full_condensing_solver_opts_calculate_size;
    config->opts_assign = &ocp_qp_full_condensing_solver_opts_assign;
    config->opts_set = &ocp_qp_full_condensing_solver_opts_set;
    config->opts_initialize_default = &ocp_qp_full_condensing_solver_opts_initialize_default;
    config->opts_update = &ocp_qp_full_condensing_solver_opts_update;
    config->memory_calculate_size = &ocp_qp_full_condensing_solver_memory_calculate_size;
    config->memory_assign = &ocp_qp_full_condensing_solver_memory_assign;
    config->workspace_calculate_size = &ocp_qp_full_condensing_solver_workspace_calculate_size;
    config->evaluate = &ocp_qp_full_condensing_solver;

    return;
}
