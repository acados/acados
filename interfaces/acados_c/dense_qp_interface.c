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

#include "acados_c/dense_qp_interface.h"

// external
#include <assert.h>
#include <stdlib.h>
#include <string.h>
// acados_c

#include "acados/utils/mem.h"

#include "acados/dense_qp/dense_qp_hpipm.h"
#ifdef ACADOS_WITH_QORE
#include "acados/dense_qp/dense_qp_qore.h"
#endif
#include "acados/dense_qp/dense_qp_qpoases.h"

qp_solver_config *dense_qp_config_create(dense_qp_solver_plan *plan)
{
    int bytes = dense_qp_solver_config_calculate_size();
    void *ptr = calloc(1, bytes);
    qp_solver_config *solver_config = dense_qp_solver_config_assign(ptr);

    dense_qp_solver_t solver_name = plan->qp_solver;

    // TODO(dimitris): cath error if solver not compiled
    // printf("\n\nSpecified solver interface not compiled with acados!\n\n");
    switch (solver_name)
    {
        case DENSE_QP_HPIPM:
            dense_qp_hpipm_config_initialize_default(solver_config);
            break;
        case DENSE_QP_QPOASES:
#ifdef ACADOS_WITH_QPOASES
            dense_qp_qpoases_config_initialize_default(solver_config);
#endif
            break;
        case DENSE_QP_QORE:
#ifdef ACADOS_WITH_QORE
            dense_qp_qore_config_initialize_default(solver_config);
#endif
            break;
    }
    return solver_config;
}

dense_qp_dims *dense_qp_dims_create()
{
    int bytes = dense_qp_dims_calculate_size();

    void *ptr = calloc(1, bytes);

    dense_qp_dims *dims = dense_qp_dims_assign(ptr);

    return dims;
}

dense_qp_in *dense_qp_in_create(qp_solver_config *config, dense_qp_dims *dims)
{
    int bytes = dense_qp_in_calculate_size(config, dims);

    void *ptr = calloc(1, bytes);

    dense_qp_in *in = dense_qp_in_assign(config, dims, ptr);

    return in;
}

dense_qp_out *dense_qp_out_create(qp_solver_config *config, dense_qp_dims *dims)
{
    int bytes = dense_qp_out_calculate_size(config, dims);

    void *ptr = calloc(1, bytes);

    dense_qp_out *out = dense_qp_out_assign(config, dims, ptr);

    return out;
}

void *dense_qp_opts_create(qp_solver_config *config, dense_qp_dims *dims)
{
    int bytes = config->opts_calculate_size(config, dims);

    void *ptr = calloc(1, bytes);

    void *opts = config->opts_assign(config, dims, ptr);

    config->opts_initialize_default(config, dims, opts);

    return opts;
}

int dense_qp_calculate_size(qp_solver_config *config, dense_qp_dims *dims, void *opts_)
{
    int bytes = sizeof(dense_qp_solver);

    bytes += config->memory_calculate_size(config, dims, opts_);
    bytes += config->workspace_calculate_size(config, dims, opts_);

    return bytes;
}

dense_qp_solver *dense_qp_assign(qp_solver_config *config, dense_qp_dims *dims, void *opts_,
                                 void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    dense_qp_solver *solver = (dense_qp_solver *) c_ptr;
    c_ptr += sizeof(dense_qp_solver);

    solver->config = config;
    solver->dims = dims;
    solver->opts = opts_;

    // TODO(dimitris): CHECK ALIGNMENT!

    solver->mem = config->memory_assign(config, dims, opts_, c_ptr);
    c_ptr += config->memory_calculate_size(config, dims, opts_);

    solver->work = (void *) c_ptr;
    c_ptr += config->workspace_calculate_size(config, dims, opts_);

    assert((char *) raw_memory + dense_qp_calculate_size(config, dims, opts_) == c_ptr);

    return solver;
}

dense_qp_solver *dense_qp_create(qp_solver_config *config, dense_qp_dims *dims, void *opts_)
{
    int bytes = dense_qp_calculate_size(config, dims, opts_);

    void *ptr = calloc(1, bytes);

    dense_qp_solver *solver = dense_qp_assign(config, dims, opts_, ptr);

    return solver;
}

int dense_qp_solve(dense_qp_solver *solver, dense_qp_in *qp_in, dense_qp_out *qp_out)
{
    return solver->config->evaluate(solver->config, qp_in, qp_out, solver->opts, solver->mem,
                                    solver->work);
}

static dense_qp_res *dense_qp_res_create(dense_qp_dims *dims)
{
    int size = dense_qp_res_calculate_size(dims);
    void *ptr = acados_malloc(size, 1);
    dense_qp_res *qp_res = dense_qp_res_assign(dims, ptr);
    return qp_res;
}

static dense_qp_res_ws *dense_qp_res_workspace_create(dense_qp_dims *dims)
{
    int size = dense_qp_res_workspace_calculate_size(dims);
    void *ptr = acados_malloc(size, 1);
    dense_qp_res_ws *res_ws = dense_qp_res_workspace_assign(dims, ptr);
    return res_ws;
}

// TODO(dimitris): better name for this wrapper?
void dense_qp_inf_norm_residuals(dense_qp_dims *dims, dense_qp_in *qp_in, dense_qp_out *qp_out,
                                 double *res)
{
    // double *residuals = malloc(4*sizeof(double));
    dense_qp_res *qp_res = dense_qp_res_create(dims);
    dense_qp_res_ws *res_ws = dense_qp_res_workspace_create(dims);
    dense_qp_res_compute(qp_in, qp_out, qp_res, res_ws);
    dense_qp_res_compute_nrm_inf(qp_res, res);
    free(qp_res);
    free(res_ws);
}
