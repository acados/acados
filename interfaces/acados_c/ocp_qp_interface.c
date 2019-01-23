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

#include "acados_c/ocp_qp_interface.h"

// external
#include <assert.h>
#include <stdlib.h>
#include <string.h>
// acados_c

#include "acados/utils/mem.h"

#include "acados/dense_qp/dense_qp_hpipm.h"
#include "acados/ocp_qp/ocp_qp_full_condensing_solver.h"
#include "acados/ocp_qp/ocp_qp_partial_condensing_solver.h"
#ifdef ACADOS_WITH_QORE
#include "acados/dense_qp/dense_qp_qore.h"
#endif
#ifdef ACADOS_WITH_QPOASES
#include "acados/dense_qp/dense_qp_qpoases.h"
#endif
#include "acados/ocp_qp/ocp_qp_hpipm.h"
#ifdef ACADOS_WITH_HPMPC
#include "acados/ocp_qp/ocp_qp_hpmpc.h"
#endif
#ifdef ACADOS_WITH_QPDUNES
#include "acados/ocp_qp/ocp_qp_qpdunes.h"
#endif
#ifdef ACADOS_WITH_OOQP
#include "acados/dense_qp/dense_qp_ooqp.h"
#include "acados/ocp_qp/ocp_qp_ooqp.h"
#endif
#ifdef ACADOS_WITH_OSQP
#include "acados/ocp_qp/ocp_qp_osqp.h"
#endif

void ocp_qp_xcond_solver_config_initialize_default(ocp_qp_solver_t solver_name,
                                                   ocp_qp_xcond_solver_config *solver_config)
{
// NOTE: this only works if solvers are ordered in the enum !!!!!!!!!!!!!!!!
if (solver_name < FULL_CONDENSING_HPIPM)
    {
        ocp_qp_partial_condensing_solver_config_initialize_default(solver_config);
    }
    else
    {
        ocp_qp_full_condensing_solver_config_initialize_default(solver_config);
    }

    switch (solver_name)
    {
        case PARTIAL_CONDENSING_HPIPM:
            ocp_qp_hpipm_config_initialize_default(solver_config->qp_solver);
            break;
#ifdef ACADOS_WITH_HPMPC
        case PARTIAL_CONDENSING_HPMPC:
            ocp_qp_hpmpc_config_initialize_default(solver_config->qp_solver);
            break;
#endif
#ifdef ACADOS_WITH_OOQP
        case PARTIAL_CONDENSING_OOQP:
            ocp_qp_ooqp_config_initialize_default(solver_config->qp_solver);
            break;
#endif
#ifdef ACADOS_WITH_OSQP
        case PARTIAL_CONDENSING_OSQP:
            ocp_qp_osqp_config_initialize_default(solver_config->qp_solver);
            break;
#endif
#ifdef ACADOS_WITH_QPDUNES
        case PARTIAL_CONDENSING_QPDUNES:
            ocp_qp_qpdunes_config_initialize_default(solver_config->qp_solver);
            break;
#endif
        case FULL_CONDENSING_HPIPM:
            dense_qp_hpipm_config_initialize_default(solver_config->qp_solver);
            break;
#ifdef ACADOS_WITH_QPOASES
        case FULL_CONDENSING_QPOASES:
            dense_qp_qpoases_config_initialize_default(solver_config->qp_solver);
            break;
#endif
#ifdef ACADOS_WITH_QORE
        case FULL_CONDENSING_QORE:
            dense_qp_qore_config_initialize_default(solver_config->qp_solver);
            break;
#endif
#ifdef ACADOS_WITH_OOQP
        case FULL_CONDENSING_OOQP:
            dense_qp_ooqp_config_initialize_default(solver_config->qp_solver);
            break;
#endif
        default:
            printf("\nerror: ocp_qp_config_create: unsupported plan->qp_solver\n");
            exit(1);
            break;
    }
}


ocp_qp_xcond_solver_config *ocp_qp_config_create(ocp_qp_solver_plan plan)
{
    int bytes = ocp_qp_xcond_solver_config_calculate_size();
    void *ptr = calloc(1, bytes);
    ocp_qp_xcond_solver_config *solver_config = ocp_qp_xcond_solver_config_assign(ptr);

    ocp_qp_xcond_solver_config_initialize_default(plan.qp_solver, solver_config);

    return solver_config;
}


void ocp_qp_config_free(void *config)
{
    free(config);
}


ocp_qp_dims *ocp_qp_dims_create(int N)
{
    int bytes = ocp_qp_dims_calculate_size(N);

    void *ptr = calloc(1, bytes);

    ocp_qp_dims *dims = ocp_qp_dims_assign(N, ptr);
    dims->N = N;

    return dims;
}

void ocp_qp_dims_free(void *dims_)
{
    free(dims_);
}


/* in */
ocp_qp_in *ocp_qp_in_create(ocp_qp_xcond_solver_config *config, ocp_qp_dims *dims)
{
    int bytes = ocp_qp_in_calculate_size(config, dims);

    void *ptr = calloc(1, bytes);

    ocp_qp_in *in = ocp_qp_in_assign(config, dims, ptr);

    return in;
}

void ocp_qp_in_free(void *in_)
{
    free(in_);
}


/* out */
ocp_qp_out *ocp_qp_out_create(ocp_qp_xcond_solver_config *config, ocp_qp_dims *dims)
{
    int bytes = ocp_qp_out_calculate_size(config, dims);

    void *ptr = calloc(1, bytes);

    ocp_qp_out *out = ocp_qp_out_assign(config, dims, ptr);

    return out;
}


void ocp_qp_out_free(void *out_)
{
    free(out_);
}

/* opts */
void *ocp_qp_opts_create(ocp_qp_xcond_solver_config *config, ocp_qp_dims *dims)
{
    int bytes = config->opts_calculate_size(config, dims);

    void *ptr = calloc(1, bytes);

    void *opts = config->opts_assign(config, dims, ptr);

    config->opts_initialize_default(config, dims, opts);

    return opts;
}


void ocp_qp_opts_free(void *opts_)
{
    free(opts_);
}


/* solver */

int ocp_qp_calculate_size(ocp_qp_xcond_solver_config *config, ocp_qp_dims *dims, void *opts_)
{
    int bytes = sizeof(ocp_qp_solver);

    bytes += config->memory_calculate_size(config, dims, opts_);
    bytes += config->workspace_calculate_size(config, dims, opts_);

    return bytes;
}

ocp_qp_solver *ocp_qp_assign(ocp_qp_xcond_solver_config *config, ocp_qp_dims *dims, void *opts_,
                             void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    ocp_qp_solver *solver = (ocp_qp_solver *) c_ptr;
    c_ptr += sizeof(ocp_qp_solver);

    solver->config = config;
    solver->dims = dims;
    solver->opts = opts_;

    // TODO(dimitris): CHECK ALIGNMENT!

    solver->mem = config->memory_assign(config, dims, opts_, c_ptr);
    c_ptr += config->memory_calculate_size(config, dims, opts_);

    solver->work = (void *) c_ptr;
    c_ptr += config->workspace_calculate_size(config, dims, opts_);

    assert((char *) raw_memory + ocp_qp_calculate_size(config, dims, opts_) == c_ptr);

    return solver;
}

ocp_qp_solver *ocp_qp_create(ocp_qp_xcond_solver_config *config, ocp_qp_dims *dims, void *opts_)
{
    config->opts_update(config, dims, opts_);

    int bytes = ocp_qp_calculate_size(config, dims, opts_);

    void *ptr = calloc(1, bytes);

    ocp_qp_solver *solver = ocp_qp_assign(config, dims, opts_, ptr);

    return solver;
}

int ocp_qp_solve(ocp_qp_solver *solver, ocp_qp_in *qp_in, ocp_qp_out *qp_out)
{
    return solver->config->evaluate(solver->config, qp_in, qp_out, solver->opts, solver->mem,
                                    solver->work);
}

static ocp_qp_res *ocp_qp_res_create(ocp_qp_dims *dims)
{
    int size = ocp_qp_res_calculate_size(dims);
    void *ptr = acados_malloc(size, 1);
    ocp_qp_res *qp_res = ocp_qp_res_assign(dims, ptr);
    return qp_res;
}

static ocp_qp_res_ws *ocp_qp_res_workspace_create(ocp_qp_dims *dims)
{
    int size = ocp_qp_res_workspace_calculate_size(dims);
    void *ptr = acados_malloc(size, 1);
    ocp_qp_res_ws *res_ws = ocp_qp_res_workspace_assign(dims, ptr);
    return res_ws;
}

// TODO(dimitris): better name for this wrapper?
void ocp_qp_inf_norm_residuals(ocp_qp_dims *dims, ocp_qp_in *qp_in, ocp_qp_out *qp_out, double *res)
{
    ocp_qp_res *qp_res = ocp_qp_res_create(dims);
    ocp_qp_res_ws *res_ws = ocp_qp_res_workspace_create(dims);
    ocp_qp_res_compute(qp_in, qp_out, qp_res, res_ws);
    ocp_qp_res_compute_nrm_inf(qp_res, res);
    free(qp_res);
    free(res_ws);
}
