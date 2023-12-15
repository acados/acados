/*
 * Copyright (c) The acados authors.
 *
 * This file is part of acados.
 *
 * The 2-Clause BSD License
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.;
 */


#include "acados_c/ocp_qp_interface.h"

// external
#include <assert.h>
#include <stdlib.h>
#include <string.h>
// acados_c

#include "acados/utils/mem.h"

#include "acados/dense_qp/dense_qp_hpipm.h"
#include "acados/ocp_qp/ocp_qp_xcond_solver.h"
#include "acados/ocp_qp/ocp_qp_partial_condensing.h"
#include "acados/ocp_qp/ocp_qp_full_condensing.h"

#ifdef ACADOS_WITH_QORE
#include "acados/dense_qp/dense_qp_qore.h"
#endif

#ifdef ACADOS_WITH_QPOASES
#include "acados/dense_qp/dense_qp_qpoases.h"
#endif

#ifdef ACADOS_WITH_DAQP
#include "acados/dense_qp/dense_qp_daqp.h"
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


// TODO: no "plan" is entering, rename?!
void ocp_qp_xcond_solver_config_initialize_from_plan(
        ocp_qp_solver_t solver_name, ocp_qp_xcond_solver_config *solver_config)
{

    switch (solver_name)
    {
        case PARTIAL_CONDENSING_HPIPM:
            ocp_qp_xcond_solver_config_initialize_default(solver_config);
            ocp_qp_hpipm_config_initialize_default(solver_config->qp_solver);
            ocp_qp_partial_condensing_config_initialize_default(solver_config->xcond);
            break;
#ifdef ACADOS_WITH_HPMPC
        case PARTIAL_CONDENSING_HPMPC:
            ocp_qp_xcond_solver_config_initialize_default(solver_config);
            ocp_qp_hpmpc_config_initialize_default(solver_config->qp_solver);
            ocp_qp_partial_condensing_config_initialize_default(solver_config->xcond);
            break;
#endif
#ifdef ACADOS_WITH_OOQP
        case PARTIAL_CONDENSING_OOQP:
            ocp_qp_xcond_solver_config_initialize_default(solver_config);
            ocp_qp_ooqp_config_initialize_default(solver_config->qp_solver);
            ocp_qp_partial_condensing_config_initialize_default(solver_config->xcond);
            break;
#endif
#ifdef ACADOS_WITH_OSQP
        case PARTIAL_CONDENSING_OSQP:
            ocp_qp_xcond_solver_config_initialize_default(solver_config);
            ocp_qp_osqp_config_initialize_default(solver_config->qp_solver);
            ocp_qp_partial_condensing_config_initialize_default(solver_config->xcond);
            break;
#endif
#ifdef ACADOS_WITH_QPDUNES
        case PARTIAL_CONDENSING_QPDUNES:
            ocp_qp_xcond_solver_config_initialize_default(solver_config);
            ocp_qp_qpdunes_config_initialize_default(solver_config->qp_solver);
            ocp_qp_partial_condensing_config_initialize_default(solver_config->xcond);
            break;
#endif
        case FULL_CONDENSING_HPIPM:
            ocp_qp_xcond_solver_config_initialize_default(solver_config);
            dense_qp_hpipm_config_initialize_default(solver_config->qp_solver);
            ocp_qp_full_condensing_config_initialize_default(solver_config->xcond);
            break;
#ifdef ACADOS_WITH_QPOASES
        case FULL_CONDENSING_QPOASES:
            ocp_qp_xcond_solver_config_initialize_default(solver_config);
            dense_qp_qpoases_config_initialize_default(solver_config->qp_solver);
            ocp_qp_full_condensing_config_initialize_default(solver_config->xcond);
            break;
#endif
#ifdef ACADOS_WITH_DAQP
        case FULL_CONDENSING_DAQP:
            ocp_qp_xcond_solver_config_initialize_default(solver_config);
            dense_qp_daqp_config_initialize_default(solver_config->qp_solver);
            ocp_qp_full_condensing_config_initialize_default(solver_config->xcond);
            break;
#endif
#ifdef ACADOS_WITH_QORE
        case FULL_CONDENSING_QORE:
            ocp_qp_xcond_solver_config_initialize_default(solver_config);
            dense_qp_qore_config_initialize_default(solver_config->qp_solver);
            ocp_qp_full_condensing_config_initialize_default(solver_config->xcond);
            break;
#endif
#ifdef ACADOS_WITH_OOQP
        case FULL_CONDENSING_OOQP:
            printf("\nocp_qp_xcond_solver: FULL_CONDENSING_OOQP.\n");
            ocp_qp_xcond_solver_config_initialize_default(solver_config);
            dense_qp_ooqp_config_initialize_default(solver_config->qp_solver);
            ocp_qp_full_condensing_config_initialize_default(solver_config->xcond);
            break;
#endif
        case INVALID_QP_SOLVER:
            printf("\nerror: ocp_qp_xcond_solver_config_initialize_from_plan: forgot to initialize plan->qp_solver\n");
            exit(1);
            break;
        default:
            printf("\nerror: ocp_qp_xcond_solver_config_initialize_from_plan: unsupported plan->qp_solver\n");
            printf("This might happen, if acados was not compiled with the specified QP solver.\n");
            exit(1);
    }

    return;
}



ocp_qp_xcond_solver_config *ocp_qp_xcond_solver_config_create(ocp_qp_solver_plan_t plan)
{
    acados_size_t bytes = ocp_qp_xcond_solver_config_calculate_size();
    void *ptr = calloc(1, bytes);
    ocp_qp_xcond_solver_config *solver_config = ocp_qp_xcond_solver_config_assign(ptr);

    ocp_qp_xcond_solver_config_initialize_from_plan(plan.qp_solver, solver_config);

    return solver_config;
}



void ocp_qp_xcond_solver_config_free(ocp_qp_xcond_solver_config *config)
{
    free(config);
}



/* dims */

ocp_qp_dims *ocp_qp_dims_create(int N)
{
    acados_size_t bytes = ocp_qp_dims_calculate_size(N);

    void *ptr = calloc(1, bytes);

    ocp_qp_dims *dims = ocp_qp_dims_assign(N, ptr);
    dims->N = N;

    return dims;
}



void ocp_qp_dims_free(void *dims_)
{
    free(dims_);
}


ocp_qp_xcond_solver_dims *ocp_qp_xcond_solver_dims_create(ocp_qp_xcond_solver_config *config, int N)
{
    acados_size_t bytes = ocp_qp_xcond_solver_dims_calculate_size(config, N);

    void *ptr = calloc(1, bytes);

    ocp_qp_xcond_solver_dims *dims = ocp_qp_xcond_solver_dims_assign(config, N, ptr);

    return dims;
}


ocp_qp_xcond_solver_dims *ocp_qp_xcond_solver_dims_create_from_ocp_qp_dims(
    ocp_qp_xcond_solver_config *config, ocp_qp_dims *dims)
{
    int N = dims->N;
    int tmp_int;

    ocp_qp_xcond_solver_dims *solver_dims = ocp_qp_xcond_solver_dims_create(config, N);

    /* alternative idea: call ocp_qp_dims_get and ocp_qp_xcond_dims_set for all fields; */
    for (int i = 0; i < N+1; i++)
    {
        ocp_qp_dims_get(config, dims, i, "nx", &tmp_int);
        ocp_qp_xcond_solver_dims_set_(config, solver_dims, i, "nx", &tmp_int);
        ocp_qp_dims_get(config, dims, i, "nu", &tmp_int);
        ocp_qp_xcond_solver_dims_set_(config, solver_dims, i, "nu", &tmp_int);
        ocp_qp_dims_get(config, dims, i, "nbx", &tmp_int);
        ocp_qp_xcond_solver_dims_set_(config, solver_dims, i, "nbx", &tmp_int);
        ocp_qp_dims_get(config, dims, i, "nbu", &tmp_int);
        ocp_qp_xcond_solver_dims_set_(config, solver_dims, i, "nbu", &tmp_int);
        ocp_qp_dims_get(config, dims, i, "ng", &tmp_int);
        ocp_qp_xcond_solver_dims_set_(config, solver_dims, i, "ng", &tmp_int);
        ocp_qp_dims_get(config, dims, i, "ns", &tmp_int);
        ocp_qp_xcond_solver_dims_set_(config, solver_dims, i, "ns", &tmp_int);

        ocp_qp_dims_get(config, dims, i, "nsbx", &tmp_int);
        ocp_qp_xcond_solver_dims_set_(config, solver_dims, i, "nsbx", &tmp_int);
        ocp_qp_dims_get(config, dims, i, "nsbu", &tmp_int);
        ocp_qp_xcond_solver_dims_set_(config, solver_dims, i, "nsbu", &tmp_int);
        ocp_qp_dims_get(config, dims, i, "nsg", &tmp_int);
        ocp_qp_xcond_solver_dims_set_(config, solver_dims, i, "nsg", &tmp_int);

        ocp_qp_dims_get(config, dims, i, "nbxe", &tmp_int);
        ocp_qp_xcond_solver_dims_set_(config, solver_dims, i, "nbxe", &tmp_int);
        ocp_qp_dims_get(config, dims, i, "nbue", &tmp_int);
        ocp_qp_xcond_solver_dims_set_(config, solver_dims, i, "nbue", &tmp_int);
        ocp_qp_dims_get(config, dims, i, "nge", &tmp_int);
        ocp_qp_xcond_solver_dims_set_(config, solver_dims, i, "nge", &tmp_int);
    }

    return solver_dims;
}


void ocp_qp_xcond_solver_dims_free(ocp_qp_xcond_solver_dims *dims)
{
    free(dims);
}



/* in */

ocp_qp_in *ocp_qp_in_create(ocp_qp_dims *dims)
{
    acados_size_t bytes = ocp_qp_in_calculate_size(dims);

    void *ptr = calloc(1, bytes);

    ocp_qp_in *in = ocp_qp_in_assign(dims, ptr);

    return in;
}


void ocp_qp_in_set(ocp_qp_xcond_solver_config *config, ocp_qp_in *in,
                   int stage, char *field, void *value)
{
    d_ocp_qp_set(field, stage, value, in);
}



void ocp_qp_in_free(void *in_)
{
    free(in_);
}



/* out */

ocp_qp_out *ocp_qp_out_create(ocp_qp_dims *dims)
{
    acados_size_t bytes = ocp_qp_out_calculate_size(dims);

    void *ptr = calloc(1, bytes);

    ocp_qp_out *out = ocp_qp_out_assign(dims, ptr);

    return out;
}



void ocp_qp_out_free(void *out_)
{
    free(out_);
}


void ocp_qp_out_get(ocp_qp_out *out, const char *field, void *value)
{
    if (!strcmp(field, "qp_info"))
    {
        qp_info **ptr = value;
        *ptr = out->misc;
    }
    else
    {
        printf("\nerror: ocp_qp_out_get: field %s not available\n", field);
        exit(1);
    }

    return;
}


/* opts */
void *ocp_qp_xcond_solver_opts_create(ocp_qp_xcond_solver_config *config, ocp_qp_xcond_solver_dims *dims)
{
    acados_size_t bytes = config->opts_calculate_size(config, dims);

    void *ptr = calloc(1, bytes);

    void *opts = config->opts_assign(config, dims, ptr);

    config->opts_initialize_default(config, dims, opts);

    return opts;
}



void ocp_qp_xcond_solver_opts_free(ocp_qp_xcond_solver_opts *opts)
{
    free(opts);
}


void ocp_qp_xcond_solver_opts_set(ocp_qp_xcond_solver_config *config,
           ocp_qp_xcond_solver_opts *opts, const char *field, void* value)
{
    ocp_qp_xcond_solver_opts_set_(config, opts, field, value);
}

/* solver */

acados_size_t ocp_qp_calculate_size(ocp_qp_xcond_solver_config *config, ocp_qp_xcond_solver_dims *dims, void *opts_)
{
    acados_size_t bytes = sizeof(ocp_qp_solver);

    bytes += config->memory_calculate_size(config, dims, opts_);
    bytes += config->workspace_calculate_size(config, dims, opts_);

    return bytes;
}



ocp_qp_solver *ocp_qp_assign(ocp_qp_xcond_solver_config *config, ocp_qp_xcond_solver_dims *dims,
                             void *opts_, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    ocp_qp_solver *solver = (ocp_qp_solver *) c_ptr;
    c_ptr += sizeof(ocp_qp_solver);

    solver->config = config;
    solver->dims = dims;
    solver->opts = opts_;

    solver->mem = config->memory_assign(config, dims, opts_, c_ptr);
    c_ptr += config->memory_calculate_size(config, dims, opts_);

    solver->work = (void *) c_ptr;
    c_ptr += config->workspace_calculate_size(config, dims, opts_);

    assert((char *) raw_memory + ocp_qp_calculate_size(config, dims, opts_) == c_ptr);

    return solver;
}



ocp_qp_solver *ocp_qp_create(ocp_qp_xcond_solver_config *config,
                             ocp_qp_xcond_solver_dims *dims, void *opts_)
{

    config->opts_update(config, dims, opts_);

    acados_size_t bytes = ocp_qp_calculate_size(config, dims, opts_);

    void *ptr = calloc(1, bytes);

    ocp_qp_solver *solver = ocp_qp_assign(config, dims, opts_, ptr);

    return solver;
}



int ocp_qp_solve(ocp_qp_solver *solver, ocp_qp_in *qp_in, ocp_qp_out *qp_out)
{
    return solver->config->evaluate(solver->config, solver->dims, qp_in, qp_out,
                                    solver->opts, solver->mem, solver->work);
}


void ocp_qp_solver_destroy(ocp_qp_solver *solver)
{
    free(solver);
}


// qp residual
static ocp_qp_res *ocp_qp_res_create(ocp_qp_dims *dims)
{
    acados_size_t size = ocp_qp_res_calculate_size(dims);
    void *ptr = acados_malloc(size, 1);
    assert(ptr != 0);
    ocp_qp_res *qp_res = ocp_qp_res_assign(dims, ptr);
    return qp_res;
}



static ocp_qp_res_ws *ocp_qp_res_workspace_create(ocp_qp_dims *dims)
{
    acados_size_t size = ocp_qp_res_workspace_calculate_size(dims);
    void *ptr = acados_malloc(size, 1);
    assert(ptr != 0);
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
