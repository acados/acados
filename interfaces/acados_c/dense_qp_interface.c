/*
 * Copyright 2019 Gianluca Frison, Dimitris Kouzoupis, Robin Verschueren,
 * Andrea Zanelli, Niels van Duijkeren, Jonathan Frey, Tommaso Sartor,
 * Branimir Novoselnik, Rien Quirynen, Rezart Qelibari, Dang Doan,
 * Jonas Koenemann, Yutao Chen, Tobias Schöls, Jonas Schlagenhauf, Moritz Diehl
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


#include "acados_c/dense_qp_interface.h"

// external
#include <assert.h>
#include <stdlib.h>
#include <string.h>

// hpipm
#include "hpipm/include/hpipm_d_dense_qp.h"

// acados_c

#include "acados/utils/mem.h"

#include "acados/dense_qp/dense_qp_hpipm.h"
#ifdef ACADOS_WITH_QORE
#include "acados/dense_qp/dense_qp_qore.h"
#endif
#ifdef ACADOS_WITH_QPOASES
#include "acados/dense_qp/dense_qp_qpoases.h"
#endif
#ifdef ACADOS_WITH_OOQP
#include "acados/dense_qp/dense_qp_ooqp.h"
#endif

qp_solver_config *dense_qp_config_create(dense_qp_solver_plan *plan)
{
    acados_size_t bytes = dense_qp_solver_config_calculate_size();
    void *ptr = calloc(1, bytes);
    qp_solver_config *solver_config = dense_qp_solver_config_assign(ptr);

    dense_qp_solver_t solver_name = plan->qp_solver;

    switch (solver_name)
    {
        case DENSE_QP_HPIPM:
            dense_qp_hpipm_config_initialize_default(solver_config);
            break;
#ifdef ACADOS_WITH_QPOASES
        case DENSE_QP_QPOASES:
            dense_qp_qpoases_config_initialize_default(solver_config);
            break;
#endif
#ifdef ACADOS_WITH_QORE
        case DENSE_QP_QORE:
            dense_qp_qore_config_initialize_default(solver_config);
            break;
#endif
#ifdef ACADOS_WITH_OOQP
        case DENSE_QP_OOQP:
            dense_qp_ooqp_config_initialize_default(solver_config);
            break;
#endif
        default:
            printf("\nerror: dense_qp_config_create: unsupported plan->qp_solver.\n");
            printf("This might happen, if acados was not compiled with the specified QP solver.\n");
            exit(1);
            break;
    }
    return solver_config;
}

dense_qp_dims *dense_qp_dims_create()
{
    acados_size_t bytes = dense_qp_dims_calculate_size();

    void *ptr = calloc(1, bytes);

    dense_qp_dims *dims = dense_qp_dims_assign(ptr);

    return dims;
}

dense_qp_in *dense_qp_in_create(qp_solver_config *config, dense_qp_dims *dims)
{
    acados_size_t bytes = dense_qp_in_calculate_size(dims);

    void *ptr = calloc(1, bytes);

    dense_qp_in *in = dense_qp_in_assign(dims, ptr);

    return in;
}

dense_qp_out *dense_qp_out_create(qp_solver_config *config, dense_qp_dims *dims)
{
    acados_size_t bytes = dense_qp_out_calculate_size(dims);

    void *ptr = calloc(1, bytes);

    dense_qp_out *out = dense_qp_out_assign(dims, ptr);

    return out;
}

void *dense_qp_opts_create(qp_solver_config *config, dense_qp_dims *dims)
{
    acados_size_t bytes = config->opts_calculate_size(config, dims);

    void *ptr = calloc(1, bytes);

    void *opts = config->opts_assign(config, dims, ptr);

    config->opts_initialize_default(config, dims, opts);

    return opts;
}

acados_size_t dense_qp_calculate_size(qp_solver_config *config, dense_qp_dims *dims, void *opts_)
{
    acados_size_t bytes = sizeof(dense_qp_solver);

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
    acados_size_t bytes = dense_qp_calculate_size(config, dims, opts_);

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
    acados_size_t size = dense_qp_res_calculate_size(dims);
    void *ptr = acados_malloc(size, 1);
    assert(ptr != 0);
    dense_qp_res *qp_res = dense_qp_res_assign(dims, ptr);
    return qp_res;
}

static dense_qp_res_ws *dense_qp_res_workspace_create(dense_qp_dims *dims)
{
    acados_size_t size = dense_qp_res_workspace_calculate_size(dims);
    void *ptr = acados_malloc(size, 1);
    assert(ptr != 0);
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

bool dense_qp_set_field_double_array(const char *field, double *arr, dense_qp_in *qp_in)
{
    if (!strcmp(field, "H"))
    {
        d_dense_qp_set_H(arr, qp_in);
    }
    else if (!strcmp(field, "g"))
    {
        d_dense_qp_set_g(arr, qp_in);
    }
    else if (!strcmp(field, "A"))
    {
        d_dense_qp_set_A(arr, qp_in);
    }
    else if (!strcmp(field, "b"))
    {
        d_dense_qp_set_b(arr, qp_in);
    }
    else if (!strcmp(field, "lb"))
    {
        d_dense_qp_set_lb(arr, qp_in);
    }
    else if (!strcmp(field, "ub"))
    {
        d_dense_qp_set_ub(arr, qp_in);
    }
    else if (!strcmp(field, "C"))
    {
        d_dense_qp_set_C(arr, qp_in);
    }
    else if (!strcmp(field, "lg"))
    {
        d_dense_qp_set_lg(arr, qp_in);
    }
    else if (!strcmp(field, "ug"))
    {
        d_dense_qp_set_ug(arr, qp_in);
    }
    else if (!strcmp(field, "Zl"))
    {
        d_dense_qp_set_Zl(arr, qp_in);
    }
    else if (!strcmp(field, "Zu"))
    {
        d_dense_qp_set_Zu(arr, qp_in);
    }
    else if (!strcmp(field, "zl"))
    {
        d_dense_qp_set_zl(arr, qp_in);
    }
    else if (!strcmp(field, "zu"))
    {
        d_dense_qp_set_zu(arr, qp_in);
    }
    else if (!strcmp(field, "ls"))
    {
        d_dense_qp_set_ls(arr, qp_in);
    }
    else if (!strcmp(field, "us"))
    {
        d_dense_qp_set_us(arr, qp_in);
    }
    else
    {
        printf("\n%s is an unknown double array field in dense_qp_in!\n", field);
        return false;
    }

    return true;
}

bool dense_qp_set_field_int_array(const char *field, int *arr, dense_qp_in *qp_in)
{
    if (!strcmp(field, "idxb"))
    {
        d_dense_qp_set_idxb(arr, qp_in);
    }
    else if (!strcmp(field, "idxs"))
    {
        d_dense_qp_set_idxs(arr, qp_in);
    }
    else
    {
        printf("\n%s is an unknown int array field in dense_qp_in!\n", field);
        return false;
    }

    return true;
}

bool dense_qp_get_field_double_array(const char *field, dense_qp_in *qp_in, double *arr)
{
    if (!strcmp(field, "H"))
    {
        d_dense_qp_get_H(qp_in, arr);
    }
    else if (!strcmp(field, "g"))
    {
        d_dense_qp_get_g(qp_in, arr);
    }
    else if (!strcmp(field, "A"))
    {
        d_dense_qp_get_A(qp_in, arr);
    }
    else if (!strcmp(field, "b"))
    {
        d_dense_qp_get_b(qp_in, arr);
    }
    else if (!strcmp(field, "lb"))
    {
        d_dense_qp_get_lb(qp_in, arr);
    }
    else if (!strcmp(field, "ub"))
    {
        d_dense_qp_get_ub(qp_in, arr);
    }
    else if (!strcmp(field, "C"))
    {
        d_dense_qp_get_C(qp_in, arr);
    }
    else if (!strcmp(field, "lg"))
    {
        d_dense_qp_get_lg(qp_in, arr);
    }
    else if (!strcmp(field, "ug"))
    {
        d_dense_qp_get_ug(qp_in, arr);
    }
    else if (!strcmp(field, "Zl"))
    {
        d_dense_qp_get_Zl(qp_in, arr);
    }
    else if (!strcmp(field, "Zu"))
    {
        d_dense_qp_get_Zu(qp_in, arr);
    }
    else if (!strcmp(field, "zl"))
    {
        d_dense_qp_get_zl(qp_in, arr);
    }
    else if (!strcmp(field, "zu"))
    {
        d_dense_qp_get_zu(qp_in, arr);
    }
    else if (!strcmp(field, "ls"))
    {
        d_dense_qp_get_ls(qp_in, arr);
    }
    else if (!strcmp(field, "us"))
    {
        d_dense_qp_get_us(qp_in, arr);
    }
    else
    {
        printf("\n%s is an unknown double array field in dense_qp_in!\n", field);
        return false;
    }

    return true;
}
bool dense_qp_get_field_int_array(const char *field, dense_qp_in *qp_in, int *arr)
{
    if (!strcmp(field, "idxb"))
    {
        d_dense_qp_get_idxb(qp_in, arr);
    }
    else if (!strcmp(field, "idxs"))
    {
        d_dense_qp_get_idxs(qp_in, arr);
    }
    else
    {
        printf("\n%s is an unknown int array field in dense_qp_in!\n", field);
        return false;
    }

    return true;
}
