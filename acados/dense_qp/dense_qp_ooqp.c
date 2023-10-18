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




#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "ooqp/cQpGenDense.h"

// blasfeo
#include "blasfeo/include/blasfeo_d_aux.h"
#include "blasfeo/include/blasfeo_d_aux_ext_dep.h"

// acados
#include "acados/dense_qp/dense_qp_ooqp.h"
#include "acados/utils/timing.h"
#include "acados/utils/print.h"
#include "acados/utils/mem.h"



static void update_gradient(const dense_qp_in *in, dense_qp_ooqp_memory *mem)
{
    dense_qp_dims *dims = in->dim;

    blasfeo_unpack_dvec(dims->nv, in->gz, 0, &mem->c[0]);
}



static void update_hessian_data(const dense_qp_in *in, dense_qp_ooqp_memory *mem)
{
    dense_qp_dims *dims = in->dim;

    // fill in the upper triangular of H in dense_qp
    blasfeo_dtrtr_l(dims->nv, in->Hv, 0, 0, in->Hv, 0, 0);

    // copy H to dQ, transpose to get row-major ordering
    blasfeo_unpack_tran_dmat(dims->nv, dims->nv, in->Hv, 0, 0, mem->dQ, dims->nv);
}



static void update_b_vector(const dense_qp_in *in, dense_qp_ooqp_memory *mem)
{
    dense_qp_dims *dims = in->dim;

    blasfeo_unpack_dvec(dims->ne, in->b, 0, mem->bA);
}



static void update_bounds(const dense_qp_in *in, dense_qp_ooqp_memory *mem)
{
    dense_qp_dims *dims = in->dim;

    for (int ii = 0; ii < dims->nv; ii++)
    {
        mem->ixlow[ii] = (char)0;
        mem->ixupp[ii] = (char)0;
        mem->xlow[ii] = 0.0;
        mem->xupp[ii] = 0.0;
    }

    for (int ii = 0; ii < dims->nb; ii++)
    {
        int idx = in->idxb[ii];
        // TODO(dimitris): check if cast is redundant
        // NOTE(dimitris): OOQP can give wrong results if there are 1e12 bounds
        if (BLASFEO_DVECEL(in->d, ii) > -1e10)
        {  // TODO(dimitris): use acados inf
            mem->ixlow[idx] = (char)1;
            mem->xlow[idx] = BLASFEO_DVECEL(in->d, ii);
        }
        if (-BLASFEO_DVECEL(in->d, ii+dims->nb+dims->ng) < 1e10)
        {  // TODO(dimitris): same here
            mem->ixupp[idx] = (char)1;
            mem->xupp[idx] = -BLASFEO_DVECEL(in->d, ii+dims->nb+dims->ng);
        }
    }
}



static void update_ineq_bounds(const dense_qp_in *in, dense_qp_ooqp_memory *mem)
{
    dense_qp_dims *dims = in->dim;

    for (int ii = 0; ii < dims->ng; ii++)
    {
        mem->iclow[ii] = (char)1;
        mem->clow[ii] = BLASFEO_DVECEL(in->d, ii+dims->nb);  // lg
        mem->icupp[ii] = (char)1;
        mem->cupp[ii] = -BLASFEO_DVECEL(in->d, ii+2*dims->nb+dims->ng);  // ug
    }
}



static void update_inequalities_data(const dense_qp_in *in, dense_qp_ooqp_memory *mem)
{
    dense_qp_dims *dims = in->dim;

    blasfeo_unpack_dmat(dims->nv, dims->ng, in->Ct, 0, 0, mem->dC, dims->nv);
}



static void dense_qp_ooqp_update_memory(const dense_qp_in *in, const dense_qp_ooqp_opts *opts,
                                        dense_qp_ooqp_memory *mem)
{
    // ------- Update objective
    update_gradient(in, mem);

    if (mem->firstRun == 1 || opts->fixHessian == 0)
    {
        update_hessian_data(in, mem);
    }

    // ------- Update equality constraints
    update_b_vector(in, mem);

    // ------- Update bounds
    update_bounds(in, mem);

    // ------- Update inequality constraints
    update_ineq_bounds(in, mem);

    if (mem->firstRun == 1 || opts->fixInequalities == 0)
    {
        update_inequalities_data(in, mem);
    }

    mem->firstRun = 0;
}



static void print_inputs(dense_qp_ooqp_memory *mem)
{
    printf("\n----------> OOQP INPUTS <----------\n\n");
    printf("NUMBER OF PRIMAL VARIABLES: %d\n", mem->nx);
    printf("NUMBER OF EQUALITY CONSTRAINTS: %d\n", mem->my);
    printf("NUMBER OF INEQUALITY CONSTRAINTS: %d\n", mem->mz);
    printf("\n-----------------------------------\n\n");

    printf("\nOBJECTIVE FUNCTION:\n");
    for (int jj = 0; jj < mem->nx; jj++)
    {
        for (int ii = 0; ii < mem->nx; ii++)
        {
            printf("%9.5f ", mem->dQ[ii+mem->nx*jj]);
        }
        printf("\n");
    }
    printf("\n");

    for (int ii = 0; ii < mem->nx; ii++)
        printf("===> c[%d] = %f\n", ii + 1, mem->c[ii]);

    printf("\nBOUNDS:\n");
    for (int ii = 0; ii < mem->nx; ii++)
    {
        printf(
            "ixlow[%d] = %d \t xlow[%d] = %4.2f \t ixupp[%d] = %d \t xupp[%d] "
            "= %4.2f\n",
            ii + 1, mem->ixlow[ii], ii + 1, mem->xlow[ii], ii + 1, mem->ixupp[ii], ii + 1,
            mem->xupp[ii]);
    }

    printf("\n");
    printf("\nINEQUALITY CONSTRAINTS:\n");
    for (int ii = 0; ii < mem->mz; ii++)
    {
        for (int jj = 0; jj < mem->nx; jj++)
        {
            printf("%9.5f ", mem->dC[ii*mem->nx + jj]);
        }
        printf("\n");
    }

    for (int ii = 0; ii < mem->mz; ii++)
    {
        printf("===> clow[%d] = %4.2f \t cupp[%d] = %4.2f\n", ii + 1, mem->clow[ii], ii + 1,
               mem->cupp[ii]);
    }
}



static void print_outputs(dense_qp_ooqp_memory *mem, dense_qp_ooqp_workspace *work,
                          int return_value)
{
    int_t ii;

    printf("\n----------> OOQP OUTPUTS <---------\n\n");
    printf("RETURN STATUS: %d\n", return_value);
    printf("OBJECTIVE VALUE: %f\n", work->objectiveValue);
    printf("FIRST AND LAST ELEMENT OF SOLUTION:\n");
    printf("x[0] = %f\n", work->x[0]);
    printf("x[%d] = %f\n", mem->nx, work->x[mem->nx - 1]);
    printf("\n----------------------------------\n\n");

    printf("\nPRIMAL SOLUTION:\n");
    for (ii = 0; ii < mem->nx; ii++)
    {
        printf("=====> x[%d] = %f\n", ii + 1, work->x[ii]);
    }

    // printf("\nGAMMA:\n");
    // for (ii = 0; ii < mem->nx; ii++)
    // {
    //     printf("=====> gamma[%d] = %f\n", ii + 1, work->gamma[ii]);
    // }

    // printf("\nPHI:\n");
    // for (ii = 0; ii < mem->nx; ii++)
    // {
    //     printf("=====> phi[%d] = %f\n", ii + 1, work->phi[ii]);
    // }

    printf("\nLAMBDA:\n");
    for (ii = 0; ii < mem->mz; ii++)
    {
        printf("=====> lambda[%d] = %f\n", ii + 1, work->lambda[ii]);
    }

    printf("\nPI:\n");
    for (ii = 0; ii < mem->mz; ii++)
    {
        printf("=====> pi[%d] = %f\n", ii + 1, work->pi[ii]);
    }
}



static void fill_in_qp_out(const dense_qp_in *in, dense_qp_out *out, dense_qp_ooqp_workspace *work)
{
    int nb = in->dim->nb;
    int ng = in->dim->ng;
    int ns = in->dim->ns;
    int nv = in->dim->nv;
    int *idxb = in->idxb;

    blasfeo_pack_dvec(nv, &work->x[0], out->v, 0);

    for (int ii = 0; ii < 2 * nb + 2 * ng + 2 * ns; ii++) out->lam->pa[ii] = 0.0;

    for (int ii = 0; ii < nb; ii++)
    {
        double delta = work->gamma[idxb[ii]] - work->phi[idxb[ii]];

        if (delta >= 0)
            out->lam->pa[ii] = delta;
        else
            out->lam->pa[nb + ng + ii] = -delta;
    }


    for (int ii = 0; ii < ng; ii++)
    {
        double delta = work->lambda[ii] - work->pi[ii];

        if (delta >= 0)
            out->lam->pa[nb+ii] = delta;
        else
            out->lam->pa[2*nb+ng+ii] = -delta;
    }

}



acados_size_t dense_qp_ooqp_opts_calculate_size(void *config_, dense_qp_dims *dims)
{
    acados_size_t size = 0;
    size += sizeof(dense_qp_ooqp_opts);
    return size;
}



void *dense_qp_ooqp_opts_assign(void *config_, dense_qp_dims *dims, void *raw_memory)
{
    dense_qp_ooqp_opts *opts;

    char *c_ptr = (char *) raw_memory;

    opts = (dense_qp_ooqp_opts *) c_ptr;
    c_ptr += sizeof(dense_qp_ooqp_opts);

    assert((char *) raw_memory + dense_qp_ooqp_opts_calculate_size(config_, dims) >= c_ptr);

    return (void *) opts;
}



void dense_qp_ooqp_opts_initialize_default(void *config_, dense_qp_dims *dims, void *opts_)
{
    dense_qp_ooqp_opts *opts = (dense_qp_ooqp_opts *) opts_;

    opts->printLevel = 0;
    opts->fixHessian = 0;
    opts->fixDynamics = 0;
    opts->fixInequalities = 0;

    return;
}



void dense_qp_ooqp_opts_update(void *config_, dense_qp_dims *dims, void *opts_)
{
    return;
}



void dense_qp_ooqp_opts_set(void *config_, void *opts_, const char *field, void *value)
{
    dense_qp_ooqp_opts *opts = opts_;

    if (!strcmp(field, "tol_stat"))
    {
        // TODO set solver exit tolerance
    }
    else if (!strcmp(field, "tol_eq"))
    {
        // TODO set solver exit tolerance
    }
    else if (!strcmp(field, "tol_ineq"))
    {
        // TODO set solver exit tolerance
    }
    else if (!strcmp(field, "tol_comp"))
    {
        // TODO set solver exit tolerance
    }
    else if (!strcmp(field, "warm_start"))
    {
        // TODO set solver warm start
    }
    else
    {
        printf("\nerror: dense_qp_ooqp_opts_set: wrong field: %s\n", field);
        exit(1);
    }

    return;
}



acados_size_t dense_qp_ooqp_memory_calculate_size(void *config_, dense_qp_dims *dims, void *opts_)
{
    int nv = dims->nv;
    int ne = dims->ne;
    int ng = dims->ng;
    // int nb = dims->nb;
    // int ns = dims->ns;
    // int nsb = dims->nsb;
    // int nsg = dims->nsg;

    acados_size_t size = 0;
    size += sizeof(dense_qp_ooqp_memory);

    size += 1 * nv * nv * sizeof(double);  // dQ
    size += 1 * nv *      sizeof(double);  // c
    size += 2 * nv *      sizeof(char);    // ixlow, ixupp
    size += 2 * nv *      sizeof(double);  // xlow, xupp
    size += 1 * nv * ne * sizeof(double);  // dA
    size += 1 * ne *      sizeof(double);  // bA
    size += 1 * nv * ng * sizeof(double);  // dC
    size += 2 * ng *      sizeof(double);  // clow, cupp
    size += 2 * ng *      sizeof(char);    // iclow, icupp

    make_int_multiple_of(8, &size);

    return size;
}



void *dense_qp_ooqp_memory_assign(void *config_, dense_qp_dims *dims, void *opts_, void *raw_memory)
{
    // dense_qp_ooqp_opts *opts = (dense_qp_ooqp_opts *)opts_;
    dense_qp_ooqp_memory *mem;

    int nv = dims->nv;
    int ne = dims->ne;
    int ng = dims->ng;
    // int nb = dims->nb;
    // int ns = dims->ns;
    // int nsb = dims->nsb;
    // int nsg = dims->nsg;

    // char pointer
    char *c_ptr = (char *) raw_memory;

    mem = (dense_qp_ooqp_memory *) c_ptr;
    c_ptr += sizeof(dense_qp_ooqp_memory);

    assert((size_t) c_ptr % 8 == 0 && "memory not 8-byte aligned!");

    assign_and_advance_double(nv * nv, &mem->dQ, &c_ptr);
    assign_and_advance_double(nv, &mem->c, &c_ptr);
    assign_and_advance_double(nv, &mem->xlow, &c_ptr);
    assign_and_advance_double(nv, &mem->xupp, &c_ptr);
    assign_and_advance_double(nv * ne, &mem->dA, &c_ptr);
    assign_and_advance_double(ne, &mem->bA, &c_ptr);
    assign_and_advance_double(nv * ng, &mem->dC, &c_ptr);
    assign_and_advance_double(ng, &mem->clow, &c_ptr);
    assign_and_advance_double(ng, &mem->cupp, &c_ptr);

    assert((size_t) c_ptr % 8 == 0 && "double not 8-byte aligned!");

    assign_and_advance_char(nv, &mem->ixlow, &c_ptr);
    assign_and_advance_char(nv, &mem->ixupp, &c_ptr);
    assign_and_advance_char(ng, &mem->iclow, &c_ptr);
    assign_and_advance_char(ng, &mem->icupp, &c_ptr);

    // initialize memory
    mem->firstRun = 1;
    mem->nx = nv;
    mem->my = ne;
    mem->mz = ng;

    assert((char *) raw_memory + dense_qp_ooqp_memory_calculate_size(config_, dims, opts_) >=
           c_ptr);

    return mem;
}



void dense_qp_ooqp_memory_get(void *config_, void *mem_, const char *field, void* value)
{
    qp_solver_config *config = config_;
    dense_qp_ooqp_memory *mem = mem_;

    if(!strcmp(field, "time_qp_solver_call"))
    {
        double *tmp_ptr = value;
        *tmp_ptr = mem->time_qp_solver_call;
    }
    else if(!strcmp(field, "iter"))
    {
        int *tmp_ptr = value;
        *tmp_ptr = mem->iter;
    }
    else
    {
        printf("\nerror: dense_qp_ooqp_memory_get: field %s not available\n", field);
        exit(1);
    }

    return;

}



/************************************************
 * workspace
 ************************************************/

acados_size_t dense_qp_ooqp_workspace_calculate_size(void *config_, dense_qp_dims *dims, void *opts_)
{
    // dense_qp_ooqp_opts *opts = (dense_qp_ooqp_opts *)opts_;

    acados_size_t size = 0;
    acados_size_t nx, my, mz;

    nx = dims->nv;
    my = dims->ne;
    mz = dims->ng;

    size += sizeof(dense_qp_ooqp_workspace);
    size += sizeof(double) * (3 * nx + my + 3 * mz);

    return size;
}



static void dense_qp_ooqp_cast_workspace(dense_qp_ooqp_workspace *work, dense_qp_ooqp_memory *mem)
{
    char *ptr = (char *)work;

    ptr += sizeof(dense_qp_ooqp_workspace);
    work->x = (double *)ptr;
    ptr += (mem->nx) * sizeof(double);
    work->gamma = (double *)ptr;
    ptr += (mem->nx) * sizeof(double);
    work->phi = (double *)ptr;
    ptr += (mem->nx) * sizeof(double);
    work->y = (double *)ptr;
    ptr += (mem->my) * sizeof(double);
    work->z = (double *)ptr;
    ptr += (mem->mz) * sizeof(double);
    work->lambda = (double *)ptr;
    ptr += (mem->mz) * sizeof(double);
    work->pi = (double *)ptr;
    ptr += (mem->mz) * sizeof(double);
}



int_t dense_qp_ooqp(void *config_, dense_qp_in *qp_in, dense_qp_out *qp_out, void *opts_,
                    void *memory_, void *work_)
{
    int ns = qp_in->dim->ns;

    if (ns > 0)
    {
        printf("\nOOQP interface can not handle ns>0 yet: what about implementing it? :)\n");
        exit(1);
    }

    acados_timer tot_timer, qp_timer, interface_timer;
    qp_info *info = (qp_info *) qp_out->misc;
    acados_tic(&tot_timer);

    dense_qp_ooqp_opts *opts = (dense_qp_ooqp_opts *)opts_;
    dense_qp_ooqp_memory *mem = (dense_qp_ooqp_memory *)memory_;
    dense_qp_ooqp_workspace *work = (dense_qp_ooqp_workspace *)work_;

    acados_tic(&interface_timer);
    // NOTE: has to be called after setting up the memory which contains the problem dimensions
    dense_qp_ooqp_cast_workspace(work, mem);
    dense_qp_ooqp_update_memory(qp_in, opts, mem);
    info->interface_time = acados_toc(&interface_timer);

    if (0)
    {
        print_dense_qp_in(qp_in);
        print_inputs(mem);
    }

    // TODO(dimitris): implement dense OOQP
    // call sparse OOQP
    acados_tic(&qp_timer);
    int ooqp_status;
    qpsolvede(mem->c,    mem->nx,      mem->dQ,
              mem->xlow, mem->ixlow,
              mem->xupp, mem->ixupp,
              mem->dA,   mem->my,      mem->bA,
              mem->dC,   mem->mz,
              mem->clow, mem->iclow,
              mem->cupp, mem->icupp,
              work->x,   work->gamma,  work->phi,
              work->y,
              work->z,   work->lambda, work->pi,
              &work->objectiveValue,
              opts->printLevel, &ooqp_status);
    info->solve_QP_time = acados_toc(&qp_timer);

    mem->time_qp_solver_call = info->solve_QP_time;
    mem->iter = -1;

    if (0) print_outputs(mem, work, ooqp_status);
    acados_tic(&interface_timer);
    fill_in_qp_out(qp_in, qp_out, work);
    dense_qp_compute_t(qp_in, qp_out);
    info->interface_time += acados_toc(&interface_timer);

    info->total_time = acados_toc(&tot_timer);
    info->num_iter = -1;
    info->t_computed = 1;

    int acados_status = ooqp_status;
    if (ooqp_status == DENSE_SUCCESSFUL_TERMINATION) acados_status = ACADOS_SUCCESS;
    if (ooqp_status == DENSE_MAX_ITS_EXCEEDED) acados_status = ACADOS_MAXITER;
    return acados_status;
}



void dense_qp_ooqp_destroy(void *mem_, void *work)
{
    free(work);
    free(mem_);
}



void dense_qp_ooqp_eval_sens(void *config_, void *qp_in, void *qp_out, void *opts_, void *mem_, void *work_)
{
    printf("\nerror: dense_qp_ooqp_eval_sens: not implemented yet\n");
    exit(1);
}

void dense_qp_ooqp_memory_reset(void *config, void *qp_in, void *qp_out, void *opts, void *mem, void *work)
{
    printf("\nerror: dense_qp_ooqp_memory_reset: not implemented yet\n");
    exit(1);
}

void dense_qp_ooqp_solver_get(void *config_, void *qp_in_, void *qp_out_, void *opts_, void *mem_, const char *field, int stage, void* value, int size1, int size2)
{
    printf("\nerror: dense_qp_ooqp_solver_get: not implemented yet\n");
    exit(1);
}

void dense_qp_ooqp_config_initialize_default(void *config_)
{
    qp_solver_config *config = config_;

    config->opts_calculate_size = (acados_size_t (*)(void *, void *)) & dense_qp_ooqp_opts_calculate_size;
    config->opts_assign = (void *(*) (void *, void *, void *) ) & dense_qp_ooqp_opts_assign;
    config->opts_initialize_default =
        (void (*)(void *, void *, void *)) & dense_qp_ooqp_opts_initialize_default;
    config->opts_update = (void (*)(void *, void *, void *)) & dense_qp_ooqp_opts_update;
    config->opts_set = &dense_qp_ooqp_opts_set;
    config->memory_calculate_size =
        (acados_size_t (*)(void *, void *, void *)) & dense_qp_ooqp_memory_calculate_size;
    config->memory_assign =
        (void *(*) (void *, void *, void *, void *) ) & dense_qp_ooqp_memory_assign;
    config->memory_get = &dense_qp_ooqp_memory_get;
    config->workspace_calculate_size =
        (acados_size_t (*)(void *, void *, void *)) & dense_qp_ooqp_workspace_calculate_size;
    config->evaluate = (int (*)(void *, void *, void *, void *, void *, void *)) & dense_qp_ooqp;
    config->eval_sens = &dense_qp_ooqp_eval_sens;
    config->memory_reset = &dense_qp_ooqp_memory_reset;
    config->solver_get = &dense_qp_ooqp_solver_get;
}
