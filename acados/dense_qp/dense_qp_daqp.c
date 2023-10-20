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


// external
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <stdio.h>
// blasfeo
#include "blasfeo/include/blasfeo_d_aux.h"
#include "blasfeo/include/blasfeo_d_blas.h"

// daqp
#define SOFT_WEIGHTS
#include "daqp/include/types.h"
#include "daqp/include/api.h"
#include "daqp/include/daqp.h"
#include "daqp/include/utils.h"

// acados
#include "acados/dense_qp/dense_qp_common.h"
#include "acados/dense_qp/dense_qp_daqp.h"
#include "acados/utils/mem.h"
#include "acados/utils/timing.h"
#include "acados/utils/print.h"
#include "acados/utils/math.h"

#include "acados_c/dense_qp_interface.h"


/************************************************
 * opts
 ************************************************/

acados_size_t dense_qp_daqp_opts_calculate_size(void *config_, dense_qp_dims *dims)
{
    acados_size_t size = 0;
    size += sizeof(dense_qp_daqp_opts);
    size += sizeof(DAQPSettings);

    make_int_multiple_of(8, &size);

    return size;
}



void *dense_qp_daqp_opts_assign(void *config_, dense_qp_dims *dims, void *raw_memory)
{
    dense_qp_daqp_opts *opts;

    char *c_ptr = (char *) raw_memory;

    opts = (dense_qp_daqp_opts *) c_ptr;
    c_ptr += sizeof(dense_qp_daqp_opts);

    opts->daqp_opts = (DAQPSettings *) c_ptr;
    c_ptr += sizeof(DAQPSettings);

    assert((char *) raw_memory + dense_qp_daqp_opts_calculate_size(config_, dims) >= c_ptr);

    return (void *) opts;
}



void dense_qp_daqp_opts_initialize_default(void *config_, dense_qp_dims *dims, void *opts_)
{
    dense_qp_daqp_opts *opts = (dense_qp_daqp_opts *) opts_;
    daqp_default_settings(opts->daqp_opts);
    opts->warm_start=1;
    return;
}



void dense_qp_daqp_opts_update(void *config_, dense_qp_dims *dims, void *opts_)
{
    return;
}



void dense_qp_daqp_opts_set(void *config_, void *opts_, const char *field, void *value)
{
    dense_qp_daqp_opts *opts = opts_;
    if (!strcmp(field, "tol_stat"))
    {
        // DAQP always "aims" at a stationary point
    }
    else if (!strcmp(field, "tol_eq"))
    {
        // Equality constraints are explicitly
        // handled by the working set
    }
    else if (!strcmp(field, "tol_ineq"))
    {
        double *tol = value;
        opts->daqp_opts->primal_tol = *tol;
    }
    else if (!strcmp(field, "tol_comp"))
    {
        // Complementary slackness is implicitly
        // handled by the worlking set
    }
    else if (!strcmp(field, "iter_max"))
    {
        int *iter_max= value;
        opts->daqp_opts->iter_limit = *iter_max;
    }
    else if (!strcmp(field, "warm_start"))
    {
        int *warm_start = value;
        opts->warm_start = *warm_start;
    }
    else
    {
        printf("\nerror: dense_qp_daqp_opts_set: wrong field: %s\n", field);
        exit(1);
    }

    return;
}



/************************************************
 * memory
 ************************************************/

static acados_size_t daqp_workspace_calculate_size(int n, int m, int ms, int ns)
{
    acados_size_t size = 0;

    size += sizeof(DAQPWorkspace);
    size += sizeof(DAQPProblem);

    size += n * n * sizeof(c_float); // H
    size += 1 * n * sizeof(c_float); // f
    size += n * (m-ms) * sizeof(c_float); // A
    size += 2 * m * sizeof(c_float); // bupper/blower

    size += n * (m-ms) * sizeof(c_float); // M
    size += 2 * m * sizeof(c_float); // dupper/dlower
    size += (n+1)*n/2 * sizeof(c_float); // Rinv
    size += n * sizeof(c_float); // v
    size += m * sizeof(c_float); // scaling

    size += 2 * n * sizeof(c_float); // x & xold
    size += 2*(n+ns+1) * sizeof(c_float); // lam & lam_star
    size += n * sizeof(c_float); // u

    size += (n+ns+2)*(n+ns+1)/2 * sizeof(c_float); // L
    size += (n+ns+1) * sizeof(c_float); // D

    size += 2*(n+ns+1) * sizeof(c_float); //xldl & zldl

    size += 4 * m * sizeof(c_float); // d_ls, d_us,  rho_ls, rho_us

    size += m * sizeof(int); // work->sense
    size += (n+ns+1) * sizeof(int); // WS

    make_int_multiple_of(8, &size);

    return size;
}


acados_size_t dense_qp_daqp_memory_calculate_size(void *config_, dense_qp_dims *dims, void *opts_)
{
    int n = dims->nv;
    int m = dims->nv + dims->ng + dims->ne;
    int ms = dims->nv;
    int nb = dims->nb;
    int ns = dims->ns;

    acados_size_t size = sizeof(dense_qp_daqp_memory);

    size += daqp_workspace_calculate_size(n, m, ms, ns);

    size += nb * 2 * sizeof(c_float); // lb_tmp & ub_tmp
    size += nb * 1 * sizeof(int); // idbx
    size += n *  1 * sizeof(int); // idxv_to_idxb;
    size += ns * 1 * sizeof(int); // idbs
    size += m  * 1 * sizeof(int); // idxdaqp_to_idxs;

    size += ns * 6 * sizeof(c_float); // Zl,Zu,zl,zu,d_ls,d_us
    make_int_multiple_of(8, &size);

    return size;
}


static void *daqp_workspace_assign(int n, int m, int ms, int ns, void *raw_memory)
{
    DAQPWorkspace *work;
    char *c_ptr = (char *) raw_memory;

    work = (DAQPWorkspace *) c_ptr;
    c_ptr += sizeof(DAQPWorkspace);

    work->qp = (DAQPProblem *) c_ptr;
    c_ptr += sizeof(DAQPProblem);

    align_char_to(8, &c_ptr);

    // double
    work->qp->H = (c_float*) c_ptr;
    c_ptr += n * n * sizeof(c_float);

    work->qp->f = (c_float*) c_ptr;
    c_ptr += 1 * n * sizeof(c_float);

    work->qp->A = (c_float*) c_ptr;
    c_ptr += n * (m-ms) * sizeof(c_float);

    work->qp->bupper = (c_float*) c_ptr;
    c_ptr += 1 * m * sizeof(c_float);

    work->qp->blower = (c_float*) c_ptr;
    c_ptr += 1 * m * sizeof(c_float);

    work->M = (c_float *) c_ptr;
    c_ptr += n * (m - ms) * sizeof(c_float);

    work->dupper = (c_float *) c_ptr;
    c_ptr += 1 * m * sizeof(c_float);

    work->dlower = (c_float *) c_ptr;
    c_ptr += 1 * m * sizeof(c_float);

    work->Rinv = (c_float *) c_ptr;
    c_ptr += (n + 1) * n / 2 * sizeof(c_float);

    work->v = (c_float *) c_ptr;
    c_ptr += n * sizeof(c_float);

    work->scaling = (c_float *) c_ptr;
    c_ptr += m * sizeof(c_float);

    work->x = (c_float *) c_ptr;
    c_ptr += n * sizeof(c_float);

    work->xold = (c_float *) c_ptr;
    c_ptr += n * sizeof(c_float);

    work->lam = (c_float *) c_ptr;
    c_ptr += (n+ns+1) * sizeof(c_float);

    work->lam_star = (c_float *) c_ptr;
    c_ptr += (n+ns+1) * sizeof(c_float);

    work->u = (c_float *) c_ptr;
    c_ptr += n * sizeof(c_float);

    work->D = (c_float *) c_ptr;
    c_ptr += (n+ns+1) * sizeof(c_float);

    work->xldl = (c_float *) c_ptr;
    c_ptr += (n+ns+1) * sizeof(c_float);

    work->zldl = (c_float *) c_ptr;
    c_ptr += (n+ns+1) * sizeof(c_float);

    work->L = (c_float *) c_ptr;
    c_ptr += (n+ns+2)*(n+ns+1)/2 * sizeof(c_float);

    work->d_ls = (c_float *) c_ptr;
    c_ptr += m * sizeof(c_float);

    work->d_us = (c_float *) c_ptr;
    c_ptr += m * sizeof(c_float);

    work->rho_ls = (c_float *) c_ptr;
    c_ptr += m * sizeof(c_float);

    work->rho_us = (c_float *) c_ptr;
    c_ptr += m * sizeof(c_float);

    // ints
    work->sense = (int *) c_ptr;
    c_ptr += m * sizeof(int);

    work->WS= (int *) c_ptr;
    c_ptr += (n+ns+1) * sizeof(int);

    // Initialize constants of workspace
    work->qp->nb = 0;
    work->qp->bin_ids = NULL;

    work->n = n;
    work->m = m;
    work->ms = ms;
    work->fval = -1;
    work->n_active = 0;
    work->iterations = 0;
    work->sing_ind  = 0;
    work->soft_slack = 0;

    work->bnb = NULL; // No need to solve MIQP

    // initialize d_ls, d_us and sense
    for (int ii=0; ii<m; ii++)
    {
        work->d_ls[ii] = 0;
        work->d_us[ii] = 0;
        work->sense[ii] = 0;
    }

    return work;
}


void *dense_qp_daqp_memory_assign(void *config_, dense_qp_dims *dims, void *opts_,
                                     void *raw_memory)
{
    dense_qp_daqp_memory *mem;

    int n = dims->nv;
    int m = dims->nv + dims->ng + dims->ne;
    int ms = dims->nv;
    int nb = dims->nb;
    int ns = dims->ns;

    // char pointer
    char *c_ptr = (char *) raw_memory;

    mem = (dense_qp_daqp_memory *) c_ptr;
    c_ptr += sizeof(dense_qp_daqp_memory);

    assert((size_t) c_ptr % 8 == 0 && "memory not 8-byte aligned!");

    // Assign raw memory to workspace
    mem->daqp_work = daqp_workspace_assign(n, m, ms, ns, c_ptr);
    c_ptr += daqp_workspace_calculate_size(n, m, ms, ns);

    assert((size_t) c_ptr % 8 == 0 && "double not 8-byte aligned!");

    mem->lb_tmp = (c_float *) c_ptr;
    c_ptr += nb * 1 * sizeof(c_float);

    mem->ub_tmp = (c_float *) c_ptr;
    c_ptr += nb * 1 * sizeof(c_float);

    mem->idxb = (int *) c_ptr;
    c_ptr += nb * 1 * sizeof(int);

    mem->idxv_to_idxb = (int *) c_ptr;
    c_ptr += n * 1 * sizeof(int);

    mem->idxs= (int *) c_ptr;
    c_ptr += ns * 1 * sizeof(int);

    mem->idxdaqp_to_idxs = (int *) c_ptr;
    c_ptr += m * 1 * sizeof(int);

    mem->Zl = (c_float *) c_ptr;
    c_ptr += ns * 1 * sizeof(c_float);

    mem->Zu = (c_float *) c_ptr;
    c_ptr += ns * 1 * sizeof(c_float);

    mem->zl = (c_float *) c_ptr;
    c_ptr += ns * 1 * sizeof(c_float);

    mem->zu = (c_float *) c_ptr;
    c_ptr += ns * 1 * sizeof(c_float);

    mem->d_ls = (c_float *) c_ptr;
    c_ptr += ns * 1 * sizeof(c_float);

    mem->d_us = (c_float *) c_ptr;
    c_ptr += ns * 1 * sizeof(c_float);

    assert((char *) raw_memory + dense_qp_daqp_memory_calculate_size(config_, dims, opts_) >=
           c_ptr);

    return mem;
}



void dense_qp_daqp_memory_get(void *config_, void *mem_, const char *field, void* value)
{
    // qp_solver_config *config = config_;
    dense_qp_daqp_memory *mem = mem_;

    if (!strcmp(field, "time_qp_solver_call"))
    {
        double *tmp_ptr = value;
        *tmp_ptr = mem->time_qp_solver_call;
    }
    else if (!strcmp(field, "iter"))
    {
        int *tmp_ptr = value;
        *tmp_ptr = mem->iter;
    }
    else
    {
        printf("\nerror: dense_qp_daqp_memory_get: field %s not available\n", field);
        exit(1);
    }

    return;

}

/************************************************
 * workspace
 ************************************************/

acados_size_t dense_qp_daqp_workspace_calculate_size(void *config_, dense_qp_dims *dims, void *opts_)
{
    return 0;
}


/************************************************
 * functions
 ************************************************/

// NOTE on transcription of acados dense QP into DAQP formulation:


// DAQP constraints are: [bounds on ALL x; linear constraints (ng); equality constraints (ne)]

// A_DAQP = [C; A]
// blower = [lower bounds on ALL x (-INF if not set); lg; b_eq]
// bupper = [upper bounds on ALL x (+INF if not set); ug; b_eq]



static void dense_qp_daqp_update_memory(dense_qp_in *qp_in, const dense_qp_daqp_opts *opts, dense_qp_daqp_memory *mem)
{
    // extract dense qp size
    DAQPWorkspace * work = mem->daqp_work;
    int nv = qp_in->dim->nv;
    int nb = qp_in->dim->nb;
    int ns = qp_in->dim->ns;
    int ng = qp_in->dim->ng;
    int ne = qp_in->dim->ne;

    // extract daqp data
    double *lb_tmp = mem->lb_tmp;
    double *ub_tmp = mem->ub_tmp;
    int *idxb = mem->idxb;
    int *idxs = mem->idxs;

    // fill in the upper triangular of H in dense_qp
    blasfeo_dtrtr_l(nv, qp_in->Hv, 0, 0, qp_in->Hv, 0, 0);

    // extract data from qp_in in row-major
    d_dense_qp_get_all_rowmaj(qp_in, work->qp->H, work->qp->f,  // objective
        work->qp->A+nv*ng, work->qp->bupper+nv+ng,  // equalities
        idxb, lb_tmp, ub_tmp,  // bounds
        work->qp->A, work->qp->blower+nv, work->qp->bupper+nv,  // general linear constraints
        mem->Zl, mem->Zu, mem->zl, mem->zu, idxs, mem->d_ls, mem->d_us  // slacks
    );

    // printf("\nDAQP: matrix A\n");
    // int m = qp_in->dim->nv + qp_in->dim->ng + qp_in->dim->ne;
    // int ms = qp_in->dim->nv;
    // d_print_exp_tran_mat(nv, m-ms, work->qp->A, nv);

    // "Unignore" all general linear inequalites (ng)
    for (int ii = nv; ii < nv+ng; ii++)
        SET_MUTABLE(ii);

    // Setup upper/lower bounds
    for (int ii = 0; ii < nv; ii++)
    {
        // "ignore" bounds that are not in acados dense QP
        work->qp->blower[ii] = -DAQP_INF;
        work->qp->bupper[ii] = +DAQP_INF;
        SET_IMMUTABLE(ii);
    }
    for (int ii = 0; ii < nb; ii++)
    {
        // "Unignore" bounds that are in acados dense QP and set bound values
        work->qp->blower[idxb[ii]] = lb_tmp[ii];
        work->qp->bupper[idxb[ii]] = ub_tmp[ii];
        SET_MUTABLE(idxb[ii]);
        mem->idxv_to_idxb[idxb[ii]] = ii;
    }

    // Mark equality constraints
    for (int ii = 0; ii < ne; ii++)
    {
        // NOTE: b_eq values are ONLY in bupper, but sense status is default upper, thus fine.
        work->sense[nv+ng+ii] &= ACTIVE+IMMUTABLE;
        // SET_ACTIVE(nv+ng+ii);
        // SET_IMMUTABLE(nv+ng+ii);
    }

    // Soft constraints
    int idxdaqp;  // index of soft constraint within DAQP ordering
    for (int ii = 0; ii < ns; ii++)
    {
        idxdaqp = idxs[ii] < nb ? idxb[idxs[ii]] : nv+idxs[ii]-nb;
        mem->idxdaqp_to_idxs[idxdaqp] = ii;

        SET_SOFT(idxdaqp);

        // Quadratic slack penalty needs to be nonzero in DAQP
        mem->Zl[ii] = MAX(1e-8,mem->Zl[ii]);
        mem->Zu[ii] = MAX(1e-8,mem->Zu[ii]);

        // Setup soft weight used in DAQP
        work->rho_ls[idxdaqp] = 1/mem->Zl[ii];
        work->rho_us[idxdaqp] = 1/mem->Zu[ii];

        // Shift QP to handle linear terms on slack
        // DAQP penalizes soft slacks s using a quadratic penalty s' s, bounded by s >= d_l
        // (instead of weighting the soft slacks in the objective as in acados,
        // the soft slacks are weighted in the constraint by rho=1/Z).
        // To remove the linear term from acados we use the transformation
        //             s_daqp = (Z*s_acados+z/Z),
        // which will shift blower/bupper and scale the nominal slack bounds with 1/Z
        work->qp->blower[idxdaqp]+=mem->zl[ii]/mem->Zl[ii];
        work->qp->bupper[idxdaqp]-=mem->zu[ii]/mem->Zu[ii];

        work->d_ls[idxdaqp] = MAX(0,mem->zl[ii]+mem->Zl[ii]*mem->d_ls[ii]);
        work->d_us[idxdaqp] = MAX(0,mem->zu[ii]+mem->Zu[ii]*mem->d_us[ii]);

        // The default state in DAQP is that the soft slacks are active at their bounds
        // => shift bupper/blower with these bounds
        work->qp->blower[idxdaqp] -= work->d_ls[idxdaqp]/mem->Zl[ii];
        work->qp->bupper[idxdaqp] += work->d_us[idxdaqp]/mem->Zu[ii];
    }
}



static void dense_qp_daqp_fill_output(dense_qp_daqp_memory *mem, const dense_qp_out *qp_out, const dense_qp_in *qp_in)
{
    int *idxv_to_idxb = mem->idxv_to_idxb;
    int *idxs = mem->idxs;
    int *idxb = mem->idxb;
    int i;
    int nv = qp_in->dim->nv;
    int nb = qp_in->dim->nb;
    int ng = qp_in->dim->ng;
    int ns = qp_in->dim->ns;
    DAQPWorkspace *work = mem->daqp_work;

    struct blasfeo_dvec *v = qp_out->v;
    struct blasfeo_dvec *lambda = qp_out->lam;

    // print DAQP solution before expansion:
    // printf("\n\nDAQP solution\n");
    // printf("------------------\n");
    // printf("\nx (primals):\n\n");
    // for (i = 0; i<nv; i++)
    //     printf("%e\t", work->x[i]);
    // printf("\nlambda (duals):\n\n");
    // for (i = 0; i<work->n_active; i++)
    //     printf("%e\t", work->lam_star[i]);
    // printf("\n\n");

    // primal variables
    blasfeo_pack_dvec(nv, work->x, 1, v, 0);


    // dual variables
    blasfeo_dvecse(2 * nb + 2 * ng + 2 * ns, 0.0, lambda, 0);
    c_float lam;
    for (i = 0; i < work->n_active; i++)
    {
        lam = work->lam_star[i];
        if (work->WS[i] < nv) // bound constraint
        {
            if (lam >= 0.0)
                BLASFEO_DVECEL(lambda, nb+ng+idxv_to_idxb[work->WS[i]]) = lam;
            else
                BLASFEO_DVECEL(lambda, idxv_to_idxb[work->WS[i]]) = -lam;
        }
        else if (work->WS[i] < nv+ng)// general constraint
        {
            if (lam >= 0.0)
                BLASFEO_DVECEL(lambda, 2*nb+ng+work->WS[i]-nv) = lam;
            else
                BLASFEO_DVECEL(lambda, nb+work->WS[i]-nv) = -lam;
        }
        else // equality constraint
            BLASFEO_DVECEL(qp_out->pi, work->WS[i]-nv-ng) = lam;
    }

    // soft slacks
    int idx;
    for (i = 0; i < ns; i++)
    {
        idx = idxs[i] < nb ? idxb[idxs[i]] : nv+idxs[i]-nb;
        // shift back QP
        work->qp->blower[idx]-=(mem->zl[i]-work->d_ls[idx]/work->scaling[idx])/mem->Zl[i];
        work->qp->bupper[idx]+=(mem->zu[i]-work->d_us[idx]/work->scaling[idx])/mem->Zu[i];

        // lower
        if (BLASFEO_DVECEL(lambda, idxs[i]) == 0) // inactive soft => active slack bound
        {
            BLASFEO_DVECEL(v, nv+i) = mem->d_ls[i];
            BLASFEO_DVECEL(lambda, 2*nb+2*ng+i) = mem->d_ls[i]/mem->Zl[i]+mem->zl[i];
        }
        else
        { // if soft active => compute slack directly from equality
            BLASFEO_DVECEL(v, nv+i) = work->qp->blower[idx];
            if (idx<nv)
                BLASFEO_DVECEL(v, nv+i) -= BLASFEO_DVECEL(v, idx);
            else
            { // general constraint
                for (int j=0, disp = (idx-nv)*nv; j < nv; j++, disp++)
                {
                    BLASFEO_DVECEL(v, nv+i) -= work->qp->A[disp] * BLASFEO_DVECEL(v, j);
                }
            }
            // compute dual variable from stationarity condition
            BLASFEO_DVECEL(lambda, 2*(nb+ng)+i) = mem->Zl[i] * BLASFEO_DVECEL(v, nv+i) + mem->zl[i]
                - BLASFEO_DVECEL(lambda, idxs[i]);
        }

        // upper
        if (BLASFEO_DVECEL(lambda, idxs[i]+nb+ng) == 0) // inactive soft => active slack bound
        {
            BLASFEO_DVECEL(v, nv+ns+i) = mem->d_us[i];
            BLASFEO_DVECEL(lambda, 2*nb+2*ng+ns+i) = mem->d_us[i]/mem->Zu[i]+mem->zu[i];
        }
        else
        { // if soft active => compute slack directly from equality
            BLASFEO_DVECEL(v, nv+ns+i) = -work->qp->bupper[idx];
            if (idx<nv)
                BLASFEO_DVECEL(v, nv+ns+i) += BLASFEO_DVECEL(v, idx);
            else
            { // general constraint
                for (int j=0, disp = (idx-nv)*nv; j < nv; j++, disp++)
                    BLASFEO_DVECEL(v, nv+ns+i) += work->qp->A[disp] * BLASFEO_DVECEL(v, j);
            }
            // compute dual variable from stationarity condition
            BLASFEO_DVECEL(lambda, 2*(nb+ng)+ns+i) = mem->Zu[i] * BLASFEO_DVECEL(v, nv+ns+i) + mem->zu[i]
                - BLASFEO_DVECEL(lambda, idxs[i]+nb+ng);
        }
    }
}



int dense_qp_daqp(void* config_, dense_qp_in *qp_in, dense_qp_out *qp_out, void *opts_, void *memory_, void *work_)
{
    qp_info *info = (qp_info *) qp_out->misc;
    acados_timer tot_timer, qp_timer, interface_timer;

    acados_tic(&tot_timer);
    acados_tic(&interface_timer);

    // cast structures
    dense_qp_daqp_opts *opts = (dense_qp_daqp_opts *) opts_;
    dense_qp_daqp_memory *memory = (dense_qp_daqp_memory *) memory_;

    // Move data into daqp workspace
    dense_qp_daqp_update_memory(qp_in,opts,memory);
    info->interface_time = acados_toc(&interface_timer);

    // Extract workspace and update settings
    DAQPWorkspace* work = memory->daqp_work;
    work->settings = opts->daqp_opts;

    // === Solve starts ===
    acados_tic(&qp_timer);
    if (opts->warm_start==0) deactivate_constraints(work);
    // setup LDP
    int update_mask,daqp_status;
    update_mask= (opts->warm_start==2) ?
        UPDATE_v+UPDATE_d: UPDATE_Rinv+UPDATE_M+UPDATE_v+UPDATE_d;
    daqp_status = update_ldp(update_mask,work);
    // if setup failed, abort
    if(daqp_status < 0)
        return daqp_status;
    // solve LDP
    if (opts->warm_start==1)
        activate_constraints(work);

    // TODO: shift active set? - not in SQP but would be nice as an option in SQP_RTI.

    daqp_status = daqp_ldp(memory->daqp_work);
    ldp2qp_solution(work);

    // extract primal and dual solution
    dense_qp_daqp_fill_output(memory,qp_out,qp_in);
    info->solve_QP_time = acados_toc(&qp_timer);

    acados_tic(&interface_timer);

    // compute slacks
    dense_qp_compute_t(qp_in, qp_out);
    info->t_computed = 1;

    // log solve info
    info->interface_time += acados_toc(&interface_timer);
    info->total_time = acados_toc(&tot_timer);
    info->num_iter = memory->daqp_work->iterations;
    memory->time_qp_solver_call = info->solve_QP_time;
    memory->iter = memory->daqp_work->iterations;

    // status
    int acados_status = daqp_status;
    if (daqp_status == EXIT_OPTIMAL || daqp_status == EXIT_SOFT_OPTIMAL)
        acados_status = ACADOS_SUCCESS;
    else if (daqp_status == EXIT_ITERLIMIT)
        acados_status = ACADOS_MAXITER;
    // NOTE: There are also:
    // EXIT_INFEASIBLE, EXIT_CYCLE, EXIT_UNBOUNDED, EXIT_NONCONVEX, EXIT_OVERDETERMINED_INITIAL

    return acados_status;
}


void dense_qp_daqp_eval_sens(void *config_, void *qp_in, void *qp_out, void *opts_, void *mem_, void *work_)
{
    printf("\nerror: dense_qp_daqp_eval_sens: not implemented yet\n");
    exit(1);
}

void dense_qp_daqp_memory_reset(void *config, void *qp_in, void *qp_out, void *opts, void *mem, void *work)
{
    printf("\nerror: dense_qp_daqp_memory_reset: not implemented yet\n");
    exit(1);
}

void dense_qp_daqp_solver_get(void *config_, void *qp_in_, void *qp_out_, void *opts_, void *mem_, const char *field, int stage, void* value, int size1, int size2)
{
    printf("\nerror: dense_qp_daqp_solver_get: not implemented yet\n");
    exit(1);
}


void dense_qp_daqp_config_initialize_default(void *config_)
{
    qp_solver_config *config = config_;

    config->opts_calculate_size = (acados_size_t (*)(void *, void *)) & dense_qp_daqp_opts_calculate_size;
    config->opts_assign = (void *(*) (void *, void *, void *) ) & dense_qp_daqp_opts_assign;
    config->opts_initialize_default =
        (void (*)(void *, void *, void *)) & dense_qp_daqp_opts_initialize_default;
    config->opts_update = (void (*)(void *, void *, void *)) & dense_qp_daqp_opts_update;
    config->opts_set = &dense_qp_daqp_opts_set;
    config->memory_calculate_size =
        (acados_size_t (*)(void *, void *, void *)) & dense_qp_daqp_memory_calculate_size;
    config->memory_assign =
        (void *(*) (void *, void *, void *, void *) ) & dense_qp_daqp_memory_assign;
    config->memory_get = &dense_qp_daqp_memory_get;
    config->workspace_calculate_size =
        (acados_size_t (*)(void *, void *, void *)) & dense_qp_daqp_workspace_calculate_size;
    config->eval_sens = &dense_qp_daqp_eval_sens;
    config->evaluate = (int (*)(void *, void *, void *, void *, void *, void *)) & dense_qp_daqp;
    config->memory_reset = &dense_qp_daqp_memory_reset;
    config->solver_get = &dense_qp_daqp_solver_get;

    return;
}
