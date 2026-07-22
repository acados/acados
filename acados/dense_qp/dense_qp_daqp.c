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
#include <math.h>
#include <string.h>
#include <stdio.h>
// blasfeo
#include "blasfeo_d_aux.h"
#include "blasfeo_d_blas.h"

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

#define DAQP_BLASFEO_MEM_ALIGNMENT 64

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
        // handled by the working set
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
    else if (!strcmp(field, "print_level"))
    {
        int *print_level = value;
        opts->print_level = *print_level;
    }
    else
    {
        printf("\nerror: dense_qp_daqp_opts_set: wrong field: %s\n", field);
        exit(1);
    }

    return;
}

void dense_qp_daqp_opts_get(void *config_, void *opts_, const char *field, void *value)
{
    // dense_qp_daqp_opts *opts = opts_;
    printf("\nerror: dense_qp_daqp_opts_get: not implemented for field: %s\n", field);
    exit(1);
}



/************************************************
 * memory
 ************************************************/

static acados_size_t daqp_workspace_calculate_size(int n, int m, int ms, int ns)
{
    acados_size_t size = 0;

    size += sizeof(DAQPWorkspace);

    size += n * (m-ms) * sizeof(c_float); // M
    size += 2 * m * sizeof(c_float); // dupper/dlower
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
    size += n * sizeof(int); // work->prox_mask
    size += (n+ns+1) * sizeof(int); // WS

    make_int_multiple_of(8, &size);

    return size;
}


acados_size_t dense_qp_daqp_memory_calculate_size(void *config_, dense_qp_dims *dims, void *opts_)
{
    int n = dims->nv;
    int m = dims->nv + dims->ng + dims->ne;
    int ms = 0;
    int nb = dims->nb;
    int ng = dims->ng;
    int ns = dims->ns;

    acados_size_t size = sizeof(dense_qp_daqp_memory);

    size += 2 * sizeof(struct blasfeo_dmat);
    size += 3 * sizeof(struct blasfeo_dvec);

    size += daqp_workspace_calculate_size(n, m, ms, ns);

    size += m * 2 * sizeof(c_float); // blower & bupper
    size += nb * 1 * sizeof(int); // idbx
    size += (nb + ng) * sizeof(int); // idxs_rev
    size += n *  1 * sizeof(int); // idxv_to_idxb;
    size += ns * 1 * sizeof(int); // idbs
    size += m  * 1 * sizeof(int); // QP-side sense snapshot

    size += ns * 6 * sizeof(c_float); // Zl,Zu,zl,zu,d_ls,d_us

    // Headroom for aligning the raw BLASFEO matrix storage below. This is
    // padding, not the size of a C object, so sizeof(...) does not apply.
    size += DAQP_BLASFEO_MEM_ALIGNMENT;
    size += blasfeo_memsize_dmat(n, n);
    size += blasfeo_memsize_dmat(n, nb + ng + dims->ne);
    size += 2 * blasfeo_memsize_dvec(n);
    size += blasfeo_memsize_dvec(m);
    make_int_multiple_of(8, &size);

    return size;
}


static void *daqp_workspace_assign(int n, int m, int ms, int ns, void *raw_memory)
{
    DAQPWorkspace *work;
    char *c_ptr = (char *) raw_memory;

    work = (DAQPWorkspace *) c_ptr;
    c_ptr += sizeof(DAQPWorkspace);

    align_char_to(8, &c_ptr);

    work->M = (c_float *) c_ptr;
    c_ptr += n * (m - ms) * sizeof(c_float);

    work->dupper = (c_float *) c_ptr;
    c_ptr += 1 * m * sizeof(c_float);

    work->dlower = (c_float *) c_ptr;
    c_ptr += 1 * m * sizeof(c_float);

    work->Rinv = NULL;

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

    work->prox_mask = (int *) c_ptr;
    c_ptr += n * sizeof(int);

    work->WS= (int *) c_ptr;
    c_ptr += (n+ns+1) * sizeof(int);

    // Initialize constants of workspace. The adapter supplies a complete LDP,
    // so DAQP never needs an accompanying DAQPProblem.
    work->qp = NULL;
    work->n = n;
    work->m = m;
    work->ms = ms;
    work->fval = -1;
    work->n_active = 0;
    work->iterations = 0;
    work->sing_ind  = 0;
    work->soft_slack = 0;

    work->RinvD = NULL;
    work->Mu = NULL;
    work->n_prox = 0;
    work->nh = 0;
    work->break_points = NULL;
    work->avi = NULL;
    work->timer = NULL;

    work->bnb = NULL; // No need to solve MIQP

    // initialize d_ls, d_us and sense
    for (int ii=0; ii<m; ii++)
    {
        work->d_ls[ii] = 0;
        work->d_us[ii] = 0;
        work->rho_ls[ii] = DAQP_DEFAULT_RHO_SOFT;
        work->rho_us[ii] = DAQP_DEFAULT_RHO_SOFT;
        work->sense[ii] = 0;
    }
    for (int ii=0; ii<n; ii++)
        work->prox_mask[ii] = 0;

    return work;
}


void *dense_qp_daqp_memory_assign(void *config_, dense_qp_dims *dims, void *opts_,
                                     void *raw_memory)
{
    dense_qp_daqp_memory *mem;

    int n = dims->nv;
    int m = dims->nv + dims->ng + dims->ne;
    int ms = 0;
    int nb = dims->nb;
    int ng = dims->ng;
    int ns = dims->ns;

    // char pointer
    char *c_ptr = (char *) raw_memory;

    mem = (dense_qp_daqp_memory *) c_ptr;
    c_ptr += sizeof(dense_qp_daqp_memory);

    mem->H_factor = (struct blasfeo_dmat *) c_ptr;
    c_ptr += sizeof(struct blasfeo_dmat);
    mem->M_factor = (struct blasfeo_dmat *) c_ptr;
    c_ptr += sizeof(struct blasfeo_dmat);
    mem->rhs_factor = (struct blasfeo_dvec *) c_ptr;
    c_ptr += sizeof(struct blasfeo_dvec);
    mem->v_factor = (struct blasfeo_dvec *) c_ptr;
    c_ptr += sizeof(struct blasfeo_dvec);
    mem->constraint_value = (struct blasfeo_dvec *) c_ptr;
    c_ptr += sizeof(struct blasfeo_dvec);
    mem->matrices_initialized = 0;

    assert((size_t) c_ptr % 8 == 0 && "memory not 8-byte aligned!");

    // Assign raw memory to workspace
    mem->daqp_work = daqp_workspace_assign(n, m, ms, ns, c_ptr);
    c_ptr += daqp_workspace_calculate_size(n, m, ms, ns);

    assert((size_t) c_ptr % 8 == 0 && "double not 8-byte aligned!");

    mem->blower = (c_float *) c_ptr;
    c_ptr += m * sizeof(c_float);

    mem->bupper = (c_float *) c_ptr;
    c_ptr += m * sizeof(c_float);

    mem->idxb = (int *) c_ptr;
    c_ptr += nb * 1 * sizeof(int);

    mem->idxs_rev = (int *) c_ptr;
    c_ptr += (nb + ng) * sizeof(int);


    mem->idxv_to_idxb = (int *) c_ptr;
    c_ptr += n * 1 * sizeof(int);

    mem->idxs= (int *) c_ptr;
    c_ptr += ns * 1 * sizeof(int);

    mem->sense = (int *) c_ptr;
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

    align_char_to(DAQP_BLASFEO_MEM_ALIGNMENT, &c_ptr);
    blasfeo_create_dmat(n, n, mem->H_factor, c_ptr);
    c_ptr += blasfeo_memsize_dmat(n, n);
    blasfeo_create_dmat(n, nb + ng + dims->ne, mem->M_factor, c_ptr);
    c_ptr += blasfeo_memsize_dmat(n, nb + ng + dims->ne);
    blasfeo_create_dvec(n, mem->rhs_factor, c_ptr);
    c_ptr += blasfeo_memsize_dvec(n);
    blasfeo_create_dvec(n, mem->v_factor, c_ptr);
    c_ptr += blasfeo_memsize_dvec(n);
    blasfeo_create_dvec(m, mem->constraint_value, c_ptr);
    c_ptr += blasfeo_memsize_dvec(m);

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


static void dense_qp_daqp_get_vectors(const dense_qp_in *qp, c_float *b, int *idxb,
        c_float *d_lg, c_float *d_ug, c_float *Zl, c_float *Zu,
        c_float *zl, c_float *zu, int *idxs, int *idxs_rev,
        c_float *d_ls, c_float *d_us)
{
    int nv = qp->dim->nv;
    int ne = qp->dim->ne;
    int nb = qp->dim->nb;
    int ng = qp->dim->ng;
    int ns = qp->dim->ns;

    if (ne > 0)
        blasfeo_unpack_dvec(ne, qp->b, 0, b, 1);
    if (nb > 0)
        for (int ii = 0; ii < nb; ii++)
            idxb[ii] = qp->idxb[ii];
    if (ng > 0)
    {
        blasfeo_unpack_dvec(ng, qp->d, nb, d_lg, 1);
        blasfeo_unpack_dvec(ng, qp->d, 2 * nb + ng, d_ug, 1);
        for (int ii = 0; ii < ng; ii++)
            d_ug[ii] = -d_ug[ii];
    }
    if (ns > 0)
    {
        for (int ii = 0; ii < nb + ng; ii++)
        {
            int idx_tmp = qp->idxs_rev[ii];
            idxs_rev[ii] = idx_tmp;
            if (idx_tmp != -1)
                idxs[idx_tmp] = ii;
        }
        blasfeo_unpack_dvec(ns, qp->Z, 0, Zl, 1);
        blasfeo_unpack_dvec(ns, qp->Z, ns, Zu, 1);
        blasfeo_unpack_dvec(ns, qp->gz, nv, zl, 1);
        blasfeo_unpack_dvec(ns, qp->gz, nv + ns, zu, 1);
        blasfeo_unpack_dvec(ns, qp->d, 2 * nb + 2 * ng, d_ls, 1);
        blasfeo_unpack_dvec(ns, qp->d, 2 * nb + 2 * ng + ns, d_us, 1);
    }
    else
    {
        for (int ii = 0; ii < nb + ng; ii++)
            idxs_rev[ii] = qp->idxs_rev[ii];
    }
}



static int dense_qp_daqp_update_memory(dense_qp_in *qp_in, const dense_qp_daqp_opts *opts, dense_qp_daqp_memory *mem)
{
    // extract dense qp size
    DAQPWorkspace * work = mem->daqp_work;
    int nv = qp_in->dim->nv;
    int nb = qp_in->dim->nb;
    int ns = qp_in->dim->ns;
    int ng = qp_in->dim->ng;
    int ne = qp_in->dim->ne;

    // extract daqp data
    double *blower = mem->blower;
    double *bupper = mem->bupper;
    int *idxb = mem->idxb;
    int *idxs = mem->idxs;
    int *sense = mem->sense;
    int update_matrices = opts->warm_start != 2 || !mem->matrices_initialized;

    // Build the QP-side sense input separately from DAQP's persistent workspace
    // state.  Only the dynamic warm-start bits are candidates for reuse; the
    // structural IMMUTABLE/SOFT bits are reconstructed from the current QP
    // below.  This is important when a bound disappears or the soft-constraint
    // mapping changes between solves.
    for (int ii = 0; ii < work->m; ii++)
        sense[ii] = opts->warm_start == 0 ? 0 :
            work->sense[ii] & (DAQP_ACTIVE | DAQP_LOWER);

    // Extract QP vectors and the compact list of actual bound indices before
    // forming the transformed constraint rows.
    dense_qp_daqp_get_vectors(qp_in,
        bupper+nv+ng,  // equalities
        idxb,  // bound indices
        blower+nv, bupper+nv,  // general linear constraints
        mem->Zl, mem->Zu, mem->zl, mem->zu, idxs, mem->idxs_rev, mem->d_ls, mem->d_us  // slacks
    );

    if (update_matrices)
    {
        // Form the LDP matrices in BLASFEO. The acados DAQP interface requires
        // the condensed Hessian to be positive definite, so H = L*L' can be
        // factorized directly without a DAQP-side fallback factorization.
        blasfeo_dpotrf_l(nv, qp_in->Hv, 0, 0, mem->H_factor, 0, 0);

        // Build only the actual constraint right-hand sides
        // [E_idxb C A'] and solve L*M' = [E_idxb C A']. This avoids forming
        // the full inverse of L when only a subset of its columns is needed.
        int n_constr = nb + ng + ne;
        if (n_constr > 0)
        {
            blasfeo_dgese(nv, n_constr, 0.0, mem->M_factor, 0, 0);
            for (int ii = 0; ii < nb; ii++)
                BLASFEO_DMATEL(mem->M_factor, idxb[ii], ii) = 1.0;
            if (ng > 0)
                blasfeo_dgecp(nv, ng, qp_in->Ct, 0, 0,
                        mem->M_factor, 0, nb);
            if (ne > 0)
                blasfeo_dgetr(ne, nv, qp_in->A, 0, 0,
                        mem->M_factor, 0, nb + ng);
            blasfeo_dtrsm_llnn(nv, n_constr, 1.0, mem->H_factor, 0, 0,
                    mem->M_factor, 0, 0, mem->M_factor, 0, 0);
        }
    }

    // "Unignore" all general linear inequalites (ng)
    for (int ii = nv; ii < nv+ng; ii++)
        sense[ii] &= ~DAQP_IMMUTABLE;

    // Setup upper/lower bounds
    for (int ii = 0; ii < nv; ii++)
    {
        // "ignore" bounds that are not in acados dense QP
        blower[ii] = -DAQP_INF;
        bupper[ii] = +DAQP_INF;
        // An absent bound must not remain active merely because this variable
        // was bounded in the previous QP.
        sense[ii] = DAQP_IMMUTABLE;
    }
    for (int ii = 0; ii < nb; ii++)
    {
        // "Unignore" bounds that are in acados dense QP and set bound values
        blower[idxb[ii]] = BLASFEO_DVECEL(qp_in->d, ii);
        bupper[idxb[ii]] = -BLASFEO_DVECEL(qp_in->d, nb + ng + ii);
        sense[idxb[ii]] &= ~DAQP_IMMUTABLE;
        mem->idxv_to_idxb[idxb[ii]] = ii;
    }

    // printf("DAQP: dmask\n");
    // blasfeo_print_tran_dvec(2*(nb+ng), qp_in->d_mask, 0);
    for (int ii = 0; ii < nb; ii++)
    {
        // NOTE: DAQP always works with double sided constraints, so currently one can not only ignore the upper bound or the lower bound.

        // "ignore" bounds that are marked as unconstrained in qp_in via dmask
        if (BLASFEO_DVECEL(qp_in->d_mask, ii) == 0.0)
        {
            blower[idxb[ii]] = -DAQP_INF;
        }
        if (BLASFEO_DVECEL(qp_in->d_mask, ii+ng+nb) == 0.0)
        {
            bupper[idxb[ii]] = +DAQP_INF;
        }
    }
    // ignore some general linear constraints.
    for (int ii = 0; ii < ng; ii++)
    {
        if (BLASFEO_DVECEL(qp_in->d_mask, nb+ii) == 0.0)
        {
            blower[ii+nv] = -DAQP_INF;
        }
        if (BLASFEO_DVECEL(qp_in->d_mask, 2*nb+ng+ii) == 0.0)
        {
            bupper[ii+nv] = +DAQP_INF;
        }
    }



    // Mark equality constraints
    for (int ii = 0; ii < ne; ii++)
    {
        // NOTE: b_eq values are ONLY in bupper, but sense status is default upper, thus fine.
        // Equalities are represented by bupper only, so always reactivate them
        // with upper sense regardless of their multiplier sign in the previous
        // solve.
        sense[nv+ng+ii] = DAQP_ACTIVE | DAQP_IMMUTABLE;
        // SET_ACTIVE(nv+ng+ii);
        // SET_IMMUTABLE(nv+ng+ii);
    }

    // Soft constraints
    int idxdaqp;  // index of soft constraint within DAQP ordering
    for (int ii = 0; ii < ns; ii++)
    {
        idxdaqp = idxs[ii] < nb ? idxb[idxs[ii]] : nv+idxs[ii]-nb;
        // DAQP's soft active-set state includes more than the bound side: it
        // also encodes whether the internal slack is fixed/free.  That state
        // is tied to the previous QP's normalized weights and slack bounds,
        // so do not reactivate soft constraints from only a partial snapshot.
        sense[idxdaqp] &= ~(DAQP_ACTIVE | DAQP_LOWER | DAQP_SLACK_FIXED);
        sense[idxdaqp] |= DAQP_SOFT;

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
        blower[idxdaqp]+=mem->zl[ii]/mem->Zl[ii];
        bupper[idxdaqp]-=mem->zu[ii]/mem->Zu[ii];

        work->d_ls[idxdaqp] = MAX(0,mem->zl[ii]+mem->Zl[ii]*mem->d_ls[ii]);
        work->d_us[idxdaqp] = MAX(0,mem->zu[ii]+mem->Zu[ii]*mem->d_us[ii]);

        // The default state in DAQP is that the soft slacks are active at their bounds
        // => shift bupper/blower with these bounds
        blower[idxdaqp] -= work->d_ls[idxdaqp]/mem->Zl[ii];
        bupper[idxdaqp] += work->d_us[idxdaqp]/mem->Zu[ii];
    }

    // From this point onward DAQP sees an LDP, not a QP. Preserve the active
    // set only for a true hot start; otherwise install the structural and
    // warm-start sense assembled above before validating the bounds.
    int do_activate = opts->warm_start != 2 || !mem->matrices_initialized;
    if (do_activate)
        memcpy(work->sense, sense, work->m * sizeof(int));

    int daqp_status = daqp_check_bounds(work, bupper, blower);
    if (daqp_status < 0)
        return daqp_status;
    do_activate |= daqp_status;

    if (update_matrices)
    {
        work->sing_ind = DAQP_EMPTY_IND;
        work->n_prox = 0;
        for (int ii = 0; ii < nv; ii++)
            work->prox_mask[ii] = 0;

        // All constraints are represented as general LDP rows. Unused bound
        // slots are immutable and never read; normalize and scatter only the
        // compact set of actual constraints.
        for (int ii = 0; ii < work->m; ii++)
            work->scaling[ii] = 1.0;

        int n_constr = nb + ng + ne;
        for (int ii = 0; ii < n_constr; ii++)
        {
            double norm_squared = 0.0;
            for (int jj = 0; jj < nv; jj++)
            {
                double value = BLASFEO_DMATEL(mem->M_factor, jj, ii);
                norm_squared += value * value;
            }
            int idx = ii < nb ? idxb[ii] : nv + ii - nb;
            if (norm_squared < work->settings->zero_tol)
            {
                work->scaling[idx] = 1.0;
#ifndef DAQP_ASSUME_VALID
                if ((bupper[idx] < -work->settings->zero_tol ||
                        blower[idx] > work->settings->zero_tol) &&
                        !DAQP_IS_IMMUTABLE(idx) && !DAQP_IS_SOFT(idx))
                    return DAQP_EXIT_INFEASIBLE;
#endif
                work->sense[idx] = DAQP_IMMUTABLE;
                continue;
            }
            work->scaling[idx] = 1.0 / sqrt(norm_squared);
            blasfeo_dcolsc(nv, work->scaling[idx], mem->M_factor, 0, ii);
        }

        for (int ii = 0; ii < nb; ii++)
            blasfeo_unpack_dmat(nv, 1, mem->M_factor, 0, ii,
                    work->M + nv * idxb[ii], nv);
        if (ng + ne > 0)
            blasfeo_unpack_dmat(nv, ng + ne, mem->M_factor, 0, nb,
                    work->M + nv * nv, nv);

        // Matrix-dependent active-set factorizations are no longer valid.
        reset_daqp_workspace(work);
        do_activate = 1;
    }

    // Transform the gradient with a triangular solve: v = inv(L)*g.
    blasfeo_dveccp(nv, qp_in->gz, 0, mem->rhs_factor, 0);
    blasfeo_dtrsv_lnn(nv, mem->H_factor, 0, 0,
            mem->rhs_factor, 0, mem->v_factor, 0);
    blasfeo_unpack_dvec(nv, mem->v_factor, 0, work->v, 1);

    // Compute the compact normalized constraint offsets G*v with BLASFEO.
    int n_constr = nb + ng + ne;
    if (n_constr > 0)
        blasfeo_dgemv_t(nv, n_constr, 1.0, mem->M_factor, 0, 0,
                mem->v_factor, 0, 0.0, mem->constraint_value, 0,
                mem->constraint_value, 0);
    for (int ii = 0; ii < work->m; ii++)
    {
        work->dupper[ii] = bupper[ii] * work->scaling[ii];
        work->dlower[ii] = blower[ii] * work->scaling[ii];
    }
    for (int ii = 0; ii < n_constr; ii++)
    {
        int idx = ii < nb ? idxb[ii] : nv + ii - nb;
        double offset = BLASFEO_DVECEL(mem->constraint_value, ii);
        work->dupper[idx] += offset;
        work->dlower[idx] += offset;
    }
    work->reuse_ind = 0;

#ifdef SOFT_WEIGHTS
    // Keep soft bounds and reciprocal quadratic weights in the normalized
    // constraint coordinates used by the LDP.
    for (int ii = 0; ii < work->m; ii++)
    {
        work->d_ls[ii] /= work->scaling[ii];
        work->d_us[ii] /= work->scaling[ii];
        work->rho_ls[ii] *= work->scaling[ii] * work->scaling[ii];
        work->rho_us[ii] *= work->scaling[ii] * work->scaling[ii];
    }
#endif

    if (do_activate)
    {
        reset_daqp_workspace(work);
        daqp_status = daqp_activate_constraints(work);
        if (daqp_status < 0)
            return daqp_status;
    }

    mem->matrices_initialized = 1;
    return 0;
}



// Map the normalized LDP solution back to the original dense-QP coordinates.
// This intentionally lives in the acados adapter: DAQP is only asked to solve
// the LDP assembled above.
static void dense_qp_daqp_ldp_to_qp_solution(dense_qp_daqp_memory *mem)
{
    DAQPWorkspace *work = mem->daqp_work;
    for (int ii = 0; ii < work->n; ii++)
        BLASFEO_DVECEL(mem->rhs_factor, ii) = work->u[ii] - work->v[ii];
    blasfeo_dtrsv_ltn(work->n, mem->H_factor, 0, 0,
            mem->rhs_factor, 0, mem->v_factor, 0);
    blasfeo_unpack_dvec(work->n, mem->v_factor, 0, work->x, 1);

    for (int ii = 0; ii < work->n_active; ii++)
        work->lam_star[ii] *= work->scaling[work->WS[ii]];
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
        mem->blower[idx]-=(mem->zl[i]-work->d_ls[idx]*work->scaling[idx])/mem->Zl[i];
        mem->bupper[idx]+=(mem->zu[i]-work->d_us[idx]*work->scaling[idx])/mem->Zu[i];

        c_float constraint_value;
        if (idx<nv)
            constraint_value = BLASFEO_DVECEL(v, idx);
        else
        {
            constraint_value = 0;
            for (int j = 0; j < nv; j++)
                constraint_value += BLASFEO_DMATEL(qp_in->Ct, j, idx-nv) * BLASFEO_DVECEL(v, j);
        }

        // Recover slacks from primal feasibility. This also handles soft
        // zero rows that DAQP can safely omit from its active-set system.
        BLASFEO_DVECEL(v, nv+i) = MAX(mem->d_ls[i], mem->blower[idx] - constraint_value);
        BLASFEO_DVECEL(lambda, 2*(nb+ng)+i) =
            mem->Zl[i] * BLASFEO_DVECEL(v, nv+i) + mem->zl[i]
            - BLASFEO_DVECEL(lambda, idxs[i]);

        BLASFEO_DVECEL(v, nv+ns+i) = MAX(mem->d_us[i], constraint_value - mem->bupper[idx]);
        BLASFEO_DVECEL(lambda, 2*(nb+ng)+ns+i) =
            mem->Zu[i] * BLASFEO_DVECEL(v, nv+ns+i) + mem->zu[i]
            - BLASFEO_DVECEL(lambda, idxs[i]+nb+ng);
    }
}



int dense_qp_daqp(void* config_, dense_qp_in *qp_in, dense_qp_out *qp_out, void *opts_, void *memory_, void *work_)
{
    qp_info *info = (qp_info *) qp_out->misc;
    acados_timer tot_timer, qp_timer, interface_timer;

    // Uncomment to print dense QPs to file for solution with matlab;
    // d_dense_qp_codegen_matlab("dense_qp_for_matlab.m", "a", qp_in->dim, qp_in);

    acados_tic(&tot_timer);
    acados_tic(&interface_timer);

    // cast structures
    dense_qp_daqp_opts *opts = (dense_qp_daqp_opts *) opts_;
    dense_qp_daqp_memory *memory = (dense_qp_daqp_memory *) memory_;

    // Extract workspace and update settings before forming the normalized LDP.
    DAQPWorkspace* work = memory->daqp_work;
    work->settings = opts->daqp_opts;
    if (opts->warm_start == 0)
        daqp_deactivate_constraints(work);

    // Form the complete normalized LDP in the acados adapter. A hot start
    // reuses its matrix part and only updates v and d.
    int daqp_status = dense_qp_daqp_update_memory(qp_in, opts, memory);
    info->interface_time = acados_toc(&interface_timer);
    if (daqp_status < 0)
        return daqp_status;

    // === Solve starts ===
    acados_tic(&qp_timer);
    daqp_status = daqp_ldp(work);
    dense_qp_daqp_ldp_to_qp_solution(memory);

    // extract primal and dual solution
    dense_qp_daqp_fill_output(memory,qp_out,qp_in);
    info->solve_QP_time = acados_toc(&qp_timer);

    // Constraint slacks are computed lazily by the residual routines when
    // needed. This matches the dense qpOASES interface and avoids an extra
    // pass over all constraints in the common solve-and-expand path.
    info->t_computed = 0;

    // log solve info
    info->total_time = acados_toc(&tot_timer);
    info->num_iter = memory->daqp_work->iterations;
    memory->time_qp_solver_call = info->solve_QP_time;
    memory->iter = memory->daqp_work->iterations;

    // status
    int acados_status = daqp_status;
    if (daqp_status == DAQP_EXIT_OPTIMAL || daqp_status == DAQP_EXIT_SOFT_OPTIMAL)
        acados_status = ACADOS_SUCCESS;
    else if (daqp_status == DAQP_EXIT_ITERLIMIT)
        acados_status = ACADOS_MAXITER;
    else if (daqp_status == DAQP_EXIT_INFEASIBLE)
        acados_status = ACADOS_INFEASIBLE;
    else if (daqp_status == DAQP_EXIT_UNBOUNDED)
        acados_status = ACADOS_UNBOUNDED;
    else
        acados_status = ACADOS_UNKNOWN;
    // NOTE: There are also:
    // EXIT_CYCLE, EXIT_UNBOUNDED, EXIT_NONCONVEX, EXIT_OVERDETERMINED_INITIAL
    return acados_status;
}


void dense_qp_daqp_eval_forw_sens(void *config_, void *qp_in, void *seed, void *qp_out, void *opts_, void *mem_, void *work_)
{
    printf("\nerror: dense_qp_daqp_eval_forw_sens: not implemented yet\n");
    exit(1);
}


void dense_qp_daqp_eval_adj_sens(void *config_, void *qp_in, void *seed, void *qp_out, void *opts_, void *mem_, void *work_)
{
    printf("\nerror: dense_qp_daqp_eval_adj_sens: not implemented yet\n");
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



void dense_qp_daqp_terminate(void *config_, void *mem_, void *work_)
{
    return;
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
    config->opts_get = &dense_qp_daqp_opts_get;
    config->memory_calculate_size =
        (acados_size_t (*)(void *, void *, void *)) & dense_qp_daqp_memory_calculate_size;
    config->memory_assign =
        (void *(*) (void *, void *, void *, void *) ) & dense_qp_daqp_memory_assign;
    config->memory_get = &dense_qp_daqp_memory_get;
    config->workspace_calculate_size =
        (acados_size_t (*)(void *, void *, void *)) & dense_qp_daqp_workspace_calculate_size;
    config->eval_forw_sens = &dense_qp_daqp_eval_forw_sens;
    config->eval_adj_sens = &dense_qp_daqp_eval_adj_sens;
    config->evaluate = (int (*)(void *, void *, void *, void *, void *, void *)) & dense_qp_daqp;
    config->memory_reset = &dense_qp_daqp_memory_reset;
    config->solver_get = &dense_qp_daqp_solver_get;
    config->terminate = &dense_qp_daqp_terminate;

    return;
}
