/*
 * Copyright (c) The acados authors.
 *
 * This file is part of acados.
 *
 * Licensed under the 2-Clause BSD License.
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

    size += n * sizeof(c_float); // x
    size += 2*(n+ns+1) * sizeof(c_float); // lam & lam_star
    size += n * sizeof(c_float); // u

    size += (n+ns+2)*(n+ns+1)/2 * sizeof(c_float); // L
    size += (n+ns+1) * sizeof(c_float); // D

    size += 2*(n+ns+1) * sizeof(c_float); //xldl & zldl

    size += 4 * m * sizeof(c_float); // w_ls, w_us, rho_ls, rho_us

    size += m * sizeof(int); // work->sense
    size += (n+ns+1) * sizeof(int); // WS

    make_int_multiple_of(8, &size);

    return size;
}


acados_size_t dense_qp_daqp_memory_calculate_size(void *config_, dense_qp_dims *dims, void *opts_)
{
    int n = dims->nv;
    int m = dims->nb + dims->ng + dims->ne;
    int ms = 0;  // (M also holds the packed rows of the simple bounds)
    int ns = dims->ns;

    acados_size_t size = sizeof(dense_qp_daqp_memory);

    size += 2 * sizeof(struct blasfeo_dmat);
    size += 3 * sizeof(struct blasfeo_dvec);

    size += daqp_workspace_calculate_size(n, m, ms, ns);

    size += m * 2 * sizeof(c_float); // blower & bupper
    size += ns * 1 * sizeof(int); // idbs

    size += ns * 6 * sizeof(c_float); // Zl,Zu,zl,zu,d_ls,d_us
    size += m * sizeof(c_float); // Mu

    // Headroom for aligning the raw BLASFEO matrix storage below. This is
    // padding, not the size of a C object, so sizeof(...) does not apply.
    size += DAQP_BLASFEO_MEM_ALIGNMENT;
    size += blasfeo_memsize_dmat(n, n);
    size += blasfeo_memsize_dmat(n, m);
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
    memset(work, 0, sizeof(*work));

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

    // The acados adapter does not call DAQP's proximal wrappers.
    work->xold = NULL;

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

    work->w_ls = (c_float *) c_ptr;
    c_ptr += m * sizeof(c_float);

    work->w_us = (c_float *) c_ptr;
    c_ptr += m * sizeof(c_float);

    work->rho_ls = (c_float *) c_ptr;
    c_ptr += m * sizeof(c_float);

    work->rho_us = (c_float *) c_ptr;
    c_ptr += m * sizeof(c_float);

    // ints
    work->sense = (int *) c_ptr;
    c_ptr += m * sizeof(int);

    work->prox_mask = NULL;

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
    work->eq = NULL;
    work->state = 0;
    work->timer = NULL;

    work->bnb = NULL; // No need to solve MIQP

    // initialize linear soft weights and sense
    for (int ii=0; ii<m; ii++)
    {
        work->w_ls[ii] = 0;
        work->w_us[ii] = 0;
        work->rho_ls[ii] = DAQP_DEFAULT_RHO_SOFT;
        work->rho_us[ii] = DAQP_DEFAULT_RHO_SOFT;
        work->sense[ii] = 0;
    }
    return work;
}


void *dense_qp_daqp_memory_assign(void *config_, dense_qp_dims *dims, void *opts_,
                                     void *raw_memory)
{
    dense_qp_daqp_memory *mem;

    int n = dims->nv;
    int m = dims->nb + dims->ng + dims->ne;
    int ms = 0;  // (M also holds the packed rows of the simple bounds)
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

    mem->daqp_work->Mu = (c_float *) c_ptr;
    c_ptr += m * sizeof(c_float);

    // (ints after the doubles, to keep these aligned)
    mem->idxs = (int *) c_ptr;
    c_ptr += ns * 1 * sizeof(int);

    mem->ldp_rows = mem->daqp_work->M;

    align_char_to(DAQP_BLASFEO_MEM_ALIGNMENT, &c_ptr);
    blasfeo_create_dmat(n, n, mem->H_factor, c_ptr);
    c_ptr += blasfeo_memsize_dmat(n, n);
    blasfeo_create_dmat(n, m, mem->M_factor, c_ptr);
    c_ptr += blasfeo_memsize_dmat(n, m);
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


// DAQP constraints are compactly ordered as:
// [actual variable bounds (nb); linear constraints (ng); equalities (ne)].
// The leading bounds with idxb[ii] == ii are DAQP's simple bounds (rows of
// Rinv); the remaining ones are general LDP rows.


static void dense_qp_daqp_get_vectors(const dense_qp_in *qp, c_float *b,
        c_float *d_lg, c_float *d_ug, c_float *Zl, c_float *Zu,
        c_float *zl, c_float *zu, int *idxs,
        c_float *d_ls, c_float *d_us)
{
    int nv = qp->dim->nv;
    int ne = qp->dim->ne;
    int nb = qp->dim->nb;
    int ng = qp->dim->ng;
    int ns = qp->dim->ns;

    if (ne > 0)
        blasfeo_unpack_dvec(ne, qp->b, 0, b, 1);
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
    int *idxb = qp_in->idxb;
    int *idxs = mem->idxs;
    int update_matrices = opts->warm_start != 2 || !mem->matrices_initialized;
    // Reuse the complete soft active-set state when its effective penalties
    // are unchanged. A free slack contributes rho to the LDL diagonal;
    // changing the penalties or soft-row mapping requires reactivation.
    int do_activate = update_matrices;
    if (!do_activate && ns > 0)
        for (int ii = 0; ii < nb + ng; ii++)
        {
            int si = qp_in->idxs_rev[ii];
            if (si < 0)
                continue;
            double Zl = MAX(1e-8, BLASFEO_DVECEL(qp_in->Z, si));
            double Zu = MAX(1e-8, BLASFEO_DVECEL(qp_in->Z, ns + si));
            double zl = BLASFEO_DVECEL(qp_in->gz, nv + si);
            double zu = BLASFEO_DVECEL(qp_in->gz, nv + ns + si);
            double sl = MAX(BLASFEO_DVECEL(qp_in->d, 2*(nb+ng) + si), -zl/Zl);
            double su = MAX(BLASFEO_DVECEL(qp_in->d, 2*(nb+ng) + ns + si), -zu/Zu);
            if (idxs[si] != ii || work->rho_ls[ii] != 1/Zl ||
                    work->rho_us[ii] != 1/Zu ||
                    work->w_ls[ii] != MAX(0, zl + Zl*sl) ||
                    work->w_us[ii] != MAX(0, zu + Zu*su))
            {
                do_activate = 1;
                break;
            }
        }

    // Retain dynamic warm-start bits before reconstructing structural state.
    // With fixed matrices, also preserve immutable zero rows.
    if (do_activate)
        for (int ii = 0; ii < work->m; ii++)
            work->sense[ii] = opts->warm_start == 0 ? 0 :
                work->sense[ii] & (DAQP_ACTIVE | DAQP_LOWER |
                        (update_matrices ? 0 : DAQP_IMMUTABLE));

    // Extract QP vectors and the compact list of actual bound indices before
    // forming the transformed constraint rows.
    dense_qp_daqp_get_vectors(qp_in,
        bupper+nb+ng,  // equalities
        blower+nb, bupper+nb,  // general linear constraints
        mem->Zl, mem->Zu, mem->zl, mem->zu, idxs, mem->d_ls, mem->d_us  // slacks
    );

    if (update_matrices)
    {
        // Simple bounds: the leading bounds on their own index
        int ms = 0;
        while (ms < nb && idxb[ms] == ms)
            ms++;
        work->ms = ms;
        work->Rinv = ms > 0 ? mem->ldp_rows : NULL;
        work->M = mem->ldp_rows + ms*nv - ms*(ms-1)/2;

        blasfeo_dpotrf_l(nv, qp_in->Hv, 0, 0, mem->H_factor, 0, 0);

        if (nb > 0)
            blasfeo_dgese(nv, nb, 0.0, mem->M_factor, 0, 0);
        for (int ii = 0; ii < nb; ii++)
            BLASFEO_DMATEL(mem->M_factor, idxb[ii], ii) = 1.0;
        if (ng > 0)
            blasfeo_dgecp(nv, ng, qp_in->Ct, 0, 0,
                    mem->M_factor, 0, nb);
        if (ne > 0)
            blasfeo_dgetr(ne, nv, qp_in->A, 0, 0,
                    mem->M_factor, 0, nb + ng);
        if (work->m > 0)
            blasfeo_dtrsm_llnn(nv, work->m, 1.0, mem->H_factor, 0, 0,
                    mem->M_factor, 0, 0, mem->M_factor, 0, 0);
    }

    // Setup upper/lower bounds
    for (int ii = 0; ii < nb; ii++)
    {
        blower[ii] = BLASFEO_DVECEL(qp_in->d, ii);
        bupper[ii] = -BLASFEO_DVECEL(qp_in->d, nb + ng + ii);
    }

    // printf("DAQP: dmask\n");
    // blasfeo_print_tran_dvec(2*(nb+ng), qp_in->d_mask, 0);
    for (int ii = 0; ii < nb; ii++)
    {
        // NOTE: DAQP always works with double sided constraints, so currently one can not only ignore the upper bound or the lower bound.

        // "ignore" bounds that are marked as unconstrained in qp_in via dmask
        if (BLASFEO_DVECEL(qp_in->d_mask, ii) == 0.0)
        {
            blower[ii] = -DAQP_INF;
        }
        if (BLASFEO_DVECEL(qp_in->d_mask, ii+ng+nb) == 0.0)
        {
            bupper[ii] = +DAQP_INF;
        }
    }
    // ignore some general linear constraints.
    for (int ii = 0; ii < ng; ii++)
    {
        if (BLASFEO_DVECEL(qp_in->d_mask, nb+ii) == 0.0)
        {
            blower[ii+nb] = -DAQP_INF;
        }
        if (BLASFEO_DVECEL(qp_in->d_mask, 2*nb+ng+ii) == 0.0)
        {
            bupper[ii+nb] = +DAQP_INF;
        }
    }



    // Mark equality constraints
    for (int ii = 0; ii < ne; ii++)
    {
        // NOTE: b_eq values are ONLY in bupper, but sense status is default upper, thus fine.
        // Equalities are represented by bupper only, so always reactivate them
        // with upper sense regardless of their multiplier sign in the previous
        // solve.
        blower[nb+ng+ii] = -DAQP_INF;
        if (do_activate)
            work->sense[nb+ng+ii] = DAQP_ACTIVE | DAQP_IMMUTABLE;
    }

    // Soft constraints
    int idxdaqp;  // index of soft constraint within DAQP ordering
    for (int ii = 0; ii < ns; ii++)
    {
        idxdaqp = idxs[ii];
        // DAQP's soft active-set state includes more than the bound side: it
        // also encodes whether the internal slack is fixed/free.  That state
        // is tied to the previous QP's weights and shifted bounds,
        // so do not reactivate soft constraints from only a partial snapshot.
        if (do_activate)
        {
            work->sense[idxdaqp] &= ~(DAQP_ACTIVE | DAQP_LOWER | DAQP_SLACK_FIXED);
            work->sense[idxdaqp] |= DAQP_SOFT;
        }

        // Quadratic slack penalty needs to be nonzero in DAQP
        mem->Zl[ii] = MAX(1e-8,mem->Zl[ii]);
        mem->Zu[ii] = MAX(1e-8,mem->Zu[ii]);

        // Setup soft weight used in DAQP
        work->rho_ls[idxdaqp] = 1/mem->Zl[ii];
        work->rho_us[idxdaqp] = 1/mem->Zu[ii];

        // Translate s >= s_min to a nonnegative DAQP slack t = s - s_min.
        // Its penalty is 0.5*Z*t^2 + (z + Z*s_min)*t, up to a constant.
        // If that linear weight is negative, shift to the unconstrained slack
        // minimizer instead; feasibility then guarantees s >= s_min.
        double sl = MAX(mem->d_ls[ii], -mem->zl[ii]/mem->Zl[ii]);
        double su = MAX(mem->d_us[ii], -mem->zu[ii]/mem->Zu[ii]);
        blower[idxdaqp] -= sl;
        bupper[idxdaqp] += su;
        work->w_ls[idxdaqp] = MAX(0, mem->zl[ii] + mem->Zl[ii]*sl);
        work->w_us[idxdaqp] = MAX(0, mem->zu[ii] + mem->Zu[ii]*su);
    }

    // remove constraints which are not present anymore, e.g. masked ones, from the warm start guess
    if (opts->warm_start == 1)
    {
        for (int ii = 0; ii < nb + ng; ii++)
        {
            if (DAQP_IS_ACTIVE(ii) && !DAQP_IS_IMMUTABLE(ii))
            {
                if (DAQP_IS_LOWER(ii) && blower[ii] <= -DAQP_INF)
                {
                    DAQP_SET_INACTIVE(ii);
                }
                else if (!DAQP_IS_LOWER(ii) && bupper[ii] >= DAQP_INF)
                {
                    DAQP_SET_INACTIVE(ii);
                }
            }
        }
    }

    int daqp_status = daqp_check_bounds(work, bupper, blower);
    if (daqp_status < 0)
        return daqp_status;
    do_activate |= daqp_status;

    if (update_matrices)
    {
        // Normalized LDP rows: packed rows of Rinv (row ii is nonzero from
        // column ii) for the simple bounds and M for the general constraints
        for (int ii = 0, disp = 0; ii < work->ms; disp += nv-ii, ii++)
        {
            c_float *row = work->Rinv + disp;
            blasfeo_unpack_dmat(nv-ii, 1, mem->M_factor, ii, ii, row, nv-ii);
            double norm_squared = 0.0;
            for (int jj = 0; jj < nv-ii; jj++)
                norm_squared += row[jj] * row[jj];
            work->scaling[ii] = 1.0 / sqrt(norm_squared);
            for (int jj = 0; jj < nv-ii; jj++)
                row[jj] *= work->scaling[ii];
        }
        if (work->m > work->ms)
            blasfeo_unpack_dmat(nv, work->m-work->ms, mem->M_factor, 0, work->ms, work->M, nv);
        for (int ii = work->ms; ii < work->m; ii++)
        {
            c_float *row = work->M + (ii-work->ms)*nv;
            double norm_squared = 0.0;
            for (int jj = 0; jj < nv; jj++)
                norm_squared += row[jj] * row[jj];
            if (norm_squared < work->settings->zero_tol)
            {
                work->scaling[ii] = 1.0;
#ifndef DAQP_ASSUME_VALID
                if ((bupper[ii] < -work->settings->zero_tol ||
                        blower[ii] > work->settings->zero_tol) &&
                        !DAQP_IS_IMMUTABLE(ii) && !DAQP_IS_SOFT(ii))
                    return DAQP_EXIT_INFEASIBLE;
#endif
                work->sense[ii] = DAQP_IMMUTABLE;
                continue;
            }
            work->scaling[ii] = 1.0 / sqrt(norm_squared);
            for (int jj = 0; jj < nv; jj++)
                row[jj] *= work->scaling[ii];
        }

        do_activate = 1;
    }

    // Transform the gradient with a triangular solve: v = inv(L)*g.
    blasfeo_dveccp(nv, qp_in->gz, 0, mem->rhs_factor, 0);
    blasfeo_dtrsv_lnn(nv, mem->H_factor, 0, 0,
            mem->rhs_factor, 0, mem->v_factor, 0);
    blasfeo_unpack_dvec(nv, mem->v_factor, 0, work->v, 1);

    // Normalized constraint offsets: d = S*(b + G*v), with G*v from the unnormalized rows
    if (work->m > 0)
        blasfeo_dgemv_t(nv, work->m, 1.0, mem->M_factor, 0, 0,
                mem->v_factor, 0, 0.0, mem->constraint_value, 0,
                mem->constraint_value, 0);
    for (int ii = 0; ii < work->m; ii++)
    {
        double offset = BLASFEO_DVECEL(mem->constraint_value, ii);
        work->dupper[ii] = (bupper[ii] + offset) * work->scaling[ii];
        work->dlower[ii] = (blower[ii] + offset) * work->scaling[ii];
    }
    // DAQP applies constraint scaling to the original-coordinate soft weights.

    if (do_activate)
    {
        reset_daqp_workspace(work);
        daqp_status = daqp_activate_constraints(work);
        if (daqp_status < 0)
            return daqp_status;
    }
    else
    {
        // The transformed RHS changed, so intermediate substitutions in the
        // active-set factorization cannot be reused.
        work->reuse_ind = 0;
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
}



static void dense_qp_daqp_fill_output(dense_qp_daqp_memory *mem, const dense_qp_out *qp_out, const dense_qp_in *qp_in)
{
    int *idxs = mem->idxs;
    int *idxb = qp_in->idxb;
    int i;
    int nv = qp_in->dim->nv;
    int nb = qp_in->dim->nb;
    int ng = qp_in->dim->ng;
    int ns = qp_in->dim->ns;
    DAQPWorkspace *work = mem->daqp_work;

    struct blasfeo_dvec *v = qp_out->v;
    struct blasfeo_dvec *lambda = qp_out->lam;

    // primal variables
    blasfeo_pack_dvec(nv, work->x, 1, v, 0);

    // dual variables
    blasfeo_dvecse(2 * nb + 2 * ng + 2 * ns, 0.0, lambda, 0);
    c_float lam;
    for (i = 0; i < work->n_active; i++)
    {
        lam = work->lam_star[i] * work->scaling[work->WS[i]];
        if (work->WS[i] < nb) // bound constraint
        {
            if (lam >= 0.0)
                BLASFEO_DVECEL(lambda, nb+ng+work->WS[i]) = lam;
            else
                BLASFEO_DVECEL(lambda, work->WS[i]) = -lam;
        }
        else if (work->WS[i] < nb+ng)// general constraint
        {
            if (lam >= 0.0)
                BLASFEO_DVECEL(lambda, nb+ng+work->WS[i]) = lam;
            else
                BLASFEO_DVECEL(lambda, work->WS[i]) = -lam;
        }
        else // equality constraint
            BLASFEO_DVECEL(qp_out->pi, work->WS[i]-nb-ng) = -lam;
    }

    // soft slacks
    int idxdaqp;
    if (ns > 0 && ng > 0)
        blasfeo_dgemv_t(nv, ng, 1.0, qp_in->Ct, 0, 0, v, 0, 0.0,
                mem->constraint_value, nb, mem->constraint_value, nb);
    for (i = 0; i < ns; i++)
    {
        idxdaqp = idxs[i];

        c_float constraint_value = idxdaqp < nb ? BLASFEO_DVECEL(v, idxb[idxdaqp]) :
            BLASFEO_DVECEL(mem->constraint_value, idxdaqp);

        // Recover slacks from primal feasibility. This also handles soft
        // zero rows that DAQP can safely omit from its active-set system.
        // Use the shifted, masked bounds so disabled constraint sides do not
        // contribute a violation to the recovered slack.
        double sl = MAX(mem->d_ls[i], -mem->zl[i]/mem->Zl[i]);
        double su = MAX(mem->d_us[i], -mem->zu[i]/mem->Zu[i]);
        BLASFEO_DVECEL(v, nv+i) = sl + MAX(0, mem->blower[idxdaqp] - constraint_value);
        if (DAQP_IS_IMMUTABLE(idxdaqp) && mem->blower[idxdaqp] > constraint_value)
            BLASFEO_DVECEL(lambda, idxdaqp) =
                mem->Zl[i] * BLASFEO_DVECEL(v, nv+i) + mem->zl[i];
        BLASFEO_DVECEL(lambda, 2*(nb+ng)+i) =
            mem->Zl[i] * BLASFEO_DVECEL(v, nv+i) + mem->zl[i]
            - BLASFEO_DVECEL(lambda, idxs[i]);

        BLASFEO_DVECEL(v, nv+ns+i) = su + MAX(0, constraint_value - mem->bupper[idxdaqp]);
        if (DAQP_IS_IMMUTABLE(idxdaqp) && constraint_value > mem->bupper[idxdaqp])
            BLASFEO_DVECEL(lambda, nb+ng+idxdaqp) =
                mem->Zu[i] * BLASFEO_DVECEL(v, nv+ns+i) + mem->zu[i];
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
