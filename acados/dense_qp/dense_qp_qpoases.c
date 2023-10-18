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
// blasfeo
#include "blasfeo/include/blasfeo_d_aux.h"
#include "blasfeo/include/blasfeo_d_blas.h"

/* Ignore compiler warnings from qpOASES */
#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wtautological-pointer-compare"
#pragma clang diagnostic ignored "-Wunused-parameter"
#pragma clang diagnostic ignored "-Wunused-function"
#include "qpOASES_e/QProblem.h"
#include "qpOASES_e/QProblemB.h"
#pragma clang diagnostic pop
#elif defined(__GNUC__)
#if __GNUC__ >= 6
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-but-set-parameter"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wunused-function"
#include "qpOASES_e/QProblem.h"
#include "qpOASES_e/QProblemB.h"
#pragma GCC diagnostic pop
#else
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wunused-function"
#include "qpOASES_e/QProblem.h"
#include "qpOASES_e/QProblemB.h"
#endif
#else
#include "qpOASES_e/QProblem.h"
#include "qpOASES_e/QProblemB.h"
#endif

// acados
#include "acados/dense_qp/dense_qp_common.h"
#include "acados/dense_qp/dense_qp_qpoases.h"
#include "acados/utils/mem.h"
#include "acados/utils/timing.h"
#include "acados/utils/print.h"
#include "acados/utils/math.h"

#include "acados_c/dense_qp_interface.h"


/************************************************
 * opts
 ************************************************/

acados_size_t dense_qp_qpoases_opts_calculate_size(void *config_, dense_qp_dims *dims)
{
    acados_size_t size = 0;
    size += sizeof(dense_qp_qpoases_opts);

    make_int_multiple_of(8, &size);

    return size;
}



void *dense_qp_qpoases_opts_assign(void *config_, dense_qp_dims *dims, void *raw_memory)
{
    dense_qp_qpoases_opts *opts;

    char *c_ptr = (char *) raw_memory;

    opts = (dense_qp_qpoases_opts *) c_ptr;
    c_ptr += sizeof(dense_qp_qpoases_opts);

    assert((char *) raw_memory + dense_qp_qpoases_opts_calculate_size(config_, dims) >= c_ptr);

    return (void *) opts;
}



void dense_qp_qpoases_opts_initialize_default(void *config_, dense_qp_dims *dims, void *opts_)
{
    dense_qp_qpoases_opts *opts = (dense_qp_qpoases_opts *) opts_;

    opts->max_cputime = 1000.0;
    opts->warm_start = 1;
    opts->max_nwsr = 1000;
    opts->use_precomputed_cholesky = 0;
    opts->hotstart = 0;
    opts->set_acado_opts = 1;
    opts->compute_t = 1;
    opts->tolerance = 1e-4;

    return;
}



void dense_qp_qpoases_opts_update(void *config_, dense_qp_dims *dims, void *opts_)
{
    //    dense_qp_qpoases_opts *opts = (dense_qp_qpoases_opts *)opts_;

    return;
}



void dense_qp_qpoases_opts_set(void *config_, void *opts_, const char *field, void *value)
{
    dense_qp_qpoases_opts *opts = opts_;

    if (!strcmp(field, "tol_stat"))
    {
        double *tol = value;
        opts->tolerance = MIN(opts->tolerance, *tol);
    }
    else if (!strcmp(field, "tol_eq"))
    {
        double *tol = value;
        opts->tolerance = MIN(opts->tolerance, *tol);
    }
    else if (!strcmp(field, "tol_ineq"))
    {
        double *tol = value;
        opts->tolerance = MIN(opts->tolerance, *tol);
    }
    else if (!strcmp(field, "tol_comp"))
    {
        double *tol = value;
        opts->tolerance = MIN(opts->tolerance, *tol);
    }
    else if (!strcmp(field, "warm_start"))
    {
        int *warm_start = value;
        if (*warm_start == 0)
        {
            opts->warm_start = 0;
            opts->hotstart = 0;
        }
        else if (*warm_start == 1)
        {
            opts->warm_start = 1;
            opts->hotstart = 0;
        }
        else if (*warm_start == 2)
        {
            opts->warm_start = 1;
            opts->hotstart = 1;
        }
        else
        {
            printf("\ndense_qp_qpoases: setting warm_start: supported values are: 0 - cold, 1 - warm, 2 - hot\n");
            exit(1);
        }
    }
    else if (!strcmp(field, "iter_max"))
    {
        int *max_iter = value;
        opts->max_nwsr = *max_iter;
    }
    else
    {
        printf("\nerror: dense_qp_qpoases_opts_set: wrong field: %s\n", field);
        exit(1);
    }

    return;
}



/************************************************
 * memory
 ************************************************/

acados_size_t dense_qp_qpoases_memory_calculate_size(void *config_, dense_qp_dims *dims, void *opts_)
{
    dense_qp_dims dims_stacked;

    int nv  = dims->nv;
    int ne  = dims->ne;
    int ng  = dims->ng;
    int nb  = dims->nb;
    int nsb = dims->nsb;
    // int nsg = dims->nsg;
    int ns  = dims->ns;

    int nv2 = nv + 2*ns;
    int ng2 = (ns > 0) ? ng + nsb : ng;
    int nb2 = nb - nsb + 2 * ns;

    // size in bytes
    acados_size_t size = sizeof(dense_qp_qpoases_memory);

    size += 1 * nv * nv * sizeof(double);      // H
    size += 1 * nv2 * nv2 * sizeof(double);    // HH
    size += 1 * nv2 * nv2 * sizeof(double);    // R
    size += 1 * nv2 * ne * sizeof(double);     // A
    size += 1 * nv * ng * sizeof(double);      // C
    size += 1 * nv2 * ng2 * sizeof(double);    // CC
    size += 1 * nv * sizeof(double);           // g
    size += 1 * nv2 * sizeof(double);          // gg
    size += 2 * nv2 * sizeof(double);          // d_lb d_ub
    size += 1 * ne * sizeof(double);           // b
    size += 2 * nb2 * sizeof(double);          // d_lb0 d_ub0
    size += 2 * ng * sizeof(double);           // d_lg0 d_ug0
    size += 2 * ng2 * sizeof(double);          // d_lg d_ug
    size += 1 * nb * sizeof(int);              // idxb
    size += 1 * nb2 * sizeof(int);             // idxb_stacked
    size += 1 * ns * sizeof(int);              // idxs
    size += 1 * nv2 * sizeof(double);          // prim_sol
    size += 1 * (nv2 + ng2) * sizeof(double);  // dual_sol
    size += 6 * ns * sizeof(double);           // Zl, Zu, zl, zu, d_ls, d_us

    if (ns > 0)
    {
        dense_qp_stack_slacks_dims(dims, &dims_stacked);
        size += dense_qp_in_calculate_size(&dims_stacked);
    }

    if (ng > 0 || ns > 0)  // QProblem
        size += QProblem_calculateMemorySize(nv2, ng2);
    else  // QProblemB
        size += QProblemB_calculateMemorySize(nv);

    make_int_multiple_of(8, &size);

    return size;
}



void *dense_qp_qpoases_memory_assign(void *config_, dense_qp_dims *dims, void *opts_,
                                     void *raw_memory)
{
    dense_qp_qpoases_memory *mem;
    dense_qp_dims dims_stacked;

    int nv  = dims->nv;
    int ne  = dims->ne;
    int ng  = dims->ng;
    int nb  = dims->nb;
    int nsb = dims->nsb;
    // int nsg = dims->nsg;
    int ns  = dims->ns;

    int nv2 = nv + 2*ns;
    int ng2 = (ns > 0) ? ng + nsb : ng;
    int nb2 = nb - nsb + 2 * ns;

    // char pointer
    char *c_ptr = (char *) raw_memory;

    mem = (dense_qp_qpoases_memory *) c_ptr;
    c_ptr += sizeof(dense_qp_qpoases_memory);

    assert((size_t) c_ptr % 8 == 0 && "memory not 8-byte aligned!");

    if (ns > 0)
    {
        dense_qp_stack_slacks_dims(dims, &dims_stacked);
        mem->qp_stacked = dense_qp_in_assign(&dims_stacked, c_ptr);
        c_ptr += dense_qp_in_calculate_size(&dims_stacked);
    }
    else
    {
        mem->qp_stacked = NULL;
    }

    assert((size_t) c_ptr % 8 == 0 && "memory not 8-byte aligned!");

    assign_and_advance_double(nv * nv, &mem->H, &c_ptr);
    assign_and_advance_double(nv2 * nv2, &mem->HH, &c_ptr);
    assign_and_advance_double(nv2 * nv2, &mem->R, &c_ptr);
    assign_and_advance_double(nv2 * ne, &mem->A, &c_ptr);
    assign_and_advance_double(nv * ng, &mem->C, &c_ptr);
    assign_and_advance_double(nv2 * ng2, &mem->CC, &c_ptr);
    assign_and_advance_double(nv, &mem->g, &c_ptr);
    assign_and_advance_double(nv2, &mem->gg, &c_ptr);
    assign_and_advance_double(ne, &mem->b, &c_ptr);
    assign_and_advance_double(nb2, &mem->d_lb0, &c_ptr);
    assign_and_advance_double(nb2, &mem->d_ub0, &c_ptr);
    assign_and_advance_double(nv2, &mem->d_lb, &c_ptr);
    assign_and_advance_double(nv2, &mem->d_ub, &c_ptr);
    assign_and_advance_double(ng, &mem->d_lg0, &c_ptr);
    assign_and_advance_double(ng, &mem->d_ug0, &c_ptr);
    assign_and_advance_double(ng2, &mem->d_lg, &c_ptr);
    assign_and_advance_double(ng2, &mem->d_ug, &c_ptr);
    assign_and_advance_double(ns, &mem->d_ls, &c_ptr);
    assign_and_advance_double(ns, &mem->d_us, &c_ptr);
    assign_and_advance_double(ns, &mem->Zl, &c_ptr);
    assign_and_advance_double(ns, &mem->Zu, &c_ptr);
    assign_and_advance_double(ns, &mem->zl, &c_ptr);
    assign_and_advance_double(ns, &mem->zu, &c_ptr);
    assign_and_advance_double(nv2, &mem->prim_sol, &c_ptr);
    assign_and_advance_double(nv2 + ng2, &mem->dual_sol, &c_ptr);

    // TODO(dimitris): update assign syntax in qpOASES
    assert((size_t) c_ptr % 8 == 0 && "double not 8-byte aligned!");

    if (ng > 0 || ns > 0)
    {  // QProblem
        QProblem_assignMemory(nv2, ng2, (QProblem **) &(mem->QP), c_ptr);
        c_ptr += QProblem_calculateMemorySize(nv2, ng2);
    }
    else
    {  // QProblemB
        QProblemB_assignMemory(nv, (QProblemB **) &(mem->QPB), c_ptr);
        c_ptr += QProblemB_calculateMemorySize(nv);
    }

    assign_and_advance_int(nb, &mem->idxb, &c_ptr);
    assign_and_advance_int(nb2, &mem->idxb_stacked, &c_ptr);
    assign_and_advance_int(ns, &mem->idxs, &c_ptr);

    assert((char *) raw_memory + dense_qp_qpoases_memory_calculate_size(config_, dims, opts_) >=
           c_ptr);

    // assign default values to fields stored in the memory
    mem->first_it = 1;  // only used if hotstart (only constant data matrices) is enabled

    return mem;
}



void dense_qp_qpoases_memory_get(void *config_, void *mem_, const char *field, void* value)
{
    // qp_solver_config *config = config_;
    dense_qp_qpoases_memory *mem = mem_;

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
        printf("\nerror: dense_qp_qpoases_memory_get: field %s not available\n", field);
        exit(1);
    }

    return;

}



/************************************************
 * workspace
 ************************************************/

acados_size_t dense_qp_qpoases_workspace_calculate_size(void *config_, dense_qp_dims *dims, void *opts_)
{
    return 0;
}



/************************************************
 * functions
 ************************************************/

int dense_qp_qpoases(void *config_, dense_qp_in *qp_in, dense_qp_out *qp_out, void *opts_,
                     void *memory_, void *work_)
{
    qp_info *info = (qp_info *) qp_out->misc;
    acados_timer tot_timer, qp_timer, interface_timer;

    acados_tic(&tot_timer);
    acados_tic(&interface_timer);
    info->t_computed = 0;

    // cast structures
    dense_qp_qpoases_opts *opts = (dense_qp_qpoases_opts *) opts_;
    dense_qp_qpoases_memory *memory = (dense_qp_qpoases_memory *) memory_;

    // extract qpoases data
    double *H = memory->H;
    double *HH = memory->HH;
    double *A = memory->A;
    double *C = memory->C;
    double *CC = memory->CC;
    double *g = memory->g;
    double *gg = memory->gg;
    double *b = memory->b;
    double *d_lb0 = memory->d_lb0;
    double *d_ub0 = memory->d_ub0;
    double *d_lb = memory->d_lb;
    double *d_ub = memory->d_ub;
    double *d_lg0 = memory->d_lg0;
    double *d_ug0 = memory->d_ug0;
    double *d_lg = memory->d_lg;
    double *d_ug = memory->d_ug;
    double *Zl = memory->Zl;
    double *Zu = memory->Zu;
    double *zl = memory->zl;
    double *zu = memory->zu;
    double *d_ls = memory->d_ls;
    double *d_us = memory->d_us;
    int *idxb = memory->idxb;
    int *idxb_stacked = memory->idxb_stacked;
    int *idxs = memory->idxs;
    double *prim_sol = memory->prim_sol;
    double *dual_sol = memory->dual_sol;
    QProblemB *QPB = memory->QPB;
    QProblem *QP = memory->QP;
    dense_qp_in *qp_stacked = memory->qp_stacked;

    // extract dense qp size
    int nv  = qp_in->dim->nv;
    // int ne  = qp_in->dim->ne;
    int ng  = qp_in->dim->ng;
    int nb  = qp_in->dim->nb;
    int nsb = qp_in->dim->nsb;
    // int nsg = qp_in->dim->nsg;
    int ns  = qp_in->dim->ns;

    int nv2 = nv + 2*ns;
    int ng2 = (ns > 0) ? ng + nsb : ng;
    int nb2 = nb - nsb + 2 * ns;

    // fill in the upper triangular of H in dense_qp
    blasfeo_dtrtr_l(nv, qp_in->Hv, 0, 0, qp_in->Hv, 0, 0);

    // extract data from qp_in in row-major
    d_dense_qp_get_all_rowmaj(qp_in, H, g, A, b, idxb, d_lb0, d_ub0, C, d_lg0, d_ug0,
                                 Zl, Zu, zl, zu, idxs, d_ls, d_us);

    // reorder box constraints bounds
    for (int ii = 0; ii < nv2; ii++)
    {
        d_lb[ii] = -QPOASES_INFTY;
        d_ub[ii] = +QPOASES_INFTY;
    }

    if (ns > 0)
    {
        dense_qp_stack_slacks(qp_in, qp_stacked);
        d_dense_qp_get_all_rowmaj(qp_stacked, HH, gg, A, b, idxb_stacked, d_lb0, d_ub0, CC, d_lg,
            d_ug, NULL, NULL, NULL, NULL, NULL, NULL, NULL);

        for (int ii = 0; ii < nb2; ii++)
        {
            d_lb[idxb_stacked[ii]] = d_lb0[ii];
            d_ub[idxb_stacked[ii]] = d_ub0[ii];
        }
    }
    else
    {
        for (int ii = 0; ii < nb; ii++)
        {
            d_lb[idxb[ii]] = d_lb0[ii];
            d_ub[idxb[ii]] = d_ub0[ii];
        }
    }

    // cholesky factorization of H
    // blasfeo_dpotrf_l(nvd, qpd->Hv, 0, 0, sR, 0, 0);

    // fill in upper triangular of R
    // blasfeo_dtrtr_l(nvd, sR, 0, 0, sR, 0, 0);

    // extract R
    // blasfeo_unpack_dmat(nvd, nvd, sR, 0, 0, R, nvd);

    info->interface_time = acados_toc(&interface_timer);
    acados_tic(&qp_timer);

    // solve dense qp
    int nwsr = opts->max_nwsr;
    double cputime = opts->max_cputime;

    int qpoases_status = 0;
    if (opts->hotstart == 1)
    {  // only to be used with fixed data matrices!
        if (ng > 0 || ns > 0)
        {  // QProblem
            if (memory->first_it == 1)
            {
                QProblemCON(QP, nv2, ng2, HST_POSDEF);
                QProblem_setPrintLevel(QP, PL_MEDIUM);
                // QProblem_setPrintLevel(QP, PL_DEBUG_ITER);
                if (opts->set_acado_opts)
                {
                    static Options options;
                    Options_setToMPC(&options);
                    options.terminationTolerance = opts->tolerance;
                    QProblem_setOptions(QP, options);
                }

                qpoases_status = (ns > 0) ?
                    QProblem_init(QP, HH, gg, CC, d_lb, d_ub, d_lg, d_ug, &nwsr, &cputime) :
                    QProblem_init(QP, H, g, C, d_lb, d_ub, d_lg0, d_ug0, &nwsr, &cputime);
                memory->first_it = 0;

                QProblem_getPrimalSolution(QP, prim_sol);
                QProblem_getDualSolution(QP, dual_sol);
            }
            else
            {
                qpoases_status = (ns > 0) ?
                    QProblem_hotstart(QP, gg, d_lb, d_ub, d_lg, d_ug, &nwsr, &cputime) :
                    QProblem_hotstart(QP, g, d_lb, d_ub, d_lg0, d_ug0, &nwsr, &cputime);

                QProblem_getPrimalSolution(QP, prim_sol);
                QProblem_getDualSolution(QP, dual_sol);
            }
        }
        else
        {
            if (memory->first_it == 1)
            {
                QProblemBCON(QPB, nv, HST_POSDEF);
                QProblemB_setPrintLevel(QPB, PL_MEDIUM);
                // QProblemB_setPrintLevel(QPB, PL_DEBUG_ITER);
                if (opts->set_acado_opts)
                {
                    static Options options;
                    Options_setToMPC(&options);
                    options.terminationTolerance = opts->tolerance;
                    QProblem_setOptions(QP, options);
                }
                QProblemB_init(QPB, H, g, d_lb, d_ub, &nwsr, &cputime);
                memory->first_it = 0;

                QProblemB_getPrimalSolution(QPB, prim_sol);
                QProblemB_getDualSolution(QPB, dual_sol);
            }
            else
            {
                QProblemB_hotstart(QPB, g, d_lb, d_ub, &nwsr, &cputime);

                QProblemB_getPrimalSolution(QPB, prim_sol);
                QProblemB_getDualSolution(QPB, dual_sol);
            }
        }
    }
    else
    {  // hotstart = 0
        if (ng > 0 || ns > 0)
        {
            QProblemCON(QP, nv2, ng2, HST_POSDEF);
            // QProblem_setPrintLevel(QP, PL_HIGH);
            QProblem_setPrintLevel(QP, PL_DEBUG_ITER);
            QProblem_printProperties(QP);
            if (opts->use_precomputed_cholesky == 1)
            {
                // NOTE(oj): why are there no options set in this case?
                // static Options options;
                // Options_setToDefault( &options );
                // options.initialStatusBounds = ST_INACTIVE;
                // options.terminationTolerance = opts->tolerance;
                // QProblem_setOptions( QP, options );


                qpoases_status = (ns > 0) ?
                    QProblem_initW(QP, HH, gg, CC, d_lb, d_ub, d_lg, d_ug, &nwsr, &cputime,
                                   /* primal_sol */ NULL, /* dual sol */ NULL,
                                   /* guessed bounds */ NULL, /* guessed constraints */ NULL,
                                   /* R */ memory->R) :
                    QProblem_initW(QP, H, g, C, d_lb, d_ub, d_lg0, d_ug0, &nwsr, &cputime,
                                   /* primal_sol */ NULL, /* dual sol */ NULL,
                                   /* guessed bounds */ NULL, /* guessed constraints */ NULL,
                                   /* R */ memory->R);  // NOTE(dimitris): can pass either NULL or 0
            }
            else
            {
                if (opts->set_acado_opts)
                {
                    static Options options;
                    Options_setToMPC(&options);
                    options.terminationTolerance = opts->tolerance;
                    QProblem_setOptions(QP, options);
                }
                if (opts->warm_start)
                {
                    qpoases_status = (ns > 0) ?
                        QProblem_initW(QP, HH, gg, CC, d_lb, d_ub, d_lg, d_ug, &nwsr, &cputime,
                                      /* primal_sol */ NULL, dual_sol, /* guessed bounds */ NULL,
                                      /* guessed constraints */ NULL, /* R */ NULL) :
                        QProblem_initW(QP, H, g, C, d_lb, d_ub, d_lg0, d_ug0, &nwsr,
                                       &cputime, NULL, dual_sol, NULL, NULL, NULL);
                }
                else
                {
                    qpoases_status = (ns > 0) ?
                        QProblem_init(QP, HH, gg, CC, d_lb, d_ub, d_lg, d_ug, &nwsr, &cputime) :
                        QProblem_init(QP, H, g, C, d_lb, d_ub, d_lg0, d_ug0, &nwsr, &cputime);
                }
            }
            QProblem_getPrimalSolution(QP, prim_sol);
            QProblem_getDualSolution(QP, dual_sol);
        }
        else
        {  // QProblemB
            QProblemBCON(QPB, nv, HST_POSDEF);
            // QProblemB_setPrintLevel(QPB, PL_MEDIUM);
            QProblemB_setPrintLevel(QPB, PL_DEBUG_ITER);
            QProblemB_printProperties(QPB);
            if (opts->use_precomputed_cholesky == 1)
            {
                // NOTE(oj): why are there no options set in this case?
                // static Options options;
                // Options_setToDefault( &options );
                // options.initialStatusBounds = ST_INACTIVE;
                // QProblemB_setOptions( QPB, options );
                qpoases_status = QProblemB_initW(QPB, H, g, d_lb, d_ub, &nwsr, &cputime,
                                                 /* primal_sol */ NULL, /* dual sol */ NULL,
                                                 /* guessed bounds */ NULL,
                                                 /* R */ memory->R);
            }
            else
            {
                if (opts->set_acado_opts)
                {
                    static Options options;
                    Options_setToMPC(&options);
                    options.terminationTolerance = opts->tolerance;
                    QProblemB_setOptions(QPB, options);
                }
                if (opts->warm_start)
                {
                    qpoases_status = QProblemB_initW(QPB, H, g, d_lb, d_ub, &nwsr, &cputime,
                                                     /* primal sol */ NULL, /* dual sol */ dual_sol,
                                                     /* guessed bounds */ NULL,
                                                     /* R */ NULL);
                }
                else
                {
                    qpoases_status = QProblemB_init(QPB, H, g, d_lb, d_ub, &nwsr, &cputime);
                }
            }
            QProblemB_getPrimalSolution(QPB, prim_sol);
            QProblemB_getDualSolution(QPB, dual_sol);
        }
    }

    // save solution statistics to memory
    memory->cputime = cputime;
    memory->nwsr = nwsr;
    info->solve_QP_time = acados_toc(&qp_timer);

    acados_tic(&interface_timer);
    // copy prim_sol and dual_sol to qp_out
    blasfeo_pack_dvec(nv2, prim_sol, 1, qp_out->v, 0);
    blasfeo_dvecse(2 * nb + 2 * ng + 2 * ns, 0.0, qp_out->lam, 0);
    for (int ii = 0; ii < nb; ii++)
    {
        if (dual_sol[idxb[ii]] >= 0.0)
            qp_out->lam->pa[ii] = dual_sol[idxb[ii]];
        else
            qp_out->lam->pa[nb + ng + ii] = -dual_sol[idxb[ii]];
    }

    for (int ii = 0; ii < ng; ii++)
    {
        if (dual_sol[nv2 + ii] >= 0.0)
            qp_out->lam->pa[nb + ii] = dual_sol[nv2 + ii];
        else
            qp_out->lam->pa[2 * nb + ng + ii] = -dual_sol[nv2 + ii];
    }

    int k = 0;
    for (int ii = 0; ii < ns; ii++)
    {
        int js = idxs[ii];

        double offset_l = 0.0;
        double offset_u = 0.0;

        if (js < nb)
        {
            if (dual_sol[nv2 + ng + k] <= 0.0)  // softened upper box constraints
            {
                qp_out->lam->pa[nb + ng + js] = -dual_sol[nv2 + ng + k];
                offset_u = -dual_sol[nv2 + ng + k];
            }
            else  // softened lower box constraints
            {
                qp_out->lam->pa[js] = dual_sol[nv2 + ng + k];
                offset_l = dual_sol[nv2 + ng + k];
            }

            k++;
        }
        else
        {
            offset_l = qp_out->lam->pa[nb+js-nb];
            offset_u = qp_out->lam->pa[2*nb+ng+js-nb];
        }

        // dual variables for sl >= d_ls
        if (dual_sol[nv + ii] >= 0)
            qp_out->lam->pa[2*nb + 2*ng + ii] = dual_sol[nv + ii] - offset_u;

        // dual variables for su >= d_us
        if (dual_sol[nv + ns + ii] >= 0)
            qp_out->lam->pa[2*nb + 2*ng + ns + ii] = dual_sol[nv + ns + ii] - offset_l;
    }

    info->interface_time += acados_toc(&interface_timer);
    info->total_time = acados_toc(&tot_timer);
    info->num_iter = nwsr;

    memory->time_qp_solver_call = info->solve_QP_time;
    memory->iter = nwsr;

    // compute slacks
    if (opts->compute_t)
    {
        dense_qp_compute_t(qp_in, qp_out);
        info->t_computed = 1;
    }

    int acados_status = qpoases_status;
    if (qpoases_status == SUCCESSFUL_RETURN) acados_status = ACADOS_SUCCESS;
    if (qpoases_status == RET_MAX_NWSR_REACHED) acados_status = ACADOS_MAXITER;
    return acados_status;
}



void dense_qp_qpoases_eval_sens(void *config_, void *qp_in, void *qp_out, void *opts_, void *mem_, void *work_)
{
    printf("\nerror: dense_qp_qpoases_eval_sens: not implemented yet\n");
    exit(1);
}


// void dense_qp_qpoases_memory_reset(void *config_, dense_qp_in *qp_in, dense_qp_out *qp_out, void *opts_,
//                      void *memory_, void *work_)
// {
//     printf("\nerror: dense_qp_qpoases_memory_reset: not implemented yet\n");
//     // exit(1);
// }



void dense_qp_qpoases_memory_reset(void *config, void *qp_in, void *qp_out, void *opts, void *mem, void *work)
{
    printf("\nerror: dense_qp_qpoases_memory_reset: not implemented yet\n");
    exit(1);
}

void dense_qp_qpoases_solver_get(void *config_, void *qp_in_, void *qp_out_, void *opts_, void *mem_, const char *field, int stage, void* value, int size1, int size2)
{
    printf("\nerror: dense_qp_qpoases_solver_get: not implemented yet\n");
    exit(1);
}


void dense_qp_qpoases_config_initialize_default(void *config_)
{
    qp_solver_config *config = config_;

    config->opts_calculate_size = (acados_size_t (*)(void *, void *)) & dense_qp_qpoases_opts_calculate_size;
    config->opts_assign = (void *(*) (void *, void *, void *) ) & dense_qp_qpoases_opts_assign;
    config->opts_initialize_default =
        (void (*)(void *, void *, void *)) & dense_qp_qpoases_opts_initialize_default;
    config->opts_update = (void (*)(void *, void *, void *)) & dense_qp_qpoases_opts_update;
    config->opts_set = &dense_qp_qpoases_opts_set;
    config->memory_calculate_size =
        (acados_size_t (*)(void *, void *, void *)) & dense_qp_qpoases_memory_calculate_size;
    config->memory_assign =
        (void *(*) (void *, void *, void *, void *) ) & dense_qp_qpoases_memory_assign;
    config->memory_get = &dense_qp_qpoases_memory_get;
    config->workspace_calculate_size =
        (acados_size_t (*)(void *, void *, void *)) & dense_qp_qpoases_workspace_calculate_size;
    config->eval_sens = &dense_qp_qpoases_eval_sens;
    // config->memory_reset = &dense_qp_qpoases_memory_reset;
    config->evaluate = (int (*)(void *, void *, void *, void *, void *, void *)) & dense_qp_qpoases;
    config->memory_reset = &dense_qp_qpoases_memory_reset;
    config->solver_get = &dense_qp_qpoases_solver_get;

    return;
}
