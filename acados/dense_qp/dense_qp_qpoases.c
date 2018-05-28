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

/************************************************
 * opts
 ************************************************/

int dense_qp_qpoases_opts_calculate_size(void *config_, dense_qp_dims *dims)
{
    int size = 0;
    size += sizeof(dense_qp_qpoases_opts);

    return size;
}

void *dense_qp_qpoases_opts_assign(void *config_, dense_qp_dims *dims, void *raw_memory)
{
    dense_qp_qpoases_opts *opts;

    char *c_ptr = (char *) raw_memory;

    opts = (dense_qp_qpoases_opts *) c_ptr;
    c_ptr += sizeof(dense_qp_qpoases_opts);

    assert((char *) raw_memory + dense_qp_qpoases_opts_calculate_size(config_, dims) == c_ptr);

    return (void *) opts;
}

void dense_qp_qpoases_opts_initialize_default(void *config_, dense_qp_dims *dims, void *opts_)
{
    dense_qp_qpoases_opts *opts = (dense_qp_qpoases_opts *) opts_;

    opts->max_cputime = 1000.0;
    opts->warm_start = 0;
    opts->max_nwsr = 1000;
    opts->use_precomputed_cholesky = 0;
    opts->hotstart = 0;
    opts->set_acado_opts = 1;
    opts->compute_t = 1;

    return;
}

void dense_qp_qpoases_opts_update(void *config_, dense_qp_dims *dims, void *opts_)
{
    //    dense_qp_qpoases_opts *opts = (dense_qp_qpoases_opts *)opts_;

    return;
}

/************************************************
 * memory
 ************************************************/

int dense_qp_qpoases_memory_calculate_size(void *config_, dense_qp_dims *dims, void *opts_)
{
    int nv  = dims->nv;
    int ne  = dims->ne;
    int ng  = dims->ng;
    int nb  = dims->nb;
    int nsb = dims->nsb;
    int nsg = dims->nsg;
    int ns  = dims->ns;

    int nv2 = nv + 2*ns;
    int ng2 = (ns > 0) ? 2*(ng + nsb) : ng;

    // size in bytes
    int size = sizeof(dense_qp_qpoases_memory);

    size += 1 * nv * nv * sizeof(double);      // H
    size += 1 * nv2 * nv2 * sizeof(double);    // HH
    size += 1 * nv2 * nv2 * sizeof(double);    // R
    size += 1 * nv * ne * sizeof(double);      // A
    size += 1 * nv * ng * sizeof(double);      // C
    size += 1 * nv2 * ng2 * sizeof(double);    // CC
    size += 1 * nv * sizeof(double);           // g
    size += 1 * nv2 * sizeof(double);          // gg
    size += 2 * nv2 * sizeof(double);          // d_lb d_ub
    size += 1 * ne * sizeof(double);           // b
    size += 2 * nb * sizeof(double);           // d_lb0 d_ub0
    size += 2 * ng * sizeof(double);           // d_lg0 d_ug0
    size += 2 * ng2 * sizeof(double);          // d_lg d_ug
    size += 1 * nb * sizeof(int);              // idxb
    size += 1 * ns * sizeof(int);              // idxs
    size += 1 * nv2 * sizeof(double);          // prim_sol
    size += 1 * (nv2 + ng2) * sizeof(double);  // dual_sol
    size += 6 * ns * sizeof(double);           // Zl, Zu, zl, zu, d_ls, d_us

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

    int nv  = dims->nv;
    int ne  = dims->ne;
    int ng  = dims->ng;
    int nb  = dims->nb;
    int nsb = dims->nsb;
    int nsg = dims->nsg;
    int ns  = dims->ns;

    int nv2 = nv + 2*ns;
    int ng2 = (ns > 0) ? 2*(ng + nsb) : ng;

    // char pointer
    char *c_ptr = (char *) raw_memory;

    mem = (dense_qp_qpoases_memory *) c_ptr;
    c_ptr += sizeof(dense_qp_qpoases_memory);

    assert((size_t) c_ptr % 8 == 0 && "double not 8-byte aligned!");

    assign_and_advance_double(nv * nv, &mem->H, &c_ptr);
    assign_and_advance_double(nv2 * nv2, &mem->HH, &c_ptr);
    assign_and_advance_double(nv2 * nv2, &mem->R, &c_ptr);
    assign_and_advance_double(nv * ne, &mem->A, &c_ptr);
    assign_and_advance_double(nv * ng, &mem->C, &c_ptr);
    assign_and_advance_double(nv2 * ng2, &mem->CC, &c_ptr);
    assign_and_advance_double(nv, &mem->g, &c_ptr);
    assign_and_advance_double(nv2, &mem->gg, &c_ptr);
    assign_and_advance_double(ne, &mem->b, &c_ptr);
    assign_and_advance_double(nb, &mem->d_lb0, &c_ptr);
    assign_and_advance_double(nb, &mem->d_ub0, &c_ptr);
    assign_and_advance_double(nv2, &mem->d_lb, &c_ptr);
    assign_and_advance_double(nv2, &mem->d_ub, &c_ptr);
    assign_and_advance_double(ng, &mem->d_lg0, &c_ptr);
    assign_and_advance_double(ng, &mem->d_ug0, &c_ptr);
    assign_and_advance_double(ng2, &mem->d_lg, &c_ptr);
    assign_and_advance_double(ng2, &mem->d_ug, &c_ptr);
    assign_and_advance_double(ns, &mem->Zl, &c_ptr);
    assign_and_advance_double(ns, &mem->Zu, &c_ptr);
    assign_and_advance_double(ns, &mem->zl, &c_ptr);
    assign_and_advance_double(ns, &mem->zu, &c_ptr);
    assign_and_advance_double(ns, &mem->d_ls, &c_ptr);
    assign_and_advance_double(ns, &mem->d_us, &c_ptr);
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
    assign_and_advance_int(ns, &mem->idxs, &c_ptr);

    assert((char *) raw_memory + dense_qp_qpoases_memory_calculate_size(config_, dims, opts_) >=
           c_ptr);

    // assign default values to fields stored in the memory
    mem->first_it = 1;  // only used if hotstart (only constant data matrices) is enabled

    return mem;
}

/************************************************
 * workspcae
 ************************************************/

int dense_qp_qpoases_workspace_calculate_size(void *config_, dense_qp_dims *dims, void *opts_)
{
    return 0;
}

/************************************************
 * functions
 ************************************************/

int dense_qp_qpoases(void *config_, dense_qp_in *qp_in, dense_qp_out *qp_out, void *opts_,
                     void *memory_, void *work_)
{
    dense_qp_info *info = (dense_qp_info *) qp_out->misc;
    acados_timer tot_timer, qp_timer, interface_timer;

    acados_tic(&tot_timer);
    acados_tic(&interface_timer);

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
    double *d_ls = memory->d_ls;
    double *d_us = memory->d_us;
    double *Zl = memory->Zl;
    double *Zu = memory->Zu;
    double *zl = memory->zl;
    double *zu = memory->zu;
    int *idxb = memory->idxb;
    int *idxs = memory->idxs;
    double *prim_sol = memory->prim_sol;
    double *dual_sol = memory->dual_sol;
    QProblemB *QPB = memory->QPB;
    QProblem *QP = memory->QP;

    // extract dense qp size
    int nv  = qp_in->dim->nv;
    int ne  = qp_in->dim->ne;
    int ng  = qp_in->dim->ng;
    int nb  = qp_in->dim->nb;
    int nsb = qp_in->dim->nsb;
    int nsg = qp_in->dim->nsg;
    int ns  = qp_in->dim->ns;

    int nv2 = nv + 2*ns;
    int ng2 = (ns > 0) ? 2*(ng + nsb) : ng;

    // fill in the upper triangular of H in dense_qp
    blasfeo_dtrtr_l(nv, qp_in->Hv, 0, 0, qp_in->Hv, 0, 0);

    // dense qp row-major
    d_cvt_dense_qp_to_rowmaj(qp_in, H, g, A, b, idxb, d_lb0, d_ub0, C, d_lg0, d_ug0,
                             Zl, Zu, zl, zu, idxs, d_ls, d_us);

    // in case ns > 0, slack variables are introduced, i.e. v2 = [v; sl; su]
    // reorder box constraints bounds
    for (int ii = 0; ii < nv2; ii++)
    {
        d_lb[ii] = -QPOASES_INFTY;
        d_ub[ii] = +QPOASES_INFTY;
    }
    for (int ii = 0; ii < nb; ii++)
    {
        d_lb[idxb[ii]] = d_lb0[ii];
        d_ub[idxb[ii]] = d_ub0[ii];
    }

    // if a box constraint is softened, remove it from box constraints bounds
    for (int ii = 0; ii < ns; ii++)
    {
        int js = idxs[ii];

        if (js < nb)
        {
            d_lb[idxb[js]] = -QPOASES_INFTY;
            d_ub[idxb[js]] = +QPOASES_INFTY;
        }
    }

    // copy d_ls and d_us to box constraints bounds
    for (int ii = 0; ii < ns; ii++)
    {
        d_lb[nv+ii] = d_ls[ii];
        d_lb[nv+ns+ii] = d_us[ii];

        d_ub[nv+ii] = +QPOASES_INFTY;
        d_ub[nv+ns+ii] = +QPOASES_INFTY;
    }

    if (ns > 0)
    {
        // in case ns > 0, HH = [H 0 0; 0 Zl 0; 0 0 Zu]

        // initialize HH tp 0.0
        for (int ii = 0; ii < nv2 * nv2; ii++) HH[ii] = 0.0;

        // copy H to upper left corner of HH
        for (int ii = 0; ii < nv; ii++)
        {
            for (int jj = 0; jj < nv; jj++)
            {
                HH[jj + nv2*ii] = H[jj + nv*ii];
            }
        }

        // copy Zl and Zu to the main diagonal of HH
        for (int ii = 0; ii < ns; ii++)
        {
            int ll = ii + nv;
            int uu = ii + nv + ns;

            HH[ll + nv2*ll] = Zl[ii];
            HH[uu + nv2*uu] = Zu[ii];
        }

        // in case ns > 0, gg = [g; zl; zu]
        // copy g to gg
        for (int ii = 0; ii < nv; ii++) gg[ii] = g[ii];

        // copy zl and zu to gg
        for (int ii = 0; ii < ns; ii++)
        {
            gg[ii+nv] = zl[ii];
            gg[ii+nv+ns] = zu[ii];
        }

        for (int ii = 0; ii < ng2; ii++)
        {
            // initialize d_ug to 0.0
            d_ug[ii] = 0.0;

            // set d_lg[i] to -INFTY
            d_lg[ii] = -QPOASES_INFTY;
        }

        // copy [d_ug0;-d_lg0] to d_ug
        for (int ii = 0; ii < ng; ii++)
        {
            d_ug[ii] = d_ug0[ii];
            d_ug[ii+ng] = -d_lg0[ii];
        }

        // in case ns > 0, CC contains softened box constraints and
        // both regular and softened general constraints

        // initialize CC to 0.0
        for (int ii = 0; ii < nv2 * ng2; ii++) CC[ii] = 0.0;

        // copy [C; -C] to upper left corner of CC
        for (int ii = 0; ii < ng; ii++)
        {
            for (int jj = 0; jj < nv; jj++)
            {
                CC[jj + nv2*ii] = C[jj + nv*ii];
                CC[jj + nv2*(ii+ng)] = -C[jj + nv*ii];
            }
        }

        // insert soft constraints into CC
        int k_sb = 0, k_sg = 0, row_b = 2*ng;
        for (int ii = 0; ii < ns; ii++)
        {
            int js = idxs[ii];

            if (js < nb)
            {
                // index of a soft box constraint
                int jx = idxb[js];

                // x_i - su_i <= ub_i
                CC[jx + nv2*row_b] = 1.0;
                CC[nv + ns + k_sb + nv2*row_b] = -1.0;
                d_ug[row_b] = d_ub0[js];


                // -x_i - sl_i <= -lb_i
                CC[jx + nv2*row_b] = -1.0;
                CC[nv + k_sb + nv2*(row_b+nsb)] = -1.0;
                d_ug[row_b+nsb] = -d_lb0[js];

                row_b++; k_sb++;
            }
            else
            {
                // index of a soft general constraint
                int row_g = js - nb;

                // C_i x - su_i <= ug_i
                CC[nv + ns + nsb + k_sg + nv2*row_g] = -1.0;

                row_g += ng;

                // -C_i x - sl_i <= -lg_i
                CC[nv + nsb + k_sg + nv2*row_g] = -1.0;

                k_sg++;
            }
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
                // static Options options;
                // Options_setToDefault( &options );
                // options.initialStatusBounds = ST_INACTIVE;
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
                    QProblem_setOptions(QP, options);
                }
                if (opts->warm_start)
                {
                    qpoases_status = (ns > 0) ?
                        QProblem_initW(QP, HH, gg, CC, d_lb, d_ub, d_lg, d_ug, &nwsr,
                                       &cputime, NULL, dual_sol, NULL, NULL, NULL) :
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
    // copy prim_sol and dual_sol to qpd_sol
    blasfeo_pack_dvec(nv2, prim_sol, qp_out->v, 0);
    for (int ii = 0; ii < 2 * nb + 2 * ng + 2 * ns; ii++) qp_out->lam->pa[ii] = 0.0;
    for (int ii = 0; ii < nb; ii++)
    {
        if (dual_sol[idxb[ii]] >= 0.0)
            qp_out->lam->pa[ii] = dual_sol[idxb[ii]];
        else
            qp_out->lam->pa[nb + ng + ii] = -dual_sol[idxb[ii]];
    }

    int k = 0;
    for (int ii = 0; ii < ns; ii++)
    {
        int js = idxs[ii];

        if (js < nb)
        {
            // dual variables for softened upper box constraints
            if (dual_sol[2*ng+k] <= 0.0)
                qp_out->lam->pa[nb+ng+js] = -dual_sol[2*ng+k];

            // dual variables for softened lower box constraints
            if (dual_sol[2*ng+nsb+k] <= 0.0)
                qp_out->lam->pa[js] = -dual_sol[2*ng+nsb+k];

            k++;
        }

        // dual variables for sl >= d_ls
        if (dual_sol[nv+ii] >= 0)
            qp_out->lam->pa[2*nb+2*ng+ii] = dual_sol[nv+ii];

        // dual variables for su >= d_us
        if (dual_sol[nv+ns+ii] >= 0)
            qp_out->lam->pa[2*nb+2*ng+ns+ii] = dual_sol[nv+ns+ii];
    }

    for (int ii = 0; ii < ng; ii++)
    {
        if (dual_sol[nv2 + ii] >= 0.0)
            qp_out->lam->pa[nb + ii] = dual_sol[nv2 + ii];
        else
            qp_out->lam->pa[2 * nb + ng + ii] = -dual_sol[nv2 + ii];
    }
    info->interface_time += acados_toc(&interface_timer);
    info->total_time = acados_toc(&tot_timer);
    info->num_iter = nwsr;

    // compute slacks
    if (opts->compute_t)
    {
        blasfeo_dvecex_sp(nb, 1.0, qp_in->idxb, qp_out->v, 0, qp_out->t, nb + ng);
        blasfeo_dgemv_t(nv, ng, 1.0, qp_in->Ct, 0, 0, qp_out->v, 0, 0.0, qp_out->t, 2 * nb + ng,
                        qp_out->t, 2 * nb + ng);
        blasfeo_dveccpsc(nb + ng, -1.0, qp_out->t, nb + ng, qp_out->t, 0);
        blasfeo_daxpy(2 * nb + 2 * ng, -1.0, qp_in->d, 0, qp_out->t, 0, qp_out->t, 0);
    }

    int acados_status = qpoases_status;
    if (qpoases_status == SUCCESSFUL_RETURN) acados_status = ACADOS_SUCCESS;
    if (qpoases_status == RET_MAX_NWSR_REACHED) acados_status = ACADOS_MAXITER;
    return acados_status;
}

void dense_qp_qpoases_config_initialize_default(void *config_)
{
    qp_solver_config *config = config_;

    config->opts_calculate_size = (int (*)(void *, void *)) & dense_qp_qpoases_opts_calculate_size;
    config->opts_assign = (void *(*) (void *, void *, void *) ) & dense_qp_qpoases_opts_assign;
    config->opts_initialize_default =
        (void (*)(void *, void *, void *)) & dense_qp_qpoases_opts_initialize_default;
    config->opts_update = (void (*)(void *, void *, void *)) & dense_qp_qpoases_opts_update;
    config->memory_calculate_size =
        (int (*)(void *, void *, void *)) & dense_qp_qpoases_memory_calculate_size;
    config->memory_assign =
        (void *(*) (void *, void *, void *, void *) ) & dense_qp_qpoases_memory_assign;
    config->workspace_calculate_size =
        (int (*)(void *, void *, void *)) & dense_qp_qpoases_workspace_calculate_size;
    config->evaluate = (int (*)(void *, void *, void *, void *, void *, void *)) & dense_qp_qpoases;

    return;
}
