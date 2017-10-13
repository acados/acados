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

#include "acados/ocp_qp/ocp_qp_condensing_qpoases.h"

#include <stdlib.h>

#include "acados/utils/math.h"

#include "blasfeo/include/blasfeo_common.h"
#include "blasfeo/include/blasfeo_d_aux_ext_dep.h"
#include "blasfeo/include/blasfeo_i_aux_ext_dep.h"
#include "blasfeo/include/blasfeo_target.h"

/* Ignore compiler warnings from qpOASES */
#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wtautological-pointer-compare"
#pragma clang diagnostic ignored "-Wunused-parameter"
#pragma clang diagnostic ignored "-Wunused-function"
#include "qpOASES_e/QProblemB.h"
#include "qpOASES_e/QProblem.h"
#pragma clang diagnostic pop
#elif defined(__GNUC__)
    #if __GNUC__ >= 6
        #pragma GCC diagnostic push
        #pragma GCC diagnostic ignored "-Wunused-but-set-parameter"
        #pragma GCC diagnostic ignored "-Wunused-parameter"
        #pragma GCC diagnostic ignored "-Wunused-function"
        #include "qpOASES_e/QProblemB.h"
        #include "qpOASES_e/QProblem.h"
        #pragma GCC diagnostic pop
    #else
        #pragma GCC diagnostic ignored "-Wunused-parameter"
        #pragma GCC diagnostic ignored "-Wunused-function"
        #include "qpOASES_e/QProblemB.h"
        #include "qpOASES_e/QProblem.h"
    #endif
#else
    #include "qpOASES_e/QProblemB.h"
    #include "qpOASES_e/QProblem.h"
#endif

#include "blasfeo/include/blasfeo_d_aux.h"
#include "blasfeo/include/blasfeo_d_blas.h"
#include "hpipm/include/hpipm_d_cond.h"
#include "hpipm/include/hpipm_d_dense_qp.h"
#include "hpipm/include/hpipm_d_dense_qp_sol.h"
#include "hpipm/include/hpipm_d_ocp_qp.h"
#include "hpipm/include/hpipm_d_ocp_qp_sol.h"



int ocp_qp_condensing_qpoases_calculate_args_size(const ocp_qp_in *qp_in) {
    int N = qp_in->N;
    int size = 0;
    size += sizeof(ocp_qp_condensing_qpoases_args);
    size += (N+1)*sizeof(int);
    return size;
}



char *ocp_qp_condensing_qpoases_assign_args(const ocp_qp_in *qp_in,
    ocp_qp_condensing_qpoases_args **args, void *mem) {

    int N = qp_in->N;

    char *c_ptr = (char *) mem;

    *args = (ocp_qp_condensing_qpoases_args *) c_ptr;
    c_ptr += sizeof(ocp_qp_condensing_qpoases_args);

    (*args)->scrapspace = c_ptr;
    c_ptr += (N+1)*sizeof(int);

    return c_ptr;
    }



static void ocp_qp_condensing_qpoases_initialize_default_args(const ocp_qp_in *qp_in,
    ocp_qp_condensing_qpoases_args *args) {

    args->cputime = 1000.0;  // maximum cpu time in seconds
    args->warm_start = 0;
    args->nwsr = 1000;

    int N = qp_in->N;

    int *ns = (int *) args->scrapspace;
    int ii;
    for (ii=0; ii < N+1; ii++) ns[ii] = 0;
    }



ocp_qp_condensing_qpoases_args *ocp_qp_condensing_qpoases_create_arguments(const ocp_qp_in *qp_in) {
    void *mem = malloc(ocp_qp_condensing_qpoases_calculate_args_size(qp_in));
    ocp_qp_condensing_qpoases_args *args;
    ocp_qp_condensing_qpoases_assign_args(qp_in, &args, mem);
    ocp_qp_condensing_qpoases_initialize_default_args(qp_in, args);

    return args;
}



int ocp_qp_condensing_qpoases_calculate_workspace_size(const ocp_qp_in *qp_in, void *args_) {
    return 0;
}



int ocp_qp_condensing_qpoases_calculate_memory_size(const ocp_qp_in *qp_in, void *args_) {
    ocp_qp_condensing_qpoases_args *args = (ocp_qp_condensing_qpoases_args *) args_;
    // extract ocp qp in size
    int N = qp_in->N;
    int *nx = (int *)qp_in->nx;
    int *nu = (int *)qp_in->nu;
    int *nb = (int *)qp_in->nb;
    int *ng = (int *)qp_in->nc;
    int **hidxb = (int **)qp_in->idxb;

    // extract ns from args
    int_t *ns = (int_t *) args->scrapspace;

    // dummy ocp qp
    struct d_ocp_qp qp;
    qp.N = N;
    qp.nx = nx;
    qp.nu = nu;
    qp.nb = nb;
    qp.ng = ng;
    qp.idxb = hidxb;
    qp.ns = ns;

    // compute dense qp size
    int nvd = 0;
    int ned = 0;
    int nbd = 0;
    int ngd = 0;
    int_t nsd = 0;
#if 0
    // [u; x] order
    d_compute_qp_size_ocp2dense(N, nx, nu, nb, hidxb, ng, &nvd, &ned, &nbd, &ngd);
#else
    // [x; u] order  // XXX update with ns !!!!!
    d_compute_qp_size_ocp2dense_rev(N, nx, nu, nb, hidxb, ng, &nvd, &ned, &nbd, &ngd);
#endif

    // dummy dense qp
    struct d_dense_qp qpd;
    qpd.nv = nvd;
    qpd.ne = ned;
    qpd.nb = nbd;
    qpd.ng = ngd;
    qpd.ns = nsd;

    // size in bytes
    int size = sizeof(ocp_qp_condensing_qpoases_memory);

    size += 1 * sizeof(struct d_ocp_qp);                       // qp
    size += 1 * sizeof(struct d_ocp_qp_sol);                   // qp_sol
    size += 1 * sizeof(struct d_dense_qp);                     // qpd
    size += 1 * sizeof(struct d_dense_qp_sol);                 // qpd_sol
    size += 1 * sizeof(struct d_cond_qp_ocp2dense_workspace);  // cond_workspace
    //  size += 1*sizeof(struct d_strmat); // sR

    size += d_memsize_ocp_qp(N, nx, nu, nb, ng, ns);
    size += d_memsize_ocp_qp_sol(N, nx, nu, nb, ng, ns);
    size += d_memsize_dense_qp(nvd, ned, nbd, ngd, nsd);
    size += d_memsize_dense_qp_sol(nvd, ned, nbd, ngd, nsd);
    size += d_memsize_cond_qp_ocp2dense(&qp, &qpd);
    size += 4 * (N + 1) * sizeof(double *);  // lam_lb lam_ub lam_lg lam_ug
    size += 1 * (N + 1) * sizeof(int *);  // hidxb_rev
    for (int ii = 0; ii <= N; ii++) {
        size += nb[ii]*sizeof(int);  // hidxb_rev
    }
    //  size += 1*d_size_strmat(nvd, nvd); // sR

    size += 1 * nvd * nvd * sizeof(double);  // H
    //  size += 1*nvd*nvd*sizeof(double); // R
    size += 1 * nvd * ned * sizeof(double);        // A
    size += 1 * nvd * ngd * sizeof(double);        // C
    size += 3 * nvd * sizeof(double);              // g d_lb d_ub
    size += 1 * ned * sizeof(double);              // b
    size += 2 * nbd * sizeof(double);              // d_lb0 d_ub0
    size += 2 * ngd * sizeof(double);              // d_lg d_ug
    size += 1 * nbd * sizeof(int);                 // idxb
    size += 1 * nvd * sizeof(double);              // prim_sol
    size += (2 * nvd + 2 * ngd) * sizeof(double);  // dual_sol

    if (ngd > 0)  // QProblem
        size += sizeof(QProblem);
    else  // QProblemB
        size += sizeof(QProblemB);

    size = (size + 63) / 64 * 64;  // make multipl of typical cache line size
    size += 1 * 64;                // align once to typical cache line size

    return size;
}



char *ocp_qp_condensing_qpoases_assign_memory(const ocp_qp_in *qp_in, void *args_,
                                             void **mem, void *raw_memory) {

    ocp_qp_condensing_qpoases_args *args = (ocp_qp_condensing_qpoases_args *) args_;
    ocp_qp_condensing_qpoases_memory **qpoases_memory = (ocp_qp_condensing_qpoases_memory **) mem;

    // extract problem size
    int N = qp_in->N;
    int *nx = (int *)qp_in->nx;
    int *nu = (int *)qp_in->nu;
    int *nb = (int *)qp_in->nb;
    int *ng = (int *)qp_in->nc;
    int **hidxb = (int **)qp_in->idxb;

    // extract ns from args
    int_t *ns = (int_t *) args->scrapspace;

    // compute dense qp size
    int nvd = 0;
    int ned = 0;
    int nbd = 0;
    int ngd = 0;
    int_t nsd = 0;
#if 0
    // [u; x] order
    d_compute_qp_size_ocp2dense(N, nx, nu, nb, hidxb, ng, &nvd, &ned, &nbd, &ngd);
#else
    // [x; u] order  // XXX update with ns !!!!!
    d_compute_qp_size_ocp2dense_rev(N, nx, nu, nb, hidxb, ng, &nvd, &ned, &nbd, &ngd);
#endif


    // char pointer
    char *c_ptr = (char *)raw_memory;

    *qpoases_memory = (ocp_qp_condensing_qpoases_memory *) c_ptr;
    c_ptr += sizeof(ocp_qp_condensing_qpoases_memory);

    //
    //  (*qpoases_memory)->sR = (struct d_strmat *) c_ptr;
    //  c_ptr += 1*sizeof(struct d_strmat);

    //
    (*qpoases_memory)->qp = (struct d_ocp_qp *)c_ptr;
    c_ptr += 1 * sizeof(struct d_ocp_qp);
    //
    (*qpoases_memory)->qp_sol = (struct d_ocp_qp_sol *)c_ptr;
    c_ptr += 1 * sizeof(struct d_ocp_qp_sol);
    //
    (*qpoases_memory)->qpd = (struct d_dense_qp *)c_ptr;
    c_ptr += 1 * sizeof(struct d_dense_qp);
    //
    (*qpoases_memory)->qpd_sol = (struct d_dense_qp_sol *)c_ptr;
    c_ptr += 1 * sizeof(struct d_dense_qp_sol);
    //
    (*qpoases_memory)->cond_workspace =
        (struct d_cond_qp_ocp2dense_workspace *)c_ptr;
    c_ptr += 1 * sizeof(struct d_cond_qp_ocp2dense_workspace);
    //
    (*qpoases_memory)->hlam_lb = (double **)c_ptr;
    c_ptr += (N + 1) * sizeof(double *);
    //
    (*qpoases_memory)->hlam_ub = (double **)c_ptr;
    c_ptr += (N + 1) * sizeof(double *);
    //
    (*qpoases_memory)->hlam_lg = (double **)c_ptr;
    c_ptr += (N + 1) * sizeof(double *);
    //
    (*qpoases_memory)->hlam_ug = (double **)c_ptr;
    c_ptr += (N + 1) * sizeof(double *);
    //
    (*qpoases_memory)->hidxb_rev = (int **)c_ptr;
    c_ptr += (N + 1) * sizeof(int *);

    //
    //  struct d_strmat *sR = (*qpoases_memory)->sR;

    //
    struct d_ocp_qp *qp = (*qpoases_memory)->qp;
    //
    struct d_ocp_qp_sol *qp_sol = (*qpoases_memory)->qp_sol;
    //
    struct d_dense_qp *qpd = (*qpoases_memory)->qpd;
    //
    struct d_dense_qp_sol *qpd_sol = (*qpoases_memory)->qpd_sol;
    //
    struct d_cond_qp_ocp2dense_workspace *cond_workspace =
        (*qpoases_memory)->cond_workspace;

    // align memory to typical cache line size
    size_t s_ptr = (size_t)c_ptr;
    s_ptr = (s_ptr + 63) / 64 * 64;
    c_ptr = (char *)s_ptr;

    //
    //  d_create_strmat(nvd, nvd, sR, c_ptr);
    //  c_ptr += sR->memory_size;

    // ocp qp structure
    d_create_ocp_qp(N, nx, nu, nb, ng, ns, qp, c_ptr);
    c_ptr += qp->memsize;
    // ocp qp sol structure
    d_create_ocp_qp_sol(N, nx, nu, nb, ng, ns, qp_sol, c_ptr);
    c_ptr += qp_sol->memsize;
    // dense qp structure
    d_create_dense_qp(nvd, ned, nbd, ngd, nsd, qpd, c_ptr);
    c_ptr += qpd->memsize;
    // dense qp sol structure
    d_create_dense_qp_sol(nvd, ned, nbd, ngd, nsd, qpd_sol, c_ptr);
    c_ptr += qpd_sol->memsize;
    // cond workspace structure
    d_create_cond_qp_ocp2dense(qp, qpd, cond_workspace, c_ptr);
    c_ptr += cond_workspace->memsize;

    // double stuff
    //
    (*qpoases_memory)->H = (double *)c_ptr;
    c_ptr += nvd * nvd * sizeof(double);
    //
    //  (*qpoases_memory)->R = (double *) c_ptr;
    //  c_ptr += nvd*nvd*sizeof(double);
    //
    (*qpoases_memory)->A = (double *)c_ptr;
    c_ptr += nvd * ned * sizeof(double);
    //
    (*qpoases_memory)->C = (double *)c_ptr;
    c_ptr += nvd * ngd * sizeof(double);
    //
    (*qpoases_memory)->g = (double *)c_ptr;
    c_ptr += nvd * sizeof(double);
    //
    (*qpoases_memory)->b = (double *)c_ptr;
    c_ptr += ned * sizeof(double);
    //
    (*qpoases_memory)->d_lb0 = (double *)c_ptr;
    c_ptr += nbd * sizeof(double);
    //
    (*qpoases_memory)->d_ub0 = (double *)c_ptr;
    c_ptr += nbd * sizeof(double);
    //
    (*qpoases_memory)->d_lb = (double *)c_ptr;
    c_ptr += nvd * sizeof(double);
    //
    (*qpoases_memory)->d_ub = (double *)c_ptr;
    c_ptr += nvd * sizeof(double);
    //
    (*qpoases_memory)->d_lg = (double *)c_ptr;
    c_ptr += ngd * sizeof(double);
    //
    (*qpoases_memory)->d_ug = (double *)c_ptr;
    c_ptr += ngd * sizeof(double);
    //
    (*qpoases_memory)->prim_sol = (double *)c_ptr;
    c_ptr += nvd * sizeof(double);
    //
    (*qpoases_memory)->dual_sol = (double *)c_ptr;
    c_ptr += (2 * nvd + 2 * ngd) * sizeof(double);
    //

    // qpOASES (HUGE!!!)
    //
    if (ngd > 0) {  // QProblem
        (*qpoases_memory)->QP = (void *)c_ptr;
        c_ptr += sizeof(QProblem);
    } else {  // QProblemB
        (*qpoases_memory)->QPB = (void *)c_ptr;
        c_ptr += sizeof(QProblemB);
    }

    // int stuff
    //
    (*qpoases_memory)->idxb = (int *)c_ptr;
    c_ptr += nbd * sizeof(int);
    //
    for (int ii = 0; ii <= N; ii++) {
        (*qpoases_memory)->hidxb_rev[ii] = (int *) c_ptr;
        c_ptr += nb[ii]*sizeof(int);
    }

    return c_ptr;
}



ocp_qp_condensing_qpoases_memory *ocp_qp_condensing_qpoases_create_memory(const ocp_qp_in *qp_in,
                                                                          void *args_) {
    ocp_qp_condensing_qpoases_args *args = (ocp_qp_condensing_qpoases_args *) args_;
    ocp_qp_condensing_qpoases_memory *mem;
    int_t memory_size = ocp_qp_condensing_qpoases_calculate_memory_size(qp_in, args);
    void *raw_memory_ptr = malloc(memory_size);
    char *ptr_end =
        ocp_qp_condensing_qpoases_assign_memory(qp_in, args, (void **) &mem, raw_memory_ptr);
    assert((char*)raw_memory_ptr + memory_size >= ptr_end); (void) ptr_end;
    return mem;
}



int ocp_qp_condensing_qpoases(const ocp_qp_in *qp_in, ocp_qp_out *qp_out, void *args_,
                              void *memory_, void *workspace_) {
    // cast structures
    ocp_qp_condensing_qpoases_args *args =
        (ocp_qp_condensing_qpoases_args *)args_;
    ocp_qp_condensing_qpoases_memory *memory =
        (ocp_qp_condensing_qpoases_memory *)memory_;

    // initialize return code
    int acados_status = ACADOS_SUCCESS;

    // loop index
    int ii, jj;

    // extract memory
    double **hlam_lb = memory->hlam_lb;
    double **hlam_ub = memory->hlam_ub;
    double **hlam_lg = memory->hlam_lg;
    double **hlam_ug = memory->hlam_ug;
    //  struct d_strmat *sR = memory->sR;
    struct d_ocp_qp *qp = memory->qp;
    struct d_ocp_qp_sol *qp_sol = memory->qp_sol;
    struct d_dense_qp *qpd = memory->qpd;
    struct d_dense_qp_sol *qpd_sol = memory->qpd_sol;
    struct d_cond_qp_ocp2dense_workspace *cond_workspace =
        memory->cond_workspace;
    double *H = memory->H;
    //  double *R = memory->R;
    double *A = memory->A;
    double *C = memory->C;
    double *g = memory->g;
    double *b = memory->b;
    double *d_lb0 = memory->d_lb0;
    double *d_ub0 = memory->d_ub0;
    double *d_lb = memory->d_lb;
    double *d_ub = memory->d_ub;
    double *d_lg = memory->d_lg;
    double *d_ug = memory->d_ug;
    int *idxb = memory->idxb;
    int **hidxb_rev = (int **) memory->hidxb_rev;
    double *prim_sol = memory->prim_sol;
    double *dual_sol = memory->dual_sol;
    QProblemB *QPB = memory->QPB;
    QProblem *QP = memory->QP;

    // extract ocp problem size
    int N = qp_in->N;
    int *nx = (int *) qp_in->nx;
    int *nu = (int *) qp_in->nu;
    int *nb = (int *)qp_in->nb;
    int *ng = (int *)qp_in->nc;

    // extract input data
    double **hA = (double **)qp_in->A;
    double **hB = (double **)qp_in->B;
    double **hb = (double **)qp_in->b;
    double **hQ = (double **)qp_in->Q;
    double **hS = (double **)qp_in->S;
    double **hR = (double **)qp_in->R;
    double **hq = (double **)qp_in->q;
    double **hr = (double **)qp_in->r;
    double **hd_lb = (double **)qp_in->lb;
    double **hd_ub = (double **)qp_in->ub;
    double **hC = (double **)qp_in->Cx;
    double **hD = (double **)qp_in->Cu;
    double **hd_lg = (double **)qp_in->lc;
    double **hd_ug = (double **)qp_in->uc;
    int **hidxb = (int **)qp_in->idxb;

    // extract output struct members
    double **hx = qp_out->x;
    double **hu = qp_out->u;
    double **hpi = qp_out->pi;
    double **hlam = qp_out->lam;

    // compute bounds indeces in order [u; x]
    for (ii = 0; ii <= N; ii++) {
        for (jj = 0; jj < nb[ii]; jj++) {
            if (hidxb[ii][jj] < nx[ii]) {  // state constraint
                hidxb_rev[ii][jj] = hidxb[ii][jj]+nu[ii];
            } else  {  // input constraint
                hidxb_rev[ii][jj] = hidxb[ii][jj]-nx[ii];
            }
        }
    }

    //
    for (ii = 0; ii <= N; ii++) {
        hlam_lb[ii] = hlam[ii];
        hlam_ub[ii] = hlam[ii] + nb[ii];
        hlam_lg[ii] = hlam[ii] + 2 * nb[ii];
        hlam_ug[ii] = hlam[ii] + 2 * nb[ii] + ng[ii];
    }

    // extract dense qp size
    int nvd = qpd->nv;
    //  int ned = qpd->ne;
    int nbd = qpd->nb;
    int ngd = qpd->ng;

    // ocp qp structure
    d_cvt_colmaj_to_ocp_qp(hA, hB, hb, hQ, hS, hR, hq, hr, hidxb_rev, hd_lb, hd_ub, hC, hD,
        hd_lg, hd_ug, NULL, NULL, NULL, NULL, NULL, qp);

    // dense qp structure
    d_cond_qp_ocp2dense(qp, qpd, cond_workspace);

#if 0
    d_print_strmat(nvd, nvd, qpd->Hg, 0, 0);
    exit(1);
#endif

    // fill in the upper triangular of H in dense_qp
    dtrtr_l_libstr(nvd, qpd->Hg, 0, 0, qpd->Hg, 0, 0);

    // dense qp row-major
    d_cvt_dense_qp_to_rowmaj(qpd, H, g, A, b, idxb, d_lb0, d_ub0, C, d_lg, d_ug,
        NULL, NULL, NULL, NULL, NULL);

    // reorder bounds
    for (ii = 0; ii < nvd; ii++) {
        d_lb[ii] = -QPOASES_INFTY;
        d_ub[ii] = +QPOASES_INFTY;
    }
    for (ii = 0; ii < nbd; ii++) {
        d_lb[idxb[ii]] = d_lb0[ii];
        d_ub[idxb[ii]] = d_ub0[ii];
    }

// cholesky factorization of H
//  dpotrf_l_libstr(nvd, qpd->Hg, 0, 0, sR, 0, 0);

//  fill in upper triangular of R
//  dtrtr_l_libstr(nvd, sR, 0, 0, sR, 0, 0);

//  extract R
//  d_cvt_strmat2mat(nvd, nvd, sR, 0, 0, R, nvd);

#if 0
    d_print_mat(nvd, nvd, H, nvd);
    d_print_mat(nvd, nvd, R, nvd);
    exit(1);
#endif

    // cold start the dual solution with no active constraints
    int warm_start = args->warm_start;
    if (!warm_start) {
        for (ii = 0; ii < 2 * nvd + 2 * ngd; ii++)
            dual_sol[ii] = 0;
    }

    // solve dense qp
    int nwsr = args->nwsr;  // max number of working set recalculations
    double cputime = args->cputime;
    int return_flag = 0;
    if (ngd > 0) {  // QProblem
        QProblemCON(QP, nvd, ngd, HST_POSDEF);
        QProblem_setPrintLevel(QP, PL_MEDIUM);
        QProblem_printProperties(QP);
        return_flag =
            QProblem_initW(QP, H, g, C, d_lb, d_ub, d_lg, d_ug, &nwsr, &cputime,
                NULL, dual_sol, NULL, NULL, NULL);  // NULL or 0
        //            NULL, NULL, NULL, NULL);
        //            NULL, NULL, NULL, R);  // to provide Cholesky factor
        QProblem_getPrimalSolution(QP, prim_sol);
        QProblem_getDualSolution(QP, dual_sol);
    } else {  // QProblemB
        QProblemBCON(QPB, nvd, HST_POSDEF);
        QProblemB_setPrintLevel(QPB, PL_MEDIUM);
        QProblemB_printProperties(QPB);
        return_flag = QProblemB_initW(QPB, H, g, d_lb, d_ub, &nwsr, &cputime,
            NULL, dual_sol, NULL, NULL);  // NULL or 0
        QProblemB_getPrimalSolution(QPB, prim_sol);
        QProblemB_getDualSolution(QPB, dual_sol);
    }

    // save solution statistics to memory
    memory->cputime = cputime;
    memory->nwsr = nwsr;

#if 0
    d_print_mat(1, nvd, prim_sol, 1);
    exit(1);
#endif

    // copy prim_sol and dual_sol to qpd_sol
    d_cvt_vec2strvec(nvd, prim_sol, qpd_sol->v, 0);
    for (ii=0; ii < 2*nbd+2*ngd; ii++) qpd_sol->lam->pa[ii] = 0.0;
    for (ii=0; ii < nbd; ii++)
        if (dual_sol[ii] >= 0.0)
            qpd_sol->lam->pa[ii] =   dual_sol[ii];
        else
            qpd_sol->lam->pa[nbd+ngd+ii] = - dual_sol[ii];
    for (ii=0; ii < ngd; ii++)
        if (dual_sol[nbd+ii] >= 0.0)
            qpd_sol->lam->pa[nbd+ii] =   dual_sol[nbd+ii];
        else
            qpd_sol->lam->pa[2*nbd+ngd+ii] = - dual_sol[nbd+ii];

    // expand solution
    d_expand_sol_dense2ocp(qp, qpd_sol, qp_sol, cond_workspace);

    // extract solution
    d_cvt_ocp_qp_sol_to_colmaj(qp, qp_sol, hu, hx, NULL, NULL, hpi,
                               hlam_lb, hlam_ub, hlam_lg, hlam_ug, NULL, NULL);

    // return
    acados_status = return_flag;
    return acados_status;
    //
}



void ocp_qp_condensing_qpoases_initialize(const ocp_qp_in *qp_in, void *args_, void **mem,
                                          void **work) {
    ocp_qp_condensing_qpoases_args *args = (ocp_qp_condensing_qpoases_args *)args_;

    *mem = ocp_qp_condensing_qpoases_create_memory(qp_in, args);
    *work = NULL;
}



void ocp_qp_condensing_qpoases_destroy(void *mem_, void *work) {
    // TODO(dimitris): replace dummy commands once interface completed
    (void) mem_;
    (void) work;
}
