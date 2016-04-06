#include <stdlib.h>
#include <math.h>

#include "ocp_qp_solver.h"

#include "target.h"
#include "block_size.h"
#include "aux_d.h"
#include "aux_s.h"
#include "blas_d.h"
#include "lqcp_solvers.h"
#include "mpc_solvers.h"

// work space size
int ocp_qp_hpmpc_workspace_size(int N, int *nx, int *nu, int *ng, int *nb,
                                struct ocp_qp_hpmpc_args *args) {
    const int bs = D_MR;
    const int ncl = D_NCL;

    int kmax = args->max_iter;

    int ii;

    int pnx[N + 1];
    for (ii = 0; ii <= N; ii++) pnx[ii] = (nx[ii] + bs - 1) / bs * bs;
    int pnz[N + 1];
    for (ii = 0; ii <= N; ii++)
        pnz[ii] = (nu[ii] + nx[ii] + 1 + bs - 1) / bs * bs;
    int pnux[N + 1];
    for (ii = 0; ii <= N; ii++) pnux[ii] = (nu[ii] + nx[ii] + bs - 1) / bs * bs;
    int pnb[N + 1];
    for (ii = 0; ii <= N; ii++) pnb[ii] = (nb[ii] + bs - 1) / bs * bs;
    int png[N + 1];
    for (ii = 0; ii <= N; ii++) png[ii] = (ng[ii] + bs - 1) / bs * bs;
    int cnx[N + 1];
    for (ii = 0; ii <= N; ii++) cnx[ii] = (nx[ii] + ncl - 1) / ncl * ncl;
    int cnux[N + 1];
    for (ii = 0; ii <= N; ii++)
        cnux[ii] = (nu[ii] + nx[ii] + ncl - 1) / ncl * ncl;
    int cng[N + 1];
    for (ii = 0; ii <= N; ii++) cng[ii] = (ng[ii] + ncl - 1) / ncl * ncl;

    int work_space = 16;  // enough to align twice to cache boundaries

    for (ii = 0; ii < N; ii++) {
        work_space += pnz[ii] * cnx[ii + 1] + pnz[ii] * cnux[ii] +
                      pnux[ii] * cng[ii] + 6 * pnb[ii] + 6 * png[ii] +
                      pnux[ii] + pnx[ii];
    }
    ii = N;
    work_space += pnz[ii] * cnux[ii] + pnux[ii] * cng[ii] + 6 * pnb[ii] +
                  6 * png[ii] + pnux[ii] + pnx[ii];

    work_space += 5 * kmax;

    work_space += d_ip2_hard_mpc_tv_work_space_size_double(
        N, nx, nu, nb, ng);  // work space needed by the IPM

    return work_space * sizeof(double);
}

/* version dealing with equality constratins: is lb=ub, then fix the variable
 * (corresponding column in A or B set to zero, and updated b) */
int ocp_qp_hpmpc(int N, int *nx, int *nu, int *nb, int *ng, double **hA,
                 double **hB, double **hb, double **hQ, double **hS,
                 double **hR, double **hq, double **hr, int **hidxb,
                 double **hlb, double **hub, double **hC, double **hD,
                 double **hlg, double **hug, double **hx, double **hu,
                 struct ocp_qp_hpmpc_args *args, double *work0) {
    //    printf("\nstart of wrapper\n");

    int acados_status = ACADOS_SUCCESS;

    const int bs = D_MR;
    const int ncl = D_NCL;

    // const int anb = nal*((2*nb+nal-1)/nal);

    double mu_tol = args->tol;
    int k_max = args->max_iter;
    double alpha_min = args->min_step;  // minimum accepted step length
    double sigma_par[] = {0.4, 0.1,
                          args->sigma_min};  // control primal-dual IP behaviour
    double mu0 = args->mu0;
    int warm_start = 0;
    int compute_mult = 1;  // compute multipliers TODO(giaf): set to zero
    int kk = -1;

    int ii, jj, ll;

    // align work space (assume cache line <= 64 byte)
    size_t addr = (((size_t)work0) + 63) / 64 * 64;
    double *ptr = (double *)addr;

    // array of pointers
    double *(hpBAbt[N]);
    double *(hpRSQrq[N + 1]);
    double *(hpDCt[N + 1]);
    double *(hd[N + 1]);
    double *(hux[N + 1]);
    double *(hpi[N + 1]);
    double *(hlam[N + 1]);
    double *(ht[N + 1]);
    // double *(hqq[N + 1]);
    // double *(hrb[N]);
    // double *(hrq[N + 1]);
    // double *(hrd[N + 1]);
    double *stat;
    double *work;

    // matrices size // TODO(giaf): pass them to the solver !!!!!
    int pnx[N + 1];
    for (ii = 0; ii <= N; ii++) pnx[ii] = (nx[ii] + bs - 1) / bs * bs;
    int pnz[N + 1];
    for (ii = 0; ii <= N; ii++)
        pnz[ii] = (nu[ii] + nx[ii] + 1 + bs - 1) / bs * bs;
    int pnux[N + 1];
    for (ii = 0; ii <= N; ii++) pnux[ii] = (nu[ii] + nx[ii] + bs - 1) / bs * bs;
    int pnb[N + 1];
    for (ii = 0; ii <= N; ii++) pnb[ii] = (nb[ii] + bs - 1) / bs * bs;
    int png[N + 1];
    for (ii = 0; ii <= N; ii++) png[ii] = (ng[ii] + bs - 1) / bs * bs;
    int cnx[N + 1];
    for (ii = 0; ii <= N; ii++) cnx[ii] = (nx[ii] + ncl - 1) / ncl * ncl;
    int cnux[N + 1];
    for (ii = 0; ii <= N; ii++)
        cnux[ii] = (nu[ii] + nx[ii] + ncl - 1) / ncl * ncl;
    int cng[N + 1];
    for (ii = 0; ii <= N; ii++) cng[ii] = (ng[ii] + ncl - 1) / ncl * ncl;

    // data structure

    ii = 0;
    hpBAbt[ii] = ptr;
    ptr += pnz[ii] * cnx[ii + 1];
    d_cvt_tran_mat2pmat(nx[ii + 1], nu[ii], hB[ii], nx[ii + 1], 0, hpBAbt[ii],
                        cnx[ii + 1]);
    d_cvt_tran_mat2pmat(
        nx[ii + 1], nx[ii], hA[ii], nx[ii + 1], nu[ii],
        hpBAbt[ii] + nu[ii] / bs * bs * cnx[ii + 1] + nu[ii] % bs, cnx[ii + 1]);
    d_cvt_tran_mat2pmat(nx[ii + 1], 1, hb[ii], nx[ii + 1], nu[ii] + nx[ii],
                        hpBAbt[ii] + (nu[ii] + nx[ii]) / bs * bs * cnx[ii + 1] +
                            (nu[ii] + nx[ii]) % bs,
                        cnx[ii + 1]);
    //    d_print_pmat(nu[ii]+nx[ii]+1, nx[ii+1], bs, hpBAbt[ii], cnx[ii+1]);
    for (ii = 1; ii < N; ii++) {
        if (hb[ii] == hb[ii - 1] && hA[ii] == hA[ii - 1] &&
            hB[ii] == hB[ii - 1]) {
            //            printf("\nsame\n");
            hpBAbt[ii] = hpBAbt[ii - 1];
        } else {
            //            printf("\ndifferent\n");
            hpBAbt[ii] = ptr;
            ptr += pnz[ii] * cnx[ii + 1];
            d_cvt_tran_mat2pmat(nx[ii + 1], nu[ii], hB[ii], nx[ii + 1], 0,
                                hpBAbt[ii], cnx[ii + 1]);
            d_cvt_tran_mat2pmat(
                nx[ii + 1], nx[ii], hA[ii], nx[ii + 1], nu[ii],
                hpBAbt[ii] + nu[ii] / bs * bs * cnx[ii + 1] + nu[ii] % bs,
                cnx[ii + 1]);
            d_cvt_tran_mat2pmat(
                nx[ii + 1], 1, hb[ii], nx[ii + 1], nu[ii] + nx[ii],
                hpBAbt[ii] + (nu[ii] + nx[ii]) / bs * bs * cnx[ii + 1] +
                    (nu[ii] + nx[ii]) % bs,
                cnx[ii + 1]);
        }
        //        d_print_pmat(nu[ii]+nx[ii]+1, nx[ii+1], bs, hpBAbt[ii],
        //        cnx[ii+1]);
    }

    ii = 0;
    hpRSQrq[ii] = ptr;
    ptr += pnz[ii] * cnux[ii];
    d_cvt_mat2pmat(nu[ii], nu[ii], hR[ii], nu[ii], 0, hpRSQrq[ii], cnux[ii]);
    d_cvt_tran_mat2pmat(nx[ii], nu[ii], hS[ii], nu[ii], nu[ii],
                        hpRSQrq[ii] + nu[ii] / bs * bs * cnux[ii] + nu[ii] % bs,
                        cnux[ii]);
    d_cvt_mat2pmat(
        nx[ii], nx[ii], hQ[ii], nx[ii], nu[ii],
        hpRSQrq[ii] + nu[ii] / bs * bs * cnux[ii] + nu[ii] % bs + nu[ii] * bs,
        cnux[ii]);
    d_cvt_tran_mat2pmat(nu[ii], 1, hr[ii], nu[ii], nu[ii] + nx[ii],
                        hpRSQrq[ii] + (nu[ii] + nx[ii]) / bs * bs * cnux[ii] +
                            (nu[ii] + nx[ii]) % bs,
                        cnux[ii]);
    d_cvt_tran_mat2pmat(nx[ii], 1, hq[ii], nx[ii], nu[ii] + nx[ii],
                        hpRSQrq[ii] + (nu[ii] + nx[ii]) / bs * bs * cnux[ii] +
                            (nu[ii] + nx[ii]) % bs + nu[ii] * bs,
                        cnux[ii]);
    //    d_print_pmat(nu[ii]+nx[ii]+1, nu[ii]+nx[ii], bs, hpRSQrq[ii],
    //    cnux[ii]);
    for (ii = 1; ii < N; ii++) {
        if (hq[ii] == hq[ii - 1] && hr[ii] == hr[ii - 1] &&
            hQ[ii] == hQ[ii - 1] && hS[ii] == hS[ii - 1] &&
            hR[ii] == hR[ii - 1]) {
            //            printf("\nsame\n");
            hpRSQrq[ii] = hpRSQrq[ii - 1];
        } else {
            //            printf("\ndifferent\n");
            hpRSQrq[ii] = ptr;
            ptr += pnz[ii] * cnux[ii];
            d_cvt_mat2pmat(nu[ii], nu[ii], hR[ii], nu[ii], 0, hpRSQrq[ii],
                           cnux[ii]);
            d_cvt_tran_mat2pmat(
                nx[ii], nu[ii], hS[ii], nu[ii], nu[ii],
                hpRSQrq[ii] + nu[ii] / bs * bs * cnux[ii] + nu[ii] % bs,
                cnux[ii]);
            d_cvt_mat2pmat(nx[ii], nx[ii], hQ[ii], nx[ii], nu[ii],
                           hpRSQrq[ii] + nu[ii] / bs * bs * cnux[ii] +
                               nu[ii] % bs + nu[ii] * bs,
                           cnux[ii]);
            d_cvt_tran_mat2pmat(nu[ii], 1, hr[ii], nu[ii], nu[ii] + nx[ii],
                                hpRSQrq[ii] +
                                    (nu[ii] + nx[ii]) / bs * bs * cnux[ii] +
                                    (nu[ii] + nx[ii]) % bs,
                                cnux[ii]);
            d_cvt_tran_mat2pmat(nx[ii], 1, hq[ii], nx[ii], nu[ii] + nx[ii],
                                hpRSQrq[ii] +
                                    (nu[ii] + nx[ii]) / bs * bs * cnux[ii] +
                                    (nu[ii] + nx[ii]) % bs + nu[ii] * bs,
                                cnux[ii]);
        }
        //        d_print_pmat(nu[ii]+nx[ii]+1, nu[ii]+nx[ii], bs, hpRSQrq[ii],
        //        cnux[ii]);
    }
    ii = N;
    hpRSQrq[ii] = ptr;
    ptr += pnz[ii] * cnux[ii];
    d_cvt_mat2pmat(nx[ii], nx[ii], hQ[ii], nx[ii], 0, hpRSQrq[ii], cnux[ii]);
    d_cvt_tran_mat2pmat(
        nx[ii], 1, hq[ii], nx[ii], nx[ii],
        hpRSQrq[ii] + (nx[ii]) / bs * bs * cnux[ii] + (nx[ii]) % bs, cnux[ii]);
    //    d_print_pmat(nu[ii]+nx[ii]+1, nu[ii]+nx[ii], bs, hpRSQrq[ii],
    //    cnux[ii]);

    ii = 0;
    hpDCt[ii] = ptr;
    ptr += pnz[ii] * cng[ii];
    d_cvt_tran_mat2pmat(ng[ii], nu[ii], hD[ii], ng[ii], 0, hpDCt[ii], cng[ii]);
    d_cvt_tran_mat2pmat(ng[ii], nx[ii], hC[ii], ng[ii], nu[ii],
                        hpDCt[ii] + nu[ii] / bs * bs * cng[ii] + nu[ii] % bs,
                        cng[ii]);
    //    d_print_pmat(nu[ii]+nx[ii], ng[ii], bs, hpDCt[ii], cng[ii]);
    for (ii = 1; ii < N; ii++) {
        if (hD[ii] == hD[ii - 1] && hC[ii] == hC[ii - 1]) {
            //            printf("\nsame\n");
            hpDCt[ii] = hpDCt[ii - 1];
        } else {
            //            printf("\ndifferent\n");
            hpDCt[ii] = ptr;
            ptr += pnux[ii] * cng[ii];
            d_cvt_tran_mat2pmat(ng[ii], nu[ii], hD[ii], ng[ii], 0, hpDCt[ii],
                                cng[ii]);
            d_cvt_tran_mat2pmat(
                ng[ii], nx[ii], hC[ii], ng[ii], nu[ii],
                hpDCt[ii] + nu[ii] / bs * bs * cng[ii] + nu[ii] % bs, cng[ii]);
        }
        //        d_print_pmat(nu[ii]+nx[ii], ng[ii], bs, hpDCt[ii], cng[ii]);
    }
    ii = N;
    hpDCt[ii] = ptr;
    ptr += pnux[ii] * cng[ii];
    d_cvt_tran_mat2pmat(ng[ii], nx[ii], hC[ii], ng[ii], 0, hpDCt[ii], cng[ii]);
    //    d_print_pmat(nu[ii]+nx[ii], ng[ii], bs, hpDCt[ii], cng[ii]);

    ii = 0;
    hd[ii] = ptr;
    ptr += 2 * pnb[ii] + 2 * png[ii];
    d_copy_mat(nb[ii], 1, hlb[ii], 1, hd[ii], 1);
    dax_mat(nb[ii], 1, -1.0, hub[ii], 1, hd[ii] + pnb[ii],
            1);  // TODO(giaf): change in solver
    d_copy_mat(ng[ii], 1, hlg[ii], 1, hd[ii] + 2 * pnb[ii], 1);
    dax_mat(ng[ii], 1, -1.0, hug[ii], 1, hd[ii] + 2 * pnb[ii] + png[ii],
            1);  // TODO(giaf): change in solver
                 //    d_print_mat(1, 2*pnb[ii]+2*png[ii], hd[ii], 1);
    for (ii = 1; ii < N; ii++) {
        if (hlb[ii] == hlb[ii - 1] && hub[ii] == hub[ii - 1] &&
            hlg[ii] == hlg[ii - 1] && hug[ii] == hug[ii - 1]) {
            //            printf("\nsame\n");
            hd[ii] = hd[ii - 1];
        } else {
            //            printf("\ndifferent\n");
            hd[ii] = ptr;
            ptr += 2 * pnb[ii] + 2 * png[ii];
            d_copy_mat(nb[ii], 1, hlb[ii], 1, hd[ii], 1);
            dax_mat(nb[ii], 1, -1.0, hub[ii], 1, hd[ii] + pnb[ii],
                    1);  // TODO(giaf): change in solver
            d_copy_mat(ng[ii], 1, hlg[ii], 1, hd[ii] + 2 * pnb[ii], 1);
            dax_mat(ng[ii], 1, -1.0, hug[ii], 1, hd[ii] + 2 * pnb[ii] + png[ii],
                    1);  // TODO(giaf): change in solver
        }
        //        d_print_mat(1, 2*pnb[ii]+2*png[ii], hd[ii], 1);
    }
    ii = N;
    hd[ii] = ptr;
    ptr += 2 * pnb[ii] + 2 * png[ii];
    d_copy_mat(nb[ii], 1, hlb[ii], 1, hd[ii], 1);
    dax_mat(nb[ii], 1, -1.0, hub[ii], 1, hd[ii] + pnb[ii],
            1);  // TODO(giaf): change in solver
    d_copy_mat(ng[ii], 1, hlg[ii], 1, hd[ii] + 2 * pnb[ii], 1);
    dax_mat(ng[ii], 1, -1.0, hug[ii], 1, hd[ii] + 2 * pnb[ii] + png[ii],
            1);  // TODO(giaf): change in solver
                 //    d_print_mat(1, 2*pnb[ii]+2*png[ii], hd[ii], 1);

    for (ii = 0; ii <= N; ii++) {
        hux[ii] = ptr;
        ptr += pnux[ii];
    }

    for (ii = 0; ii <= N; ii++) {
        hpi[ii] = ptr;
        ptr += pnx[ii];
    }

    for (ii = 0; ii <= N; ii++) {
        ht[ii] = ptr;
        hlam[ii] = ptr + 2 * pnb[ii] + 2 * png[ii];
        ptr += 4 * pnb[ii] + 4 * png[ii];
    }

    stat = ptr;
    ptr += 5 * k_max;

    // align work space (again) (assume cache line <= 64 byte)
    addr = (((size_t)ptr) + 63) / 64 * 64;
    // size_t offset = addr % 64;
    // double *ptr = work0 + offset / 8;
    work = (double *)addr;

    // estimate mu0 if not user-provided
    if (mu0 <= 0) {
        // first stage
        ii = 0;
        for (jj = 0; jj < nu[ii]; jj++)
            for (ll = 0; ll < nu[ii]; ll++)
                mu0 = fmax(mu0, fabs(hR[ii][jj * nu[ii] + ll]));
        for (jj = 0; jj < nx[ii]; jj++)
            for (ll = 0; ll < nu[ii]; ll++)
                mu0 = fmax(mu0, fabs(hS[ii][jj * nu[ii] + ll]));
        for (jj = 0; jj < nx[ii]; jj++)
            for (ll = 0; ll < nx[ii]; ll++)
                mu0 = fmax(mu0, fabs(hQ[ii][jj * nx[ii] + ll]));
        for (jj = 0; jj < nu[ii]; jj++) mu0 = fmax(mu0, fabs(hr[ii][jj]));
        for (jj = 0; jj < nx[ii]; jj++) mu0 = fmax(mu0, fabs(hq[ii][jj]));
        // middle stages
        for (jj = 1; jj < N; jj++) {
            if (hq[ii] == hq[ii - 1] && hr[ii] == hr[ii - 1] &&
                hQ[ii] == hQ[ii - 1] && hS[ii] == hS[ii - 1] &&
                hR[ii] == hR[ii - 1]) {
                for (jj = 0; jj < nu[ii]; jj++)
                    for (ll = 0; ll < nu[ii]; ll++)
                        mu0 = fmax(mu0, fabs(hR[ii][jj * nu[ii] + ll]));
                for (jj = 0; jj < nx[ii]; jj++)
                    for (ll = 0; ll < nu[ii]; ll++)
                        mu0 = fmax(mu0, fabs(hS[ii][jj * nu[ii] + ll]));
                for (jj = 0; jj < nx[ii]; jj++)
                    for (ll = 0; ll < nx[ii]; ll++)
                        mu0 = fmax(mu0, fabs(hQ[ii][jj * nx[ii] + ll]));
                for (jj = 0; jj < nu[ii]; jj++)
                    mu0 = fmax(mu0, fabs(hr[ii][jj]));
                for (jj = 0; jj < nx[ii]; jj++)
                    mu0 = fmax(mu0, fabs(hq[ii][jj]));
            }
        }
        // last stage
        ii = N;
        for (jj = 0; jj < nx[ii]; jj++)
            for (ll = 0; ll < nx[ii]; ll++)
                mu0 = fmax(mu0, fabs(hQ[ii][jj * nx[ii] + ll]));
        for (jj = 0; jj < nx[ii]; jj++) mu0 = fmax(mu0, fabs(hq[ii][jj]));
    }

    // TODO(giaf): check for equality constraints in the inputs

    // call the solver
    int hpmpc_status = d_ip2_hard_mpc_tv(
        &kk, k_max, mu0, mu_tol, alpha_min, warm_start, sigma_par, stat, N, nx,
        nu, nb, hidxb, ng, hpBAbt, hpRSQrq, hpDCt, hd, hux, compute_mult, hpi,
        hlam, ht, work);

    if (hpmpc_status == 1) acados_status = ACADOS_MAXITER;

    if (hpmpc_status == 2) acados_status = ACADOS_MINSTEP;

    // copy back inputs and states
    for (ii = 0; ii < N; ii++)
        for (jj = 0; jj < nu[ii]; jj++) hu[ii][jj] = hux[ii][jj];

    for (ii = 0; ii <= N; ii++)
        for (jj = 0; jj < nx[ii]; jj++) hx[ii][jj] = hux[ii][nu[ii] + jj];

#if 0
    printf("\n");
    for (ii = 0; ii < kk; ii++)
        printf("%d %e %e %e %e %e\n", ii, stat[0+ii*5], stat[1+ii*5],
               stat[2+ii*5], stat[3+ii*5], stat[4+ii*5]);
#endif

    // TODO(giaf): check for equality constraints in the inputs

    // TODO(giaf): compute residuals ?????

    // return

    return acados_status;
}
