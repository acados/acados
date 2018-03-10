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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// flush denormals to zero
#if defined(TARGET_X64_AVX2) || defined(TARGET_X64_AVX) ||  \
    defined(TARGET_X64_SSE3) || defined(TARGET_X86_ATOM) || \
    defined(TARGET_AMD_SSE3)
#include <xmmintrin.h>  // needed to flush to zero sub-normals with _MM_SET_FLUSH_ZERO_MODE (_MM_FLUSH_ZERO_ON); in the main()
#endif

#include "hpmpc/include/aux_d.h"
#include "hpmpc/include/mpc_solvers.h"

#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/ocp_qp/ocp_qp_hpmpc.h"
#include "acados/utils/timing.h"
#include "acados/utils/math.h"
#include "acados/utils/types.h"

// define IP solver arguments && number of repetitions
#define NREP 1
#define MAXITER 50
#define TOL 1e-8
#define MINSTEP 1e-8

/************************************************
Mass-spring system: nx/2 masses connected each other with springs (in a row),
and the first and the last one to walls. nu (<=nx) controls act on the first nu
masses. The system is sampled with sampling time Ts.
************************************************/
void mass_spring_system(double Ts, int nx, int nu, double *A, double *B,
                        double *b, double *x0) {
    int nx2 = nx * nx;

    int info = 0;

    int pp = nx / 2;  // number of masses

    /************************************************
     * build the continuous time system
     ************************************************/

    double *T;
    d_zeros(&T, pp, pp);
    int ii;
    for (ii = 0; ii < pp; ii++) T[ii * (pp + 1)] = -2;
    for (ii = 0; ii < pp - 1; ii++) T[ii * (pp + 1) + 1] = 1;
    for (ii = 1; ii < pp; ii++) T[ii * (pp + 1) - 1] = 1;

    double *Z;
    d_zeros(&Z, pp, pp);
    double *I;
    d_zeros(&I, pp, pp);
    for (ii = 0; ii < pp; ii++) I[ii * (pp + 1)] = 1.0;  // = eye(pp);
    double *Ac;
    d_zeros(&Ac, nx, nx);
    dmcopy(pp, pp, Z, pp, Ac, nx);
    dmcopy(pp, pp, T, pp, Ac + pp, nx);
    dmcopy(pp, pp, I, pp, Ac + pp * nx, nx);
    dmcopy(pp, pp, Z, pp, Ac + pp * (nx + 1), nx);
    free(T);
    free(Z);
    free(I);

    d_zeros(&I, nu, nu);
    for (ii = 0; ii < nu; ii++) I[ii * (nu + 1)] = 1.0;  // I = eye(nu);
    double *Bc;
    d_zeros(&Bc, nx, nu);
    dmcopy(nu, nu, I, nu, Bc + pp, nx);
    free(I);

    /************************************************
     * compute the discrete time system
     ************************************************/

    double *bb;
    d_zeros(&bb, nx, 1);
    dmcopy(nx, 1, bb, nx, b, nx);

    dmcopy(nx, nx, Ac, nx, A, nx);
    dscal_3l(nx2, Ts, A);
    expm(nx, A);

    d_zeros(&T, nx, nx);
    d_zeros(&I, nx, nx);
    for (ii = 0; ii < nx; ii++) I[ii * (nx + 1)] = 1.0;  // I = eye(nx);
    dmcopy(nx, nx, A, nx, T, nx);
    daxpy_3l(nx2, -1.0, I, T);
    dgemm_nn_3l(nx, nu, nx, T, nx, Bc, nx, B, nx);

    int *ipiv = (int *)malloc(nx * sizeof(int));
    dgesv_3l(nx, nu, Ac, nx, ipiv, B, nx, &info);
    free(ipiv);

    free(Ac);
    free(Bc);
    free(bb);

    /************************************************
     * initial state
     ************************************************/

    if (nx == 4) {
        x0[0] = 5;
        x0[1] = 10;
        x0[2] = 15;
        x0[3] = 20;
    } else {
        int jj;
        for (jj = 0; jj < nx; jj++) x0[jj] = 1;
    }
}

int main() {
    printf("\n");
    printf("\n");
    printf("\n");
    printf(
        " HPMPC -- Library for High-Performance implementation of solvers for "
        "MPC.\n");
    printf(
        " Copyright (C) 2014-2015 by Technical University of Denmark. All "
        "rights reserved.\n");
    printf("\n");
    printf(" HPMPC is distributed in the hope that it will be useful,\n");
    printf(" but WITHOUT ANY WARRANTY; without even the implied warranty of\n");
    printf(" MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n");
    printf(" See the GNU Lesser General Public License for more details.\n");
    printf("\n");
    printf("\n");
    printf("\n");

#if defined(TARGET_X64_AVX2) || defined(TARGET_X64_AVX) ||  \
    defined(TARGET_X64_SSE3) || defined(TARGET_X86_ATOM) || \
    defined(TARGET_AMD_SSE3)
    _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);  // flush to zero subnormals !!!
                                                 // works only with one thread
                                                 // !!!
#endif

    int ii, jj;

    int rep, nrep = NREP;

    int nx = 8;  // number of states (it has to be even for the mass-spring
                 // system test problem)
    int nu = 3;  // number of inputs (controllers) (it has to be at least 1 and
                 // at most nx/2 for the mass-spring system test problem)
    int N = 20;  // horizon length
    int M = 2;
    // int nb = 11;  // number of box constrained inputs and states
    int ng = 0;   // 4;  // number of general constraints
    int ngN = 0;  // 4;  // number of general constraints at the last stage

    // int nbu = nu < nb ? nu : nb;
    // int nbx = nb - nu > 0 ? nb - nu : 0;
    int nbu = 3;
    int nbx = 4;
    // stage-wise variant size
    int nxx[N + 1];
    nxx[0] = 0;
    for (ii = 1; ii <= N; ii++) nxx[ii] = nx;

    int nuu[N + 1];
    for (ii = 0; ii < N; ii++) nuu[ii] = nu;
    nuu[N] = 0;

    int nbb[N + 1];
    // nbb[0] = nbu;
    // for (ii = 1; ii < N; ii++) nbb[ii] = nb;
    // nbb[N] = nbx;

    // Andrea XXX change this back to 11 bounds, changed to debug the strmat
    // interface
    for (ii = 0; ii < M; ii++)  // XXX not M !!!
        nbb[ii] = nuu[ii] + nxx[ii] / 2;

    for (; ii <= N; ii++) nbb[ii] = 0;

    int ngg[N + 1];
    for (ii = 0; ii < N; ii++) ngg[ii] = ng;
    ngg[N] = ngN;

    printf(
        " Test problem: mass-spring system with %d masses and %d controls.\n",
        nx / 2, nu);
    printf("\n");
    printf(
        " MPC problem size: %d states, %d inputs, %d horizon length, %d "
        "two-sided box constraints, %d two-sided general constraints.\n",
        nx, nu, N, 7, ng);
    printf("\n");
    printf(
        " IP method parameters: predictor-corrector IP, double precision, %d "
        "maximum iterations, %5.1e exit tolerance in duality measure.\n",
        MAXITER, TOL);
    printf("\n");
#if defined(TARGET_X64_AVX2)
    printf(" HPMPC built for the AVX2 architecture\n");
#endif
#if defined(TARGET_X64_AVX)
    printf(" HPMPC built for the AVX architecture\n");
#endif
    printf("\n");

    /************************************************
     * dynamical system
     ************************************************/

    // state space matrices & initial state
    double *A;
    d_zeros(&A, nx, nx);  // states update matrix
    double *B;
    d_zeros(&B, nx, nu);  // inputs matrix
    double *b;
    d_zeros(&b, nx, 1);  // states offset
    double *x0;
    d_zeros(&x0, nx, 1);  // initial state

    // mass-spring system
    double Ts = 0.5;  // sampling time
    mass_spring_system(Ts, nx, nu, A, B, b, x0);

    for (jj = 0; jj < nx; jj++) b[jj] = 0.1;

    for (jj = 0; jj < nx; jj++) x0[jj] = 0;
    x0[0] = 2.5;
    x0[1] = 2.5;

    //    d_print_mat(nx, nx, A, nx);
    //    d_print_mat(nx, nu, B, nx);
    //    d_print_mat(nx, 1, b, nx);
    //    d_print_mat(nx, 1, x0, nx);

    // compute b0 = b + A*x0
    double *b0;
    d_zeros(&b0, nx, 1);
    dcopy_3l(nx, b, 1, b0, 1);
    dgemv_n_3l(nx, nx, A, nx, x0, b0);
    //    d_print_mat(nx, 1, b, nx);
    //    d_print_mat(nx, 1, b0, nx);

    // then A0 is a matrix of size 0x0
    double *A0;
    d_zeros(&A0, 0, 0);

    /************************************************
     * general constraints
     ************************************************/

    double *C;
    d_zeros(&C, ng, nx);
    double *D;
    d_zeros(&D, ng, nu);
    double *lg;
    d_zeros(&lg, ng, 1);
    double *ug;
    d_zeros(&ug, ng, 1);

    double *CN;
    d_zeros(&CN, ngN, nx);
    // for (ii = 0; ii < ngN; ii++) CN[ii * (ngN + 1)] = 1.0;
    //    d_print_mat(ngN, nx, CN, ngN);
    double *lgN;
    d_zeros(&lgN, ngN, 1);  // force all states to 0 at the last stage
    double *ugN;
    d_zeros(&ugN, ngN, 1);  // force all states to 0 at the last stage

    /************************************************
     * box & general constraints
     ************************************************/

    int *idxb0;
    int_zeros(&idxb0, nbb[0], 1);
    // double *d0; d_zeros(&d0, 2*nb[0]+2*ng[0], 1);
    double *lb0;
    d_zeros(&lb0, nbb[1], 1);
    double *ub0;
    d_zeros(&ub0, nbb[1], 1);
    for (ii = 0; ii < nbb[0]; ii++) {
        if (ii < nuu[0]) {
            lb0[ii] = -0.5;  // umin
            ub0[ii] = 0.5;   // umax
        } else {
            lb0[ii] = -4.0;  // xmin
            ub0[ii] = 4.0;   // xmax
        }
        idxb0[ii] = ii;
    }

    // for(ii=0; ii<ng; ii++)  //Andrea: no general constraints atm
    // {
    // // d0[2*nb[0]+ii]       = - 100.0; // dmin
    // // d0[2*nb[0]+ng[0]+ii] =   100.0; // dmax
    // }
    // i_print_mat(1, nb[0], idxb0, 1);
    // d_print_mat(1, 2*nb[0]+2*ng[0], d0, 1);

    int *idxb1;
    int_zeros(&idxb1, nbb[1], 1);
    // double *d1; d_zeros(&d1, 2*nb[1]+2*ng[1], 1);
    int_zeros(&idxb1, nbb[1], 1);
    double *lb1;
    d_zeros(&lb1, nbb[1], 1);
    double *ub1;
    d_zeros(&ub1, nbb[1], 1);
    for (ii = 0; ii < nbb[1]; ii++) {
        if (ii < nuu[1]) {   // input
            lb1[ii] = -0.5;  // umin
            ub1[ii] = 0.5;   // umax
        } else {             // state
            lb1[ii] = -4.0;  // xmin
            ub1[ii] = 4.0;   // xmax
        }
        idxb1[ii] = ii;
    }

    // for(ii=0; ii<ng[1]; ii++)  //Andrea: no general constraints atm
    // {
    // // d1[2*nb[1]+ii]       = - 100.0; // dmin
    // // d1[2*nb[1]+ng[1]+ii] =   100.0; // dmax
    // }
    // i_print_mat(1, nb[1], idxb1, 1);
    // d_print_mat(1, 2*nb[1]+2*ng[1], d1, 1);

    int *idxbN;
    int_zeros(&idxbN, nbb[N], 1);
    // double *dN; d_zeros(&dN, 2*nb[N]+2*ng[N], 1);
    int_zeros(&idxbN, nbb[N], 1);
    double *lbN;
    d_zeros(&lbN, nbb[N], 1);
    double *ubN;
    d_zeros(&ubN, nbb[N], 1);
    for (ii = 0; ii < nbb[N]; ii++) {
        if (ii < nuu[N]) {
            lbN[ii] = -0.5;  // umin
            ubN[ii] = 0.5;   // umax
        } else {
            lbN[ii] = -4.0;  // xmin
            ubN[ii] = 4.0;   // xmax
        }
        idxbN[ii] = ii;
    }
    // for(ii=0; ii<ng[N]; ii++)//Andrea: no general constraints atm
    // {
    // // dN[2*nb[N]+ii]       = - 100.0; // dmin
    // // dN[2*nb[N]+ng[N]+ii] =   100.0; // dmax
    // }
    // i_print_mat(1, nb[N], idxbN, 1);
    // d_print_mat(1, 2*nb[N]+2*ng[N], dN, 1);

    /************************************************
     * cost function
     ************************************************/

    double *Q;
    d_zeros(&Q, nx, nx);
    for (ii = 0; ii < nx; ii++) Q[ii * (nx + 1)] = 1.0;

    double *R;
    d_zeros(&R, nu, nu);
    for (ii = 0; ii < nu; ii++) R[ii * (nu + 1)] = 2.0;

    double *S;
    d_zeros(&S, nu, nx);

    double *q;
    d_zeros(&q, nx, 1);
    for (ii = 0; ii < nx; ii++) q[ii] = 0.1;

    double *r;
    d_zeros(&r, nu, 1);
    for (ii = 0; ii < nu; ii++) r[ii] = 0.2;

    // Q0 and q0 are matrices of size 0
    double *Q0;
    d_zeros(&Q0, 0, 0);
    double *q0;
    d_zeros(&q0, 0, 1);

    // compute r0 = r + S*x0
    double *r0;
    d_zeros(&r0, nu, 1);
    dcopy_3l(nu, r, 1, r0, 1);
    dgemv_n_3l(nu, nx, S, nu, x0, r0);

    // then S0 is a matrix of size nux0
    double *S0;
    d_zeros(&S0, nu, 0);

    /************************************************
     * problems data
     ************************************************/

    double *hA[N];
    double *hB[N];
    double *hb[N];
    double *hQ[N + 1];
    double *hS[N];
    double *hR[N];
    double *hq[N + 1];
    double *hr[N];
    double *hlb[N + 1];
    double *hub[N + 1];
    int *hidxb[N + 1];
    double *hC[N + 1];
    double *hD[N];
    double *hlg[N + 1];
    double *hug[N + 1];

    hA[0] = A0;
    hB[0] = B;
    hb[0] = b0;
    hQ[0] = Q0;
    hS[0] = S0;
    hR[0] = R;
    hq[0] = q0;
    hr[0] = r0;
    hlb[0] = lb0;
    hub[0] = ub0;
    hidxb[0] = idxb0;
    hC[0] = C;
    hD[0] = D;
    hlg[0] = lg;
    hug[0] = ug;
    for (ii = 1; ii < N; ii++) {
        hA[ii] = A;
        hB[ii] = B;
        hb[ii] = b;
        hQ[ii] = Q;
        hS[ii] = S;
        hR[ii] = R;
        hq[ii] = q;
        hr[ii] = r;
        hlb[ii] = lb1;
        hub[ii] = ub1;
        hidxb[ii] = idxb1;
        hC[ii] = C;
        hD[ii] = D;
        hlg[ii] = lg;
        hug[ii] = ug;
    }
    hQ[N] = Q;  // or maybe initialize to the solution of the DARE???
    hq[N] = q;  // or maybe initialize to the solution of the DARE???
    hlb[N] = lbN;
    hub[N] = ubN;
    hidxb[N] = idxbN;
    hC[N] = CN;
    hlg[N] = lgN;
    hug[N] = ugN;

    /************************************************
     * solution
     ************************************************/

    double *hx[N + 1];
    double *hu[N];
    double *hpi[N];
    double *hlam[N + 1];
    double *ht[N + 1];

    double *lam_in[N + 1];
    double *t_in[N + 1];

    for (ii = 0; ii < N; ii++) {
        d_zeros(&hx[ii], nxx[ii], 1);
        d_zeros(&hu[ii], nuu[ii], 1);
        d_zeros(&hpi[ii], nxx[ii + 1], 1);
        d_zeros(&hlam[ii], 2 * nbb[ii] + 2 * ngg[ii], 1);
        d_zeros(&ht[ii], 2 * nbb[ii] + 2 * ngg[ii], 1);
        d_zeros(&lam_in[ii], 2 * nbb[ii] + 2 * ngg[ii], 1);
        d_zeros(&t_in[ii], 2 * nbb[ii] + 2 * ngg[ii], 1);

        // Init multipliers and slacks
        for (jj = 0; jj < 2 * nbb[ii] + 2 * ngg[ii]; jj++) {
            lam_in[ii][jj] = 1.0;
            t_in[ii][jj] = 1.0;
        }
    }
    d_zeros(&hx[N], nxx[N], 1);
    d_zeros(&hlam[N], 2 * nbb[N] + 2 * ngg[N], 1);
    d_zeros(&ht[N], 2 * nbb[N] + 2 * ngg[N], 1);
    d_zeros(&lam_in[N], 2 * nbb[N] + 2 * ngg[N], 1);
    d_zeros(&t_in[N], 2 * nbb[N] + 2 * ngg[N], 1);

    // Init multipliers and slacks
    for (jj = 0; jj < 2 * nbb[ii] + 2 * ngg[ii]; jj++) {
        lam_in[N][jj] = 1.0;
        t_in[N][jj] = 1.0;
    }

    /************************************************
     * Solver arguments
     ************************************************/

    // solver arguments
    ocp_qp_hpmpc_opts hpmpc_args;
    hpmpc_args.tol = TOL;
    hpmpc_args.max_iter = MAXITER;
    //  hpmpc_args.min_step = MINSTEP;
    hpmpc_args.mu0 = 0.1;
    //  hpmpc_args.sigma_min = 1e-3;
    hpmpc_args.warm_start = 0;
    hpmpc_args.N2 = N;
    hpmpc_args.lam0 = lam_in;
    hpmpc_args.t0 = t_in;

    /************************************************
     * work space
     ************************************************/

    int work_space_size =
        d_ip2_res_mpc_hard_work_space_size_bytes_libstr(N, nxx, nuu, nbb, ngg);

    // Adding memory for data
    for (int ii = 0; ii <= N; ii++) {
        work_space_size += blasfeo_memsize_dmat(nuu[ii] + nxx[ii] + 1, nxx[ii + 1]);
        work_space_size += blasfeo_memsize_dvec(nxx[ii + 1]);
        work_space_size +=
            blasfeo_memsize_dmat(nuu[ii] + nxx[ii] + 1, nuu[ii] + nxx[ii]);
        work_space_size += blasfeo_memsize_dvec(nuu[ii] + nxx[ii]);
        work_space_size += blasfeo_memsize_dmat(nuu[ii] + nxx[ii] + 1, ngg[ii]);
        work_space_size += blasfeo_memsize_dvec(2 * nbb[ii] + 2 * ngg[ii]);
        work_space_size += blasfeo_memsize_dvec(nuu[ii] + nxx[ii]);
        work_space_size += blasfeo_memsize_dvec(nxx[ii + 1]);
        work_space_size += blasfeo_memsize_dvec(2 * nbb[ii] + 2 * ngg[ii]);
        work_space_size += blasfeo_memsize_dvec(2 * nbb[ii] + 2 * ngg[ii]);
        work_space_size += blasfeo_memsize_dvec(2 * nbb[ii] + 2 * ngg[ii]);
        work_space_size += blasfeo_memsize_dvec(2 * nbb[ii] + 2 * ngg[ii]);
    }

    work_space_size += 10000 * sizeof(double) * (N + 1);
    void *workspace;

    v_zeros_align(&workspace, work_space_size);

    /************************************************
     * create the in and out struct
     ************************************************/

    ocp_qp_in qp_in;
    qp_in.N = N;
    qp_in.nx = (const int *)nxx;
    qp_in.nu = (const int *)nuu;
    qp_in.nb = (const int *)nbb;
    qp_in.nc = (const int *)ngg;
    qp_in.A = (const double **)hA;
    qp_in.B = (const double **)hB;
    qp_in.b = (const double **)hb;
    qp_in.Q = (const double **)hQ;
    qp_in.S = (const double **)hS;
    qp_in.R = (const double **)hR;
    qp_in.q = (const double **)hq;
    qp_in.r = (const double **)hr;
    qp_in.idxb = (const int **)hidxb;
    qp_in.lb = (const double **)hlb;
    qp_in.ub = (const double **)hub;
    qp_in.Cx = (const double **)hC;
    qp_in.Cu = (const double **)hD;
    qp_in.lc = (const double **)hlg;
    qp_in.uc = (const double **)hug;

    ocp_qp_out qp_out;
    qp_out.x = hx;
    qp_out.u = hu;
    qp_out.pi = hpi;
    qp_out.lam = hlam;

    /************************************************
     * call the solver
     ************************************************/

    int return_value;

    acados_timer tv0;
    acados_tic(&tv0);  // start

    // call the QP OCP solver
    // return_value = ocp_qp_hpmpc_libstr(&qp_in, &qp_out, &hpmpc_args,
    // workspace);
    return_value =
        ocp_qp_hpmpc_libstr_pt(&qp_in, &qp_out, &hpmpc_args, M, 0.1, workspace);

    double time = acados_toc(&tv0);  // stop

    if (return_value == ACADOS_SUCCESS)
        printf("\nACADOS status: solution found\n");

    if (return_value == ACADOS_MAXITER)
        printf("\nACADOS status: maximum number of iterations reached\n");

    if (return_value == ACADOS_MINSTEP)
        printf("\nACADOS status: below minimum step size length\n");

    printf("\nu = \n");
    for (ii = 0; ii < N; ii++) d_print_mat(1, nuu[ii], hu[ii], 1);

    printf("\nx = \n");
    for (ii = 0; ii <= N; ii++) d_print_mat(1, nxx[ii], hx[ii], 1);

    printf("\n");
    printf(" Average solution time over %d runs: %5.2e seconds\n", nrep, time);
    printf("\n\n");

    /************************************************
     * free memory
     ************************************************/

    // d_free(A);
    // d_free(B);
    // d_free(b);
    // d_free(x0);
    // d_free(A0);
    // d_free(b0);
    // d_free(Q);
    // d_free(S);
    // d_free(R);
    // d_free(q);
    // d_free(r);
    // d_free(Q0);
    // d_free(S0);
    // d_free(q0);
    // d_free(r0);
    // i_free(idxb0);
    // d_free(lb0);
    // d_free(ub0);
    // i_free(idxb1);
    // d_free(lb1);
    // d_free(ub1);
    // i_free(idxbN);
    // d_free(lbN);
    // d_free(ubN);
    // d_free(C);
    // d_free(D);
    // d_free(lg);
    // d_free(ug);
    // d_free(CN);
    // d_free(lgN);
    // d_free(ugN);

    // for (ii = 0; ii < N; ii++) {
    //     d_free(hx[ii]);
    //     d_free(hu[ii]);
    // }
    // d_free(hx[N]);

    // free(workspace);

    return 0;
}
