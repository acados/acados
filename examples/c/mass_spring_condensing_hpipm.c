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
#include <sys/time.h>

#include "blasfeo/include/blasfeo_target.h"
#include "blasfeo/include/blasfeo_d_aux_ext_dep.h"
#include "blasfeo/include/blasfeo_i_aux_ext_dep.h"

// flush denormals to zero
#if defined(TARGET_X64_INTEL_HASWELL) ||      \
    defined(TARGET_X64_INTEL_SABDY_BRIDGE) || \
    defined(TARGET_X64_INTEL_CORE) || defined(TARGET_X86_AMD_BULLDOZER)
#include <xmmintrin.h>  // needed to flush to zero sub-normals with _MM_SET_FLUSH_ZERO_MODE (_MM_FLUSH_ZERO_ON); in the main()
#endif

#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/ocp_qp/ocp_qp_condensing_hpipm.h"
#include "acados/utils/math.h"
#include "acados/utils/types.h"
#include "acados/utils/timing.h"

// define IP solver arguments && number of repetitions
#define NREP 1000
#define MAXITER 20
#define TOL 1e-8
#define MINSTEP 1e-8

#define ELIMINATE_X0

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

#if defined(TARGET_X64_INTEL_HASWELL) ||      \
    defined(TARGET_X64_INTEL_SABDY_BRIDGE) || \
    defined(TARGET_X64_INTEL_CORE) || defined(TARGET_X86_AMD_BULLDOZER)
    _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);  // flush to zero subnormals !!!
                                                 // works only with one thread
                                                 // !!!
#endif

    int ii, jj;

    int rep, nrep = NREP;

    int nx = 8;   // number of states (it has to be even for the mass-spring
                  // system test problem)
    int nu = 3;   // number of inputs (controllers) (it has to be at least 1 and
                  // at most nx/2 for the mass-spring system test problem)
    int N = 15;   // horizon length
    int nb = 11;  // number of box constrained inputs and states
    int ng = 0;   // 4;  // number of general constraints
    int ngN = 4;  // 4;  // number of general constraints at the last stage

    //    int N2 = 3;   // horizon length of partially condensed problem

    int nbu = nu < nb ? nu : nb;
    int nbx = nb - nu > 0 ? nb - nu : 0;

    // stage-wise variant size
    int nxx[N + 1];
#if defined(ELIMINATE_X0)
    nxx[0] = 0;
#else
    nxx[0] = nx;
#endif
    for (ii = 1; ii <= N; ii++) nxx[ii] = nx;

    int nuu[N + 1];
    for (ii = 0; ii < N; ii++) nuu[ii] = nu;
    nuu[N] = 0;

    int nbb[N + 1];
#if defined(ELIMINATE_X0)
    nbb[0] = nbu;
#else
    nbb[0] = nb;
#endif
    for (ii = 1; ii < N; ii++) nbb[ii] = nb;
    nbb[N] = nbx;

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
        nx, nu, N, nb, ng);
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

#if defined(ELIMINATE_X0)
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
#endif

    /************************************************
     * box constraints
     ************************************************/

    int jj_end;

    int *idxb0;
    int_zeros(&idxb0, nbb[0], 1);
    double *lb0;
    d_zeros(&lb0, nbb[0], 1);
    double *ub0;
    d_zeros(&ub0, nbb[0], 1);
#if defined(ELIMINATE_X0)
    for (jj = 0; jj < nbb[0]; jj++) {
        lb0[jj] = -0.5;  // umin
        ub0[jj] = +0.5;  // umin
        idxb0[jj] = jj;
    }
#else
    jj_end = nbx < nbb[0] ? nbx : nbb[0];
    for (jj = 0; jj < jj_end; jj++) {
//        lb0[jj] = x0[jj - nbu];  // initial state
//        ub0[jj] = x0[jj - nbu];  // initial state
        lb0[jj] = x0[jj];  // initial state
        ub0[jj] = x0[jj];  // initial state
        idxb0[jj] = jj;
    }
    for (; jj < nbb[0]; jj++) {
        lb0[jj] = -0.5;  // umin
        ub0[jj] = +0.5;  // umax
        idxb0[jj] = jj;
    }
#endif
    //    int_print_mat(nbb[0], 1, idxb0, nbb[0]);
    //    d_print_mat(nbb[0], 1, lb0, nbb[0]);

    int *idxb1;
    int_zeros(&idxb1, nbb[1], 1);
    double *lb1;
    d_zeros(&lb1, nbb[1], 1);
    double *ub1;
    d_zeros(&ub1, nbb[1], 1);
    jj_end = nbx < nbb[1] ? nbx : nbb[1];
    for (jj = 0; jj < jj_end; jj++) {
        lb1[jj] = -4.0;  // xmin
        ub1[jj] = +4.0;  // xmax
        idxb1[jj] = jj;
    }
    for (; jj < nbb[1]; jj++) {
        lb1[jj] = -0.5;  // umin
        ub1[jj] = +0.5;  // umax
        idxb1[jj] = jj;
    }
    //    int_print_mat(nbb[1], 1, idxb1, nbb[1]);
    //    d_print_mat(nbb[1], 1, lb1, nbb[1]);

    int *idxbN;
    int_zeros(&idxbN, nbb[N], 1);
    double *lbN;
    d_zeros(&lbN, nbb[N], 1);
    double *ubN;
    d_zeros(&ubN, nbb[N], 1);
    jj_end = nbx < nbb[N] ? nbx : nbb[N];
    for (jj = 0; jj < jj_end; jj++) {
        lbN[jj] = -4.0;  // xmin
        ubN[jj] = +4.0;  // xmax
        idxbN[jj] = jj;
    }
    for (; jj < nbb[N]; jj++) {
        lbN[jj] = -0.5;  // umin
        ubN[jj] = +0.5;  // umax
        idxbN[jj] = jj;
    }
    //    int_print_mat(nbb[N], 1, idxbN, nbb[N]);
    //    d_print_mat(nbb[N], 1, lbN, nbb[N]);

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
    for (ii = 0; ii < ngN; ii++) CN[ii * (ngN + 1)] = 1.0;
    //    d_print_mat(ngN, nx, CN, ngN);
    double *lgN;
    d_zeros(&lgN, ngN, 1);  // force all states to 0 at the last stage
    double *ugN;
    d_zeros(&ugN, ngN, 1);  // force all states to 0 at the last stage

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

#if defined(ELIMINATE_X0)
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
#endif

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

#if defined(ELIMINATE_X0)
    hA[0] = A0;
    hb[0] = b0;
    hQ[0] = Q0;
    hS[0] = S0;
    hq[0] = q0;
    hr[0] = r0;
#else
    hA[0] = A;
    hb[0] = b;
    hQ[0] = Q;
    hS[0] = S;
    hq[0] = q;
    hr[0] = r;
#endif
    hB[0] = B;
    hR[0] = R;
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

    for (ii = 0; ii < N; ii++) {
        d_zeros(&hx[ii], nxx[ii], 1);
        d_zeros(&hu[ii], nuu[ii], 1);
        d_zeros(&hpi[ii], nxx[ii + 1], 1);
        d_zeros(&hlam[ii], 2 * nbb[ii] + 2 * ngg[ii], 1);
        d_zeros(&ht[ii], 2 * nbb[ii] + 2 * ngg[ii], 1);
    }
    d_zeros(&hx[N], nxx[N], 1);
    d_zeros(&hlam[N], 2 * nbb[N] + 2 * ngg[N], 1);
    d_zeros(&ht[N], 2 * nbb[N] + 2 * ngg[N], 1);

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
    qp_out.t = ht;  // XXX why also the slack variables ???

    /************************************************
     * solver arguments (fully sparse)
     ************************************************/

    // solver arguments
    ocp_qp_condensing_hpipm_args *hpipm_args = ocp_qp_condensing_hpipm_create_arguments(&qp_in);
//    hpipm_args->mu_max = TOL;
//    hpipm_args->iter_max = MAXITER;
//    hpipm_args->alpha_min = MINSTEP;
    hpipm_args->mu0 = 1.0;  // 0.0

    /************************************************
     * work space (fully sparse)
     ************************************************/

    int work_space_size =
        ocp_qp_condensing_hpipm_calculate_workspace_size(&qp_in, hpipm_args);
    printf("\nwork space size: %d bytes\n", work_space_size);
    void *workspace = malloc(work_space_size);

    //    void *mem;
    //    ocp_qp_hpipm_create_memory(&qp_in, hpipm_args, &mem);
    int memory_size =
        ocp_qp_condensing_hpipm_calculate_memory_size(&qp_in, hpipm_args);
    printf("\nmemory: %d bytes\n", memory_size);
    void *memory = malloc(memory_size);

    ocp_qp_condensing_hpipm_memory *hpipm_memory =
        ocp_qp_condensing_hpipm_create_memory(&qp_in, hpipm_args);

    /************************************************
     * call the solver (fully sparse)
     ************************************************/

    int return_value;

    acados_timer timer;
    acados_tic(&timer);

    //  nrep = 1;
    for (rep = 0; rep < nrep; rep++) {
        // call the QP OCP solver
        //        return_value = ocp_qp_hpipm(&qp_in, &qp_out, hpipm_args,
        //        workspace);
        return_value =
            ocp_qp_condensing_hpipm(&qp_in, &qp_out, hpipm_args, hpipm_memory, workspace);
    }

    real_t time = acados_toc(&timer)/nrep;

    if (return_value == ACADOS_SUCCESS)
        printf("\nACADOS status: solution found in %d iterations\n",
               hpipm_memory->iter);

    if (return_value == ACADOS_MAXITER)
        printf("\nACADOS status: maximum number of iterations reached\n");

    if (return_value == ACADOS_MINSTEP)
        printf("\nACADOS status: below minimum step size length\n");

    printf("\nu = \n");
    for (ii = 0; ii < N; ii++) d_print_mat(1, nuu[ii], hu[ii], 1);

    printf("\nx = \n");
    for (ii = 0; ii <= N; ii++) d_print_mat(1, nxx[ii], hx[ii], 1);

    printf("\npi = \n");
    for (ii = 0; ii < N; ii++) d_print_mat(1, nxx[ii+1], hpi[ii], 1);

    printf("\nlam = \n");
    for (ii = 0; ii <= N; ii++) d_print_mat(1, 2*nbb[ii]+2*ngg[ii], hlam[ii], 1);

    printf("\n");
    printf(" inf norm res: %e, %e, %e, %e, %e\n", hpipm_memory->inf_norm_res[0],
           hpipm_memory->inf_norm_res[1], hpipm_memory->inf_norm_res[2],
           hpipm_memory->inf_norm_res[3], hpipm_memory->inf_norm_res[4]);
    printf("\n");
    printf(
        " Solution time for %d IPM iterations, averaged over %d runs: %5.2e "
        "seconds\n", hpipm_memory->iter, nrep, time);
    printf("\n\n");

    /************************************************
     * free memory
     ************************************************/

    d_free(A);
    d_free(B);
    d_free(b);
    d_free(x0);
    d_free(Q);
    d_free(S);
    d_free(R);
    d_free(q);
    d_free(r);
#if defined(ELIMINATE_X0)
    d_free(A0);
    d_free(b0);
    d_free(Q0);
    d_free(S0);
    d_free(q0);
    d_free(r0);
#endif
    int_free(idxb0);
    d_free(lb0);
    d_free(ub0);
    int_free(idxb1);
    d_free(lb1);
    d_free(ub1);
    int_free(idxbN);
    d_free(lbN);
    d_free(ubN);
    d_free(C);
    d_free(D);
    d_free(lg);
    d_free(ug);
    d_free(CN);
    d_free(lgN);
    d_free(ugN);

    for (ii = 0; ii < N; ii++) {
        d_free(hx[ii]);
        d_free(hu[ii]);
        d_free(hpi[ii]);
        d_free(hlam[ii]);
        d_free(ht[ii]);
    }
    d_free(hx[N]);
    d_free(hlam[N]);
    d_free(ht[N]);

    free(workspace);
    free(memory);

    return 0;
}
