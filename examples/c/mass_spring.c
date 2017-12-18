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
#include <stdio.h>
#include <stdlib.h>
// blasfeo
#include "blasfeo/include/blasfeo_target.h"
#include "blasfeo/include/blasfeo_common.h"
#include "blasfeo/include/blasfeo_d_aux_ext_dep.h"
#include "blasfeo/include/blasfeo_i_aux_ext_dep.h"
// hpipm
#include "hpipm/include/hpipm_d_ocp_qp.h"
// acados
#include "acados/utils/math.h"
// acados_c
#include <acados_c/ocp_qp.h>


/*****************************************************************************************
* Mass-spring system: nx/2 masses connected each other with springs (in a row),
* and the first and the last one to walls. nu (<=nx) controls act on the first nu
* masses. The system is sampled with sampling time Ts.
******************************************************************************************/

void mass_spring_system(double Ts, int nx, int nu, double *A, double *B, double *b) {

    int nx2 = nx * nx;

    int info = 0;

    int pp = nx / 2;  // number of masses

    /***********************************************
    * build the continuous time system
    ***********************************************/

    double *T;
    d_zeros(&T, pp, pp);

    for (int ii = 0; ii < pp; ii++) T[ii * (pp + 1)] = -2;
    for (int ii = 0; ii < pp - 1; ii++) T[ii * (pp + 1) + 1] = 1;
    for (int ii = 1; ii < pp; ii++) T[ii * (pp + 1) - 1] = 1;

    double *Z;
    d_zeros(&Z, pp, pp);
    double *I;
    d_zeros(&I, pp, pp);
    for (int ii = 0; ii < pp; ii++) I[ii * (pp + 1)] = 1.0;  // I = eye(pp);
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
    for (int ii = 0; ii < nu; ii++) I[ii * (nu + 1)] = 1.0;  // I = eye(nu);
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
    for (int ii = 0; ii < nx; ii++) I[ii * (nx + 1)] = 1.0;  // I = eye(nx);
    dmcopy(nx, nx, A, nx, T, nx);
    daxpy_3l(nx2, -1.0, I, T);
    dgemm_nn_3l(nx, nu, nx, T, nx, Bc, nx, B, nx);
    free(T);
    free(I);

    int *ipiv = (int *)malloc(nx * sizeof(int));
    dgesv_3l(nx, nu, Ac, nx, ipiv, B, nx, &info);
    free(ipiv);

    free(Ac);
    free(Bc);
    free(bb);
}



ocp_qp_in *create_ocp_qp_in_mass_spring( ) {

    int nx_ = 8;   // number of states (it has to be even for the mass-spring system test problem)

    int nu_ = 3;   // number of inputs (controllers) (it has to be at least 1 and
                   // at most nx_/2 for the mass-spring system test problem)

    int N = 15;    // horizon length
    int nb_ = 11;  // number of box constrained inputs and states
    int ng_ = 0;   // 4;  // number of general constraints

    int num_of_stages_equal_to_zero = 4;  // number of states to be enforced to zero at last stage

    #ifdef GENERAL_CONSTRAINT_AT_TERMINAL_STAGE
    int ngN = num_of_stages_equal_to_zero;  // number of general constraints at the last stage
    #else
    int ngN = 0;
    #endif

    int nbu = nu_ < nb_ ? nu_ : nb_;
    int nbx = nb_ - nu_ > 0 ? nb_ - nu_ : 0;

    int nx[N+1];
#if defined(ELIMINATE_X0)
    nx[0] = 0;
#else
    nx[0] = nx_;
#endif
    for (int ii = 1; ii <= N; ii++) {
        nx[ii] = nx_;
    }

    int nu[N+1];
    for (int ii = 0; ii < N; ii++) {
        nu[ii] = nu_;
    }
    nu[N] = 0;

    int nb[N+1];
#if defined(ELIMINATE_X0)
    nb[0] = nbu;
#else
    nb[0] = nb_;
#endif
    for (int ii = 1; ii < N; ii++) {
        nb[ii] = nb_;
    }
    nb[N] = nbx;

    int ng[N+1];
    for (int ii = 0; ii < N; ii++) {
        ng[ii] = ng_;
    }
    ng[N] = ngN;

	int ns[N+1];
	for (int ii = 0; ii <= N; ii++) {
        ns[ii] = 0;
    }

    printf("Test problem: mass-spring system with %d masses and %d controls.\n\n", nx_ / 2, nu_);
    printf("MPC problem size: %d states, %d inputs, %d horizon length, %d two-sided box "
           "constraints, %d two-sided general constraints.\n\n", nx_, nu_, N, nb_, ng_);

    /************************************************
    * dynamical system
    ************************************************/

    // state space matrices & initial state
    double *A;
    d_zeros(&A, nx_, nx_);  // states update matrix
    double *B;
    d_zeros(&B, nx_, nu_);  // inputs matrix
    double *b;
    d_zeros(&b, nx_, 1);  // states offset
    double *x0;
    d_zeros(&x0, nx_, 1);  // initial state

    // mass-spring system
    double Ts = 0.5;
    mass_spring_system(Ts, nx_, nu_, A, B, b);

    // TODO(dimitris): @giaf, why do we overwrite b here?
    for (int jj = 0; jj < nx_; jj++) b[jj] = 0.1;

    // initial state
    for (int jj = 0; jj < nx_; jj++) x0[jj] = 0;
    x0[0] = 2.5;
    x0[1] = 2.5;

//    d_print_mat(nx_, nx_, A, nx_);
//    d_print_mat(nx_, nu_, B, nx_);
//    d_print_mat(nx_, 1, b, nx_);
//    d_print_mat(nx_, 1, x0, nx_);

#if defined(ELIMINATE_X0)
    // compute b0 = b + A*x0
    double *b0;
    d_zeros(&b0, nx_, 1);
    dcopy_3l(nx_, b, 1, b0, 1);
    dgemv_n_3l(nx_, nx_, A, nx_, x0, b0);
    //    d_print_mat(nx_, 1, b, nx_);
    //    d_print_mat(nx_, 1, b0, nx_);

    // then A0 is a matrix of size 0x0
    double *A0;
    d_zeros(&A0, 0, 0);
#endif

    /************************************************
    * box constraints
    ************************************************/

    int jj_end;

    int *idxb0;
    int_zeros(&idxb0, nb[0], 1);
    double *lb0;
    d_zeros(&lb0, nb[0], 1);
    double *ub0;
    d_zeros(&ub0, nb[0], 1);
#if defined(ELIMINATE_X0)
    for (int jj = 0; jj < nb[0]; jj++) {
        lb0[jj] = -0.5;  // umin
        ub0[jj] = +0.5;  // umin
        idxb0[jj] = jj;
    }
#else
    jj_end = nu[0] < nb[0] ? nu[0] : nb[0];
    for (int jj = 0; jj < jj_end; jj++) {
        lb0[jj] = -0.5;  // umin
        ub0[jj] =  0.5;  // umax
        idxb0[jj] = jj;
    }
    for (int jj = jj_end; jj < nb[0]; jj++) {
        lb0[jj] = x0[jj-jj_end];  // initial state
        ub0[jj] = x0[jj-jj_end];  // initial state
        idxb0[jj] = jj;
    }
#endif

    int *idxb1;
    int_zeros(&idxb1, nb[1], 1);
    double *lb1;
    d_zeros(&lb1, nb[1], 1);
    double *ub1;
    d_zeros(&ub1, nb[1], 1);
    jj_end = nu[1] < nb[1] ? nu[1] : nb[1];
    for (int jj = 0; jj < jj_end; jj++) {
        lb1[jj] = -0.5;  // umin
        ub1[jj] = +0.5;  // umax
        idxb1[jj] = jj;
    }
    for (int jj = jj_end; jj < nb[1]; jj++) {
        lb1[jj] = -4.0;  // xmin
        ub1[jj] = +4.0;  // xmax
        idxb1[jj] = jj;
    }
    //    int_print_mat(nb[1], 1, idxb1, nb[1]);
    //    d_print_mat(nb[1], 1, lb1, nb[1]);

    int *idxbN;
    int_zeros(&idxbN, nb[N], 1);
    double *lbN;
    d_zeros(&lbN, nb[N], 1);
    double *ubN;
    d_zeros(&ubN, nb[N], 1);
    jj_end = nu[N] < nb[N] ? nu[N] : nb[N];
    for (int jj = 0; jj < jj_end; jj++) {
        lbN[jj] = -0.5;  // umin
        ubN[jj] = +0.5;  // umax
        idxbN[jj] = jj;
    }
    for (int jj = jj_end; jj < nb[N]; jj++)
    {
        lbN[jj] = -4.0;  // xmin
        ubN[jj] = +4.0;  // xmax
        idxbN[jj] = jj;
    }

    #ifndef GENERAL_CONSTRAINT_AT_TERMINAL_STAGE
    for (int jj = nu[N]; jj < num_of_stages_equal_to_zero; jj++)
    {
        lbN[jj] = 0.0;
        ubN[jj] = 0.0;
        idxbN[jj] = jj;
    }
    #endif

    //    int_print_mat(nb[N], 1, idxbN, nb[N]);
    //    d_print_mat(nb[N], 1, lbN, nb[N]);

    /************************************************
    * general constraints
    ************************************************/

    double *C;
    d_zeros(&C, ng_, nx_);
    double *D;
    d_zeros(&D, ng_, nu_);
    double *lg;
    d_zeros(&lg, ng_, 1);
    double *ug;
    d_zeros(&ug, ng_, 1);

    double *CN;
    d_zeros(&CN, ngN, nx_);
    for (int ii = 0; ii < ngN; ii++) CN[ii * (ngN + 1)] = 1.0;
    //    d_print_mat(ngN, nx_, CN, ngN);
    double *lgN;
    d_zeros(&lgN, ngN, 1);  // force all states to 0 at the last stage
    double *ugN;
    d_zeros(&ugN, ngN, 1);  // force all states to 0 at the last stage

    /************************************************
    * cost function
    ************************************************/

    double *Q;
    d_zeros(&Q, nx_, nx_);
    for (int ii = 0; ii < nx_; ii++) Q[ii * (nx_ + 1)] = 1.0;

    double *R;
    d_zeros(&R, nu_, nu_);
    for (int ii = 0; ii < nu_; ii++) R[ii * (nu_ + 1)] = 2.0;

    double *S;
    d_zeros(&S, nu_, nx_);

    double *q;
    d_zeros(&q, nx_, 1);
    for (int ii = 0; ii < nx_; ii++) q[ii] = 0.1;

    double *r;
    d_zeros(&r, nu_, 1);
    for (int ii = 0; ii < nu_; ii++) r[ii] = 0.2;

#if defined(ELIMINATE_X0)
    // Q0 and q0 are matrices of size 0
    double *Q0;
    d_zeros(&Q0, 0, 0);
    double *q0;
    d_zeros(&q0, 0, 1);

    // compute r0 = r + S*x0
    double *r0;
    d_zeros(&r0, nu_, 1);
    dcopy_3l(nu_, r, 1, r0, 1);
    dgemv_n_3l(nu_, nx_, S, nu_, x0, r0);

    // then S0 is a matrix of size nux0
    double *S0;
    d_zeros(&S0, nu_, 0);
#endif

    /************************************************
    * problem data
    ************************************************/

    double *hA[N];
    double *hB[N];
    double *hb[N];
    double *hQ[N+1];
    double *hS[N+1];
    double *hR[N+1];
    double *hq[N+1];
    double *hr[N+1];
    double *hlb[N+1];
    double *hub[N+1];
    int *hidxb[N+1];
    double *hC[N+1];
    double *hD[N+1];
    double *hlg[N+1];
    double *hug[N+1];

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
    for (int ii = 1; ii < N; ii++) {
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

    int nbx_vec[N+1];
    int nbu_vec[N+1];
    for (int ii = 0; ii < N+1; ii++)
    {
        nbx_vec[ii] = 0;
        nbu_vec[ii] = 0;
        for (int jj = 0; jj < nb[ii]; jj++)
            {
            if (hidxb[ii][jj] < nu[ii])
                {
                nbu_vec[ii]++;
                }
            else
                {
                nbx_vec[ii]++;
                }
            }
    }

    ocp_qp_dims dims;

    dims.N = N;
    dims.nx = nx;
    dims.nu = nu;
    dims.nb = nb;
    dims.ng = ng;
    dims.ns = ns;
    dims.nbu = nbu_vec;
    dims.nbx = nbx_vec;

    ocp_qp_in *qp_in = create_ocp_qp_in(&dims);
	d_cvt_colmaj_to_ocp_qp(hA, hB, hb, hQ, hS, hR, hq, hr, hidxb, hlb, hub, hC, hD, hlg, hug, NULL, NULL, NULL, NULL, NULL, qp_in);

    // free objective
    free(Q);
    free(S);
    free(R);
    free(q);
    free(r);

#if defined(ELIMINATE_X0)
    free(Q0);
    free(q0);
    free(r0);
    free(S0);
#endif

    // free dynamics
    free(A);
    free(B);
    free(b);
    free(x0);

#if defined(ELIMINATE_X0)
    free(b0);
    free(A0);
#endif

    // free constraints
    free(C);
    free(D);
    free(lg);
    free(ug);
    free(CN);
    free(lgN);
    free(ugN);

    free(idxb0);
    free(idxb1);
    free(idxbN);
    free(lb0);
    free(lb1);
    free(lbN);
    free(ub0);
    free(ub1);
    free(ubN);

    return qp_in;
}