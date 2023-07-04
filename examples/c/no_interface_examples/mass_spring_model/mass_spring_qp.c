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
#include <stdio.h>
#include <stdlib.h>

// blasfeo
#include "blasfeo/include/blasfeo_d_aux_ext_dep.h"
#include "blasfeo/include/blasfeo_i_aux_ext_dep.h"

// hpipm
#include "hpipm/include/hpipm_d_ocp_qp.h"

// acados
#include "acados/utils/math.h"
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/ocp_qp/ocp_qp_xcond_solver.h"

// acados_c
#include <acados_c/ocp_qp_interface.h>



/*****************************************************************************************
* Mass-spring system: nx/2 masses connected each other with springs (in a row),
* and the first and the last one to walls. nu (<=nx) controls act on the first nu
* masses. The system is sampled with sampling time Ts.
******************************************************************************************/

static void mass_spring_system(double Ts, int nx, int nu, double *A, double *B, double *b) {

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



ocp_qp_xcond_solver_dims *create_ocp_qp_dims_mass_spring(ocp_qp_xcond_solver_config *config, int N, int nx_, int nu_, int nb_, int ng_, int ngN)
{

    int nbu_ = nu_<nb_ ? nu_ : nb_;
    int nbx_ = nb_ - nu_ > 0 ? nb_ - nu_ : 0;

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

    int nbu[N+1];
    for (int ii = 0; ii < N; ii++)
    {
        nbu[ii] = nbu_;
    }
    nbu[N] = 0;

    int nbx[N+1];
#if defined(ELIMINATE_X0)
    nbx[0] = 0;
#else
    nbx[0] = nx_;
#endif
    for (int ii = 1; ii <= N; ii++)
    {
        nbx[ii] = nbx_;
    }

    // int nb[N+1];
    // for (int ii = 0; ii <= N; ii++) {
    //     nb[ii] = nbu[ii]+nbx[ii];
    // }

    int ng[N+1];
    for (int ii = 0; ii < N; ii++) {
        ng[ii] = ng_;
    }
    ng[N] = ngN;

    int ns[N+1];
    for (int ii = 0; ii <= N; ii++) {
        ns[ii] = 0;
    }

    int nbxe[N+1];
#if defined(ELIMINATE_X0)
    nbxe[0] = 0;
#else
    nbxe[0] = nx_;
#endif
    for (int ii = 1; ii <= N; ii++)
    {
        nbxe[ii] = 0;
    }


    // dims
    int dims_size = ocp_qp_xcond_solver_dims_calculate_size(config, N);
    void *dims_mem = malloc(dims_size);
    ocp_qp_xcond_solver_dims *dims = ocp_qp_xcond_solver_dims_assign(config, N, dims_mem);

    for (int ii=0; ii<=N; ii++)
    {
		config->dims_set(config, dims, ii, "nx", &nx[ii]);
		config->dims_set(config, dims, ii, "nu", &nu[ii]);
		config->dims_set(config, dims, ii, "nbx", &nbx[ii]);
		config->dims_set(config, dims, ii, "nbu", &nbu[ii]);
		config->dims_set(config, dims, ii, "ng", &ng[ii]);
		config->dims_set(config, dims, ii, "ns", &ns[ii]);
    }
	config->dims_set(config, dims, 0, "nbxe", &nbxe[0]);

//d_ocp_qp_dim_print(dims->orig_dims);

    return dims;

}



ocp_qp_in *create_ocp_qp_in_mass_spring(ocp_qp_dims *dims)
{

    /************************************************
    * extract dims
    ************************************************/

    int N = dims->N;
    int *nx = dims->nx;
    int *nu = dims->nu;
    int *nb = dims->nb;
    int *ng = dims->ng;
    // int *ns = dims->ns;

    int nx_ = nx[1];
    int nu_ = nu[1];
    int ng_ = ng[1];
    int ngN = ng[N];

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
    int *idxbxe0;
    int_zeros(&idxbxe0, nx[0], 1);
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
	for(int jj=0; jj<nx[0]; jj++)
		idxbxe0[jj] = jj;

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
    for (int jj = nu[N]; jj < ngN; jj++)
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


    ocp_qp_in *qp_in = ocp_qp_in_create(dims);

//    d_cvt_colmaj_to_ocp_qp(hA, hB, hb, hQ, hS, hR, hq, hr, hidxb, hlb, hub, hC, hD, hlg, hug, NULL, NULL, NULL, NULL, NULL, NULL, NULL, qp_in);
	int ii;
	for(ii=0; ii<N; ii++)
	{
		d_ocp_qp_set_A(ii, hA[ii], qp_in);
		d_ocp_qp_set_B(ii, hB[ii], qp_in);
		d_ocp_qp_set_b(ii, hb[ii], qp_in);
		d_ocp_qp_set_R(ii, hR[ii], qp_in);
		d_ocp_qp_set_S(ii, hS[ii], qp_in);
		d_ocp_qp_set_Q(ii, hQ[ii], qp_in);
		d_ocp_qp_set_r(ii, hr[ii], qp_in);
		d_ocp_qp_set_q(ii, hq[ii], qp_in);
		d_ocp_qp_set_idxb(ii, hidxb[ii], qp_in);
		d_ocp_qp_set_lb(ii, hlb[ii], qp_in);
		d_ocp_qp_set_ub(ii, hub[ii], qp_in);
		d_ocp_qp_set_C(ii, hC[ii], qp_in);
		d_ocp_qp_set_D(ii, hD[ii], qp_in);
		d_ocp_qp_set_lg(ii, hlg[ii], qp_in);
		d_ocp_qp_set_ug(ii, hug[ii], qp_in);
	}
	d_ocp_qp_set_R(ii, hR[ii], qp_in);
	d_ocp_qp_set_S(ii, hS[ii], qp_in);
	d_ocp_qp_set_Q(ii, hQ[ii], qp_in);
	d_ocp_qp_set_r(ii, hr[ii], qp_in);
	d_ocp_qp_set_q(ii, hq[ii], qp_in);
	d_ocp_qp_set_idxb(ii, hidxb[ii], qp_in);
	d_ocp_qp_set_lb(ii, hlb[ii], qp_in);
	d_ocp_qp_set_ub(ii, hub[ii], qp_in);
	d_ocp_qp_set_C(ii, hC[ii], qp_in);
	d_ocp_qp_set_D(ii, hD[ii], qp_in);
	d_ocp_qp_set_lg(ii, hlg[ii], qp_in);
	d_ocp_qp_set_ug(ii, hug[ii], qp_in);
	//
	d_ocp_qp_set_idxbxe(0, idxbxe0, qp_in);

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
    free(idxbxe0);

    return qp_in;
}



ocp_qp_dims *create_ocp_qp_dims_mass_spring_soft_constr(int N, int nx_, int nu_, int nb_, int ng_, int ngN)
{

    int nbu_ = nu_<nb_ ? nu_ : nb_;
    int nbx_ = nb_ - nu_ > 0 ? nb_ - nu_ : 0;

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

    int nbu[N+1];
    for (int ii = 0; ii < N; ii++)
    {
        nbu[ii] = nbu_;
    }
    nbu[N] = 0;

    int nbx[N+1];
#if defined(ELIMINATE_X0)
    nbx[0] = 0;
#else
    nbx[0] = nx_;
#endif
    for (int ii = 1; ii <= N; ii++)
    {
        nbx[ii] = nbx_;
    }

    int nb[N+1];
    for (int ii = 0; ii <= N; ii++) {
        nb[ii] = nbu[ii]+nbx[ii];
    }

    int ng[N+1];
    for (int ii = 0; ii < N; ii++) {
        ng[ii] = ng_;
    }
    ng[N] = ngN;

    int ns[N+1];
    for (int ii = 0; ii <= N; ii++) {
        ns[ii] = nbx[ii]+ng[ii];
    }



    // dims
    int dims_size = ocp_qp_dims_calculate_size(N);
    void *dims_mem = malloc(dims_size);
    ocp_qp_dims *dims = ocp_qp_dims_assign(N, dims_mem);

    dims->N = N;
    for (int ii=0; ii<=N; ii++)
    {
        dims->nx[ii] = nx[ii];
        dims->nu[ii] = nu[ii];
        dims->nb[ii] = nb[ii];
        dims->ng[ii] = ng[ii];
        dims->ns[ii] = ns[ii];
        dims->nsbx[ii] = nbx[ii];
        dims->nsbu[ii] = 0;
        dims->nsg[ii] = ng[ii];
        dims->nbu[ii] = nbu[ii];
        dims->nbx[ii] = nbx[ii];
    }

    return dims;

}



ocp_qp_in *create_ocp_qp_in_mass_spring_soft_constr(ocp_qp_dims *dims)
{

    int ii;

    /************************************************
    * extract dims
    ************************************************/

    int N = dims->N;
    int *nx = dims->nx;
    int *nu = dims->nu;
    int *nb = dims->nb;
    int *ng = dims->ng;
    int *ns = dims->ns;

    int nx_ = nx[1];
    int nu_ = nu[1];
    int ng_ = ng[1];
    int ngN = ng[N];

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
    for (int jj = 0; jj < nx_; jj++) b[jj] = 0.0;

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
        lb1[jj] = -1.0;  // xmin
        ub1[jj] = +1.0;  // xmax
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
        lbN[jj] = -1.0;  // xmin
        ubN[jj] = +1.0;  // xmax
        idxbN[jj] = jj;
    }

    #ifndef GENERAL_CONSTRAINT_AT_TERMINAL_STAGE
    for (int jj = nu[N]; jj < ngN; jj++)
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
    * soft constraints
    ************************************************/

    double *Zl0; d_zeros(&Zl0, ns[0], 1);
    for(ii=0; ii<ns[0]; ii++)
        Zl0[ii] = 1e3;
    double *Zu0; d_zeros(&Zu0, ns[0], 1);
    for(ii=0; ii<ns[0]; ii++)
        Zu0[ii] = 1e3;
    double *zl0; d_zeros(&zl0, ns[0], 1);
    for(ii=0; ii<ns[0]; ii++)
        zl0[ii] = 1e2;
    double *zu0; d_zeros(&zu0, ns[0], 1);
    for(ii=0; ii<ns[0]; ii++)
        zu0[ii] = 1e2;
    int *idxs0; int_zeros(&idxs0, ns[0], 1);
    for(ii=0; ii<ns[0]; ii++)
        idxs0[ii] = nu[0]+ii;
    double *d_ls0; d_zeros(&d_ls0, ns[0], 1);
    for(ii=0; ii<ns[0]; ii++)
        d_ls0[ii] = 0.0;
    double *d_us0; d_zeros(&d_us0, ns[0], 1);
    for(ii=0; ii<ns[0]; ii++)
        d_us0[ii] = 0.0;

    double *Zl1; d_zeros(&Zl1, ns[1], 1);
    for(ii=0; ii<ns[1]; ii++)
        Zl1[ii] = 1e3;
    double *Zu1; d_zeros(&Zu1, ns[1], 1);
    for(ii=0; ii<ns[1]; ii++)
        Zu1[ii] = 1e3;
    double *zl1; d_zeros(&zl1, ns[1], 1);
    for(ii=0; ii<ns[1]; ii++)
        zl1[ii] = 1e2;
    double *zu1; d_zeros(&zu1, ns[1], 1);
    for(ii=0; ii<ns[1]; ii++)
        zu1[ii] = 1e2;
    int *idxs1; int_zeros(&idxs1, ns[1], 1);
    for(ii=0; ii<ns[1]; ii++)
        idxs1[ii] = nu[1]+ii;
    double *d_ls1; d_zeros(&d_ls1, ns[1], 1);
    for(ii=0; ii<ns[1]; ii++)
        d_ls1[ii] = 0.0;
    double *d_us1; d_zeros(&d_us1, ns[1], 1);
    for(ii=0; ii<ns[1]; ii++)
        d_us1[ii] = 0.0;

    double *ZlN; d_zeros(&ZlN, ns[N], 1);
    for(ii=0; ii<ns[N]; ii++)
        ZlN[ii] = 1e3;
    double *ZuN; d_zeros(&ZuN, ns[N], 1);
    for(ii=0; ii<ns[N]; ii++)
        ZuN[ii] = 1e3;
    double *zlN; d_zeros(&zlN, ns[N], 1);
    for(ii=0; ii<ns[N]; ii++)
        zlN[ii] = 1e2;
    double *zuN; d_zeros(&zuN, ns[N], 1);
    for(ii=0; ii<ns[N]; ii++)
        zuN[ii] = 1e2;
    int *idxsN; int_zeros(&idxsN, ns[N], 1);
    for(ii=0; ii<ns[N]; ii++)
        idxsN[ii] = nu[N]+ii;
    double *d_lsN; d_zeros(&d_lsN, ns[N], 1);
    for(ii=0; ii<ns[N]; ii++)
        d_lsN[ii] = 0.0;
    double *d_usN; d_zeros(&d_usN, ns[N], 1);
    for(ii=0; ii<ns[N]; ii++)
        d_usN[ii] = 0.0;

    /************************************************
    * cost function
    ************************************************/

    double *Q;
    d_zeros(&Q, nx_, nx_);
    for (int ii = 0; ii < nx_; ii++) Q[ii * (nx_ + 1)] = 0.0;

    double *R;
    d_zeros(&R, nu_, nu_);
    for (int ii = 0; ii < nu_; ii++) R[ii * (nu_ + 1)] = 2.0;

    double *S;
    d_zeros(&S, nu_, nx_);

    double *q;
    d_zeros(&q, nx_, 1);
    for (int ii = 0; ii < nx_; ii++) q[ii] = 0.0;

    double *r;
    d_zeros(&r, nu_, 1);
    for (int ii = 0; ii < nu_; ii++) r[ii] = 0.0;

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
    double *hZl[N+1];
    double *hZu[N+1];
    double *hzl[N+1];
    double *hzu[N+1];
    int *hidxs[N+1];

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
    hZl[0] = Zl0;
    hZu[0] = Zu0;
    hzl[0] = zl0;
    hzu[0] = zu0;
    hidxs[0] = idxs0;
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
        hZl[ii] = Zl1;
        hZu[ii] = Zu1;
        hzl[ii] = zl1;
        hzu[ii] = zu1;
        hidxs[ii] = idxs1;
    }
    hQ[N] = Q;  // or maybe initialize to the solution of the DARE???
    hq[N] = q;  // or maybe initialize to the solution of the DARE???
    hlb[N] = lbN;
    hub[N] = ubN;
    hidxb[N] = idxbN;
    hC[N] = CN;
    hlg[N] = lgN;
    hug[N] = ugN;
    hZl[N] = ZlN;
    hZu[N] = ZuN;
    hzl[N] = zlN;
    hzu[N] = zuN;
    hidxs[N] = idxsN;


    ocp_qp_in *qp_in = ocp_qp_in_create(dims);

//    d_cvt_colmaj_to_ocp_qp(hA, hB, hb, hQ, hS, hR, hq, hr, hidxb, hlb, hub, hC, hD, hlg, hug, hZl, hZu, hzl, hzu, hidxs, hd_ls, hd_us, qp_in);
	for(ii=0; ii<N; ii++)
	{
		d_ocp_qp_set_A(ii, hA[ii], qp_in);
		d_ocp_qp_set_B(ii, hB[ii], qp_in);
		d_ocp_qp_set_b(ii, hb[ii], qp_in);
		d_ocp_qp_set_R(ii, hR[ii], qp_in);
		d_ocp_qp_set_S(ii, hS[ii], qp_in);
		d_ocp_qp_set_Q(ii, hQ[ii], qp_in);
		d_ocp_qp_set_r(ii, hr[ii], qp_in);
		d_ocp_qp_set_q(ii, hq[ii], qp_in);
		d_ocp_qp_set_idxb(ii, hidxb[ii], qp_in);
		d_ocp_qp_set_lb(ii, hlb[ii], qp_in);
		d_ocp_qp_set_ub(ii, hub[ii], qp_in);
		d_ocp_qp_set_C(ii, hC[ii], qp_in);
		d_ocp_qp_set_D(ii, hD[ii], qp_in);
		d_ocp_qp_set_lg(ii, hlg[ii], qp_in);
		d_ocp_qp_set_ug(ii, hug[ii], qp_in);
		d_ocp_qp_set_Zl(ii, hZl[ii], qp_in);
		d_ocp_qp_set_Zu(ii, hZu[ii], qp_in);
		d_ocp_qp_set_zl(ii, hzl[ii], qp_in);
		d_ocp_qp_set_zu(ii, hzu[ii], qp_in);
		d_ocp_qp_set_idxb(ii, hidxs[ii], qp_in);
	}
	d_ocp_qp_set_R(ii, hR[ii], qp_in);
	d_ocp_qp_set_S(ii, hS[ii], qp_in);
	d_ocp_qp_set_Q(ii, hQ[ii], qp_in);
	d_ocp_qp_set_r(ii, hr[ii], qp_in);
	d_ocp_qp_set_q(ii, hq[ii], qp_in);
	d_ocp_qp_set_idxb(ii, hidxb[ii], qp_in);
	d_ocp_qp_set_lb(ii, hlb[ii], qp_in);
	d_ocp_qp_set_ub(ii, hub[ii], qp_in);
	d_ocp_qp_set_C(ii, hC[ii], qp_in);
	d_ocp_qp_set_D(ii, hD[ii], qp_in);
	d_ocp_qp_set_lg(ii, hlg[ii], qp_in);
	d_ocp_qp_set_ug(ii, hug[ii], qp_in);
	d_ocp_qp_set_Zl(ii, hZl[ii], qp_in);
	d_ocp_qp_set_Zu(ii, hZu[ii], qp_in);
	d_ocp_qp_set_zl(ii, hzl[ii], qp_in);
	d_ocp_qp_set_zu(ii, hzu[ii], qp_in);
	d_ocp_qp_set_idxb(ii, hidxs[ii], qp_in);

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

    d_free(Zl0);
    d_free(Zu0);
    d_free(zl0);
    d_free(zu0);
    int_free(idxs0);
    d_free(d_ls0);
    d_free(d_us0);
    d_free(Zl1);
    d_free(Zu1);
    d_free(zl1);
    d_free(zu1);
    int_free(idxs1);
    d_free(d_ls1);
    d_free(d_us1);
    d_free(ZlN);
    d_free(ZuN);
    d_free(zlN);
    d_free(zuN);
    int_free(idxsN);
    d_free(d_lsN);
    d_free(d_usN);

    return qp_in;
}
