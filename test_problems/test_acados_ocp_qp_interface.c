/**************************************************************************************************
*                                                                                                 *
* This file is part of HPMPC.                                                                     *
*                                                                                                 *
* HPMPC -- Library for High-Performance implementation of solvers for MPC.                        *
* Copyright (C) 2014-2015 by Technical University of Denmark. All rights reserved.                *
*                                                                                                 *
* HPMPC is free software; you can redistribute it and/or                                          *
* modify it under the terms of the GNU Lesser General Public                                      *
* License as published by the Free Software Foundation; either                                    *
* version 2.1 of the License, or (at your option) any later version.                              *
*                                                                                                 *
* HPMPC is distributed in the hope that it will be useful,                                        *
* but WITHOUT ANY WARRANTY; without even the implied warranty of                                  *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                                            *
* See the GNU Lesser General Public License for more details.                                     *
*                                                                                                 *
* You should have received a copy of the GNU Lesser General Public                                *
* License along with HPMPC; if not, write to the Free Software                                    *
* Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA                  *
*                                                                                                 *
* Author: Gianluca Frison, giaf (at) dtu.dk                                                       *
*                                                                                                 *
**************************************************************************************************/

// define problems size
#define NX 8
#define NU 3
#define NN 10
#define NB 11
#define NG 0
#define NGN 8

// define IP solver arguments && number of repetitions
#define NREP 100
#define ITMAX 10
#define TOL 1e-6
#define STEPMIN 1e-8

// system headers
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>

// HPMPC headers
#include "aux_d.h"
#include "mpc_solvers.h"
#include "target.h"
#include "block_size.h"
#include "tools.h"
#include "c_interface.h"

// flush denormals to zero
#if defined(TARGET_X64_AVX2) || defined(TARGET_X64_AVX) || defined(TARGET_X64_SSE3) || defined(TARGET_X86_ATOM) || defined(TARGET_AMD_SSE3)
#include <xmmintrin.h> // needed to flush to zero sub-normals with _MM_SET_FLUSH_ZERO_MODE (_MM_FLUSH_ZERO_ON); in the main()
#endif

// structure of solver arguments
struct ocp_qp_solver_args
	{
	double tol;
	int max_iter;
	double min_step;
	double mu0;
	double sigma_min;
	};

// header of the TODO function returning the work space size
	int ocp_qp_hpmpc_workspace_double(int N, int *nxx, int *nuu, int *nbb, int *ngg, struct ocp_qp_solver_args *args);



/************************************************ 
Mass-spring system: nx/2 masses connected each other with springs (in a row), and the first and the last one to walls. nu (<=nx) controls act on the first nu masses. The system is sampled with sampling time Ts. 
************************************************/
void mass_spring_system(double Ts, int nx, int nu, int N, double *A, double *B, double *b, double *x0)
	{

	int nx2 = nx*nx;

	int info = 0;

	int pp = nx/2; // number of masses
	
/************************************************
* build the continuous time system 
************************************************/
	
	double *T; d_zeros(&T, pp, pp);
	int ii;
	for(ii=0; ii<pp; ii++) T[ii*(pp+1)] = -2;
	for(ii=0; ii<pp-1; ii++) T[ii*(pp+1)+1] = 1;
	for(ii=1; ii<pp; ii++) T[ii*(pp+1)-1] = 1;

	double *Z; d_zeros(&Z, pp, pp);
	double *I; d_zeros(&I, pp, pp); for(ii=0; ii<pp; ii++) I[ii*(pp+1)]=1.0; // = eye(pp);
	double *Ac; d_zeros(&Ac, nx, nx);
	dmcopy(pp, pp, Z, pp, Ac, nx);
	dmcopy(pp, pp, T, pp, Ac+pp, nx);
	dmcopy(pp, pp, I, pp, Ac+pp*nx, nx);
	dmcopy(pp, pp, Z, pp, Ac+pp*(nx+1), nx); 
	free(T);
	free(Z);
	free(I);
	
	d_zeros(&I, nu, nu); for(ii=0; ii<nu; ii++) I[ii*(nu+1)]=1.0; //I = eye(nu);
	double *Bc; d_zeros(&Bc, nx, nu);
	dmcopy(nu, nu, I, nu, Bc+pp, nx);
	free(I);
	
/************************************************
* compute the discrete time system 
************************************************/

	double *bb; d_zeros(&bb, nx, 1);
	dmcopy(nx, 1, bb, nx, b, nx);
		
	dmcopy(nx, nx, Ac, nx, A, nx);
	dscal_3l(nx2, Ts, A);
	expm(nx, A);
	
	d_zeros(&T, nx, nx);
	d_zeros(&I, nx, nx); for(ii=0; ii<nx; ii++) I[ii*(nx+1)]=1.0; //I = eye(nx);
	dmcopy(nx, nx, A, nx, T, nx);
	daxpy_3l(nx2, -1.0, I, T);
	dgemm_nn_3l(nx, nu, nx, T, nx, Bc, nx, B, nx);
	
	int *ipiv = (int *) malloc(nx*sizeof(int));
	dgesv_3l(nx, nu, Ac, nx, ipiv, B, nx, &info);
	free(ipiv);

	free(Ac);
	free(Bc);
	free(bb);
	
			
/************************************************
* initial state 
************************************************/
	
	if(nx==4)
		{
		x0[0] = 5;
		x0[1] = 10;
		x0[2] = 15;
		x0[3] = 20;
		}
	else
		{
		int jj;
		for(jj=0; jj<nx; jj++)
			x0[jj] = 1;
		}

	}



int main()
	{
	
	printf("\n");
	printf("\n");
	printf("\n");
	printf(" HPMPC -- Library for High-Performance implementation of solvers for MPC.\n");
	printf(" Copyright (C) 2014-2015 by Technical University of Denmark. All rights reserved.\n");
	printf("\n");
	printf(" HPMPC is distributed in the hope that it will be useful,\n");
	printf(" but WITHOUT ANY WARRANTY; without even the implied warranty of\n");
	printf(" MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n");
	printf(" See the GNU Lesser General Public License for more details.\n");
	printf("\n");
	printf("\n");
	printf("\n");
	
#if defined(TARGET_X64_AVX2) || defined(TARGET_X64_AVX) || defined(TARGET_X64_SSE3) || defined(TARGET_X86_ATOM) || defined(TARGET_AMD_SSE3)
	_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON); // flush to zero subnormals !!! works only with one thread !!!
#endif

	int ii, jj;
	
	int rep, nrep=NREP;

	int nx = NX; // number of states (it has to be even for the mass-spring system test problem)
	int nu = NU; // number of inputs (controllers) (it has to be at least 1 and at most nx/2 for the mass-spring system test problem)
	int N  = NN; // horizon lenght
	int nb  = NB; // number of box constrained inputs and states
	int ng  = NG; //4;  // number of general constraints
	int ngN  = NGN; //4;  // number of general constraints at the last stage
	
	int nbu = nu<nb ? nu : nb ;
	int nbx = nb-nu>0 ? nb-nu : 0;


	// stage-wise variant size
	int nxx[N+1];
	nxx[0] = 0;
	for(ii=1; ii<=N; ii++)
		nxx[ii] = nx;

	int nuu[N+1];
	for(ii=0; ii<N; ii++)
		nuu[ii] = nu;
	nuu[N] = 0;

	int nbb[N+1];
	nbb[0] = nbu;
	for(ii=1; ii<N; ii++)
		nbb[ii] = nb;
	nbb[N] = nbx;

	int ngg[N+1];
	for(ii=0; ii<N; ii++)
		ngg[ii] = ng;
	ngg[N] = ngN;
	



	printf(" Test problem: mass-spring system with %d masses and %d controls.\n", nx/2, nu);
	printf("\n");
	printf(" MPC problem size: %d states, %d inputs, %d horizon length, %d two-sided box constraints, %d two-sided general constraints.\n", nx, nu, N, nb, ng);
	printf("\n");
	printf(" IP method parameters: predictor-corrector IP, double precision, %d maximum iterations, %5.1e exit tolerance in duality measure.\n", ITMAX, TOL);
	printf("\n");
#if defined(TARGET_X64_AVX2)
	printf(" HPMPC built for the AVX2 architecture\n");
#endif
#if defined(TARGET_X64_AVX)
	printf(" HPMPC built for the AVX architecture\n");
#endif
	printf("\n");

	int info = 0;
		

/************************************************
* dynamical system
************************************************/	

	// state space matrices & initial state
	double *A; d_zeros(&A, nx, nx); // states update matrix
	double *B; d_zeros(&B, nx, nu); // inputs matrix
	double *b; d_zeros(&b, nx, 1); // states offset
	double *x0; d_zeros_align(&x0, nx, 1); // initial state

	// mass-spring system
	double Ts = 0.5; // sampling time
	mass_spring_system(Ts, nx, nu, N, A, B, b, x0);
	
	for(jj=0; jj<nx; jj++)
		b[jj] = 0.1;
	
	for(jj=0; jj<nx; jj++)
		x0[jj] = 0;
	x0[0] = 2.5;
	x0[1] = 2.5;

//	d_print_mat(nx, nx, A, nx);
//	d_print_mat(nx, nu, B, nx);
//	d_print_mat(nx, 1, b, nx);
//	d_print_mat(nx, 1, x0, nx);

	// compute b0 = b + A*x0
	double *b0; d_zeros(&b0, nx, 1);
	dcopy_3l(nx, b, 1, b0, 1);
	dgemv_n_3l(nx, nx, A, nx, x0, b0);
//	d_print_mat(nx, 1, b, nx);
//	d_print_mat(nx, 1, b0, nx);

	// then A0 is a matrix of size 0x0
	double *A0; d_zeros(&A0, 0, 0);

/************************************************
* box constraints
************************************************/	

	int *idx0; i_zeros(&idx0, nbb[0], 1);
	double *lb0; d_zeros(&lb0, nbb[0], 1);
	double *ub0; d_zeros(&ub0, nbb[0], 1);
	for(jj=0; jj<nbu; jj++)
		{
		lb0[jj] = - 0.5;   //   umin
		ub0[jj] =   0.5;   //   umax
		idx0[jj] = jj;
		}
//	i_print_mat(nbb[0], 1, idx0, nbb[0]);

	int *idx1; i_zeros(&idx1, nbb[1], 1);
	double *lb1; d_zeros(&lb1, nbb[1], 1);
	double *ub1; d_zeros(&ub1, nbb[1], 1);
	for(jj=0; jj<nbu; jj++)
		{
		lb1[jj] = - 0.5;   //   umin
		ub1[jj] =   0.5;   //   umax
		idx1[jj] = jj;
		}
	for(; jj<nb; jj++)
		{
		lb1[jj] = - 4.0;   //   umin
		ub1[jj] =   4.0;   //   umax
		idx1[jj] = jj;
		}
//	i_print_mat(nbb[1], 1, idx1, nbb[1]);

	int *idxN; i_zeros(&idxN, nbb[N], 1);
	double *lbN; d_zeros(&lbN, nbb[N], 1);
	double *ubN; d_zeros(&ubN, nbb[N], 1);
	for(jj=0; jj<nbx; jj++)
		{
		lbN[jj] = - 4.0;   //   umin
		ubN[jj] =   4.0;   //   umax
		idxN[jj] = jj;
		}
//	i_print_mat(nbb[N], 1, idx1, nbb[N]);
	
/************************************************
* general constraints
************************************************/	

	double *C; d_zeros(&C, ng, nx);
	double *D; d_zeros(&D, ng, nu);
	double *lg; d_zeros(&lg, ng, 1);
	double *ug; d_zeros(&ug, ng, 1);

	double *CN; d_zeros(&CN, ngN, nx);
	for(ii=0; ii<ngN; ii++)
		CN[ii*(ngN+1)] = 1.0;
//	d_print_mat(ngN, nx, CN, ngN);
	double *lgN; d_zeros(&lgN, ngN, 1); // force all states to 0 at the last stage
	double *ugN; d_zeros(&ugN, ngN, 1); // force all states to 0 at the last stage
	

/************************************************
* cost function
************************************************/	
	
	double *Q; d_zeros(&Q, nx, nx);
	for(ii=0; ii<nx; ii++) Q[ii*(nx+1)] = 1.0;

	double *R; d_zeros(&R, nu, nu);
	for(ii=0; ii<nu; ii++) R[ii*(nu+1)] = 2.0;

	double *S; d_zeros(&S, nu, nx);

	double *q; d_zeros(&q, nx, 1);
	for(ii=0; ii<nx; ii++) q[ii] = 0.1;

	double *r; d_zeros(&r, nu, 1);
	for(ii=0; ii<nu; ii++) r[ii] = 0.2;

	// Q0 and q0 are matrices of size 0
	double *Q0; d_zeros(&Q0, 0, 0);
	double *q0; d_zeros(&q0, 0, 1);

	// compute r0 = r + S*x0
	double *r0; d_zeros(&r0, nu, 1);
	dcopy_3l(nu, r, 1, r0, 1);
	dgemv_n_3l(nu, nx, S, nu, x0, r0);

	// then S0 is a matrix of size nux0
	double *S0; d_zeros(&S0, nu, 0);

/************************************************
* problems data
************************************************/	

	double *hA[N];
	double *hB[N];
	double *hb[N];
	double *hQ[N];
	double *hS[N];
	double *hR[N];
	double *hq[N];
	double *hr[N];
	double *hlb[N];
	double *hub[N];
	int *hidx[N];
	double *hC[N];
	double *hD[N];
	double *hlg[N];
	double *hug[N];

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
	hidx[0] = idx0;
	hC[0] = C;
	hD[0] = D;
	hlg[0] = lg;
	hug[0] = ug;
	for(ii=1; ii<N; ii++)
		{
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
		hidx[ii] = idx1;
		hC[ii] = C;
		hD[ii] = D;
		hlg[ii] = lg;
		hug[ii] = ug;
		}
	hQ[N] = Q; // or maybe initialize to the solution of the DARE???
	hq[N] = q; // or maybe initialize to the solution of the DARE???
	hlb[N] = lbN;
	hub[N] = ubN;
	hidx[N] = idxN;
	hC[N] = CN;
	hlg[N] = lgN;
	hug[N] = ugN;

/************************************************
* solution
************************************************/	

	double *hx[N+1];
	double *hu[N+1];

	for(ii=0; ii<N; ii++)
		{
		d_zeros(&hx[ii], nxx[ii], 1);
		d_zeros(&hu[ii], nuu[ii], 1);
		}
	d_zeros(&hx[N], nxx[N], 1);

/************************************************
* solver arguments
************************************************/	

	// solver arguments
	struct ocp_qp_solver_args args;
	args.tol = 1e-6;
	args.max_iter = 10;
	args.min_step = 1e-8;
	args.mu0 = 0.0;
	args.sigma_min = 1e-3;

/************************************************
* work space
************************************************/	

//	int work_space_size = ocp_qp_hpmpc_workspace_double(N, nxx, nuu, nbb, ngg, args);

//	double *work; d_zeros_align(&work, work_space_size, 1); // aligned memory allocation

/************************************************
* call the solver
************************************************/	

	int return_value;


	struct timeval tv0, tv1, tv2;
	gettimeofday(&tv0, NULL); // stop

	for(rep=0; rep<nrep; rep++)
		{

		// call the QP OCP solver
//		return_value = ocp_qp_solver(N, nxx, nuu, nbb, ngg, hA, hB, hb, hQ, hS, hR, hq, hr, hidx, hlb, hub, hC, hD, hld, hud, hx, hu, args, work);	
		
		}

	gettimeofday(&tv1, NULL); // stop


	double time = (tv1.tv_sec-tv0.tv_sec)/(nrep+0.0)+(tv1.tv_usec-tv0.tv_usec)/(nrep*1e6);

	printf("\n");
	printf(" Average solution time over %d runs: %5.2e seconds\n", nrep, time);
	printf("\n\n");

/************************************************
* free memory
************************************************/	

	d_free(A);
	d_free(B);
	d_free(b);
	d_free(x0);
	d_free(A0);
	d_free(b0);
	d_free(Q);
	d_free(S);
	d_free(R);
	d_free(q);
	d_free(r);
	d_free(Q0);
	d_free(S0);
	d_free(q0);
	d_free(r0);
	i_free(idx0);
	d_free(lb0);
	d_free(ub0);
	i_free(idx1);
	d_free(lb1);
	d_free(ub1);
	i_free(idxN);
	d_free(lbN);
	d_free(ubN);
	d_free(C);
	d_free(D);
	d_free(lg);
	d_free(ug);
	d_free(CN);
	d_free(lgN);
	d_free(ugN);

	for(ii=0; ii<N; ii++)
		{
		d_free(hx[ii]);
		d_free(hu[ii]);
		}
	d_free(hx[N]);

//	d_free_align(work); // aligned memory free
	
	return 0;
	
	}


