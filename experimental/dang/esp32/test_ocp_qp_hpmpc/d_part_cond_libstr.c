/**************************************************************************************************
* acados/external/hpmpc/lqcp_solvers/d_part_cond_libstr.c                                                                                                *
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

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "aux_d.h"
#include "blas_d.h"
#include "block_size.h"
#include "lqcp_aux.h"
// 2017.03.13 Dang add target.h
#include "target.h"

#ifdef BLASFEO

#include <blasfeo_target.h>
#include <blasfeo_common.h>
#include <blasfeo_d_blas.h>
#include <blasfeo_d_aux.h>



void d_cond_BAbt_libstr(int N, int *nx, int *nu, struct d_strmat *hsBAbt, struct d_strmat *hsGamma, struct d_strmat *sBAbt2, void *work_space)
	{

	int ii, jj;

	struct d_strmat sA;

	int nu_tmp;

	nu_tmp = 0;
	ii = 0;
	// B & A & b
	dgecp_libstr(nu[0]+nx[0]+1, nx[1], 1.0, &hsBAbt[0], 0, 0, &hsGamma[0], 0, 0);
	//
	nu_tmp += nu[0];
	ii++;

	for(ii=1; ii<N; ii++)
		{
		// TODO check for equal pointers and avoid copy

		d_create_strmat(nx[ii+1], nx[ii], &sA, work_space);

		// pA in work space
		dgetr_libstr(nx[ii], nx[ii+1], 1.0, &hsBAbt[ii], nu[ii], 0, &sA, 0, 0); // pA in work // TODO avoid copy for LA_BLAS and LA_REFERENCE

		// Gamma * A^T
		dgemm_nt_libstr(nu_tmp+nx[0]+1, nx[ii+1], nx[ii], 1.0, &hsGamma[ii-1], 0, 0, &sA, 0, 0, 0.0, &hsGamma[ii], nu[ii], 0, &hsGamma[ii], nu[ii], 0); // Gamma * A^T

		dgecp_libstr(nu[ii], nx[ii+1], 1.0, &hsBAbt[ii], 0, 0, &hsGamma[ii], 0, 0);

		nu_tmp += nu[ii];

		dgead_libstr(1, nx[ii+1], 1.0, &hsBAbt[ii], nu[ii]+nx[ii], 0, &hsGamma[ii], nu_tmp+nx[0], 0);
		}

	// B & A & b
	dgecp_libstr(nu_tmp+nx[0]+1, nx[N], 1.0, &hsGamma[N-1], 0, 0, sBAbt2, 0, 0);

	return;

	}



void d_cond_RSQrq_libstr(int N, int *nx, int *nu, struct d_strmat *hsBAbt, struct d_strmat *hsRSQrq, struct d_strmat *hsGamma, struct d_strmat *sRSQrq2, void *work_space, int *work_space_sizes)
	{

	// early return
	if(N<1)
		return;

	int nn;

	int nu2[N+1];
	int nu3[N+1];
	nu2[0]= 0; // sum
	nu3[0]= 0; // reverse sum

	for(nn=0; nn<=N; nn++)
		{
		nu2[nn+1] = nu2[nn] + nu[nn];
		nu3[nn+1] = nu3[nn] + nu[N-nn-1];
		}

	struct d_strmat sL;
	struct d_strmat sM;
	struct d_strmat sLx;
	struct d_strmat sBAbtL;



	// early return
	if(N==1)
		{
		dgecp_libstr(nu[N-1]+nx[N-1]+1, nu[N-1]+nx[N-1], 1.0, &hsRSQrq[N-1], 0, 0, sRSQrq2, 0, 0);
		return;
		}



	char *c_ptr[4];
	c_ptr[0] = (char *) work_space;
	c_ptr[1] = c_ptr[0] + work_space_sizes[0];
	c_ptr[2] = c_ptr[1] + work_space_sizes[1];
	c_ptr[3] = c_ptr[2] + work_space_sizes[2];

	// final stage
	d_create_strmat(nu[N-1]+nx[N-1]+1, nu[N-1]+nx[N-1], &sL, (void *) c_ptr[0]);
	d_create_strmat(nu[N-1], nx[N-1], &sM, (void *) c_ptr[1]);

	dgecp_libstr(nu[N-1]+nx[N-1]+1, nu[N-1]+nx[N-1], 1.0, &hsRSQrq[N-1], 0, 0, &sL, 0, 0);

	// D
//	dgecp_libstr(nu[N-1], nu[N-1], 1.0, &sL, 0, 0, sRSQrq2, nu3[0], nu3[0]);
	dtrcp_l_libstr(nu[N-1], 1.0, &sL, 0, 0, sRSQrq2, nu3[0], nu3[0]);

	// M
	dgetr_libstr(nx[N-1], nu[N-1], 1.0, &sL, nu[N-1], 0, &sM, 0, 0);

	dgemm_nt_libstr(nu2[N-1]+nx[0]+1, nu[N-1], nx[N-1], 1.0, &hsGamma[N-2], 0, 0, &sM, 0, 0, 0.0, sRSQrq2, nu3[1], nu3[0], sRSQrq2, nu3[1], nu3[0]);

	// m
	dgead_libstr(1, nu[N-1], 1.0, &sL, nu[N-1]+nx[N-1], 0, sRSQrq2, nu2[N]+nx[0], nu3[0]);



	// middle stages
	for(nn=1; nn<N-1; nn++)
		{

		d_create_strmat(nx[N-nn]+1, nx[N-nn], &sLx, (void *) c_ptr[2]);
		d_create_strmat(nu[N-nn-1]+nx[N-nn-1]+1, nx[N-nn], &sBAbtL, (void *) c_ptr[3]);

		dgecp_libstr(nx[N-nn]+1, nx[N-nn], 1.0, &sL, nu[N-nn], nu[N-nn], &sLx, 0, 0);

		dpotrf_l_libstr(nx[N-nn]+1, nx[N-nn], &sLx, 0, 0, &sLx, 0, 0);

		dtrtr_l_libstr(nx[N-nn], 1.0, &sLx, 0, 0, &sLx, 0, 0);

		dtrmm_rutn_libstr(nu[N-nn-1]+nx[N-nn-1]+1, nx[N-nn], 1.0, &hsBAbt[N-nn-1], 0, 0, &sLx, 0, 0, 0.0, &sBAbtL, 0, 0, &sBAbtL, 0, 0);

		dgead_libstr(1, nx[N-nn], 1.0, &sLx, nx[N-nn], 0, &sBAbtL, nu[N-nn-1]+nx[N-nn-1], 0);

		d_create_strmat(nu[N-nn-1]+nx[N-nn-1]+1, nu[N-nn-1]+nx[N-nn-1], &sL, (void *) c_ptr[0]);
		d_create_strmat(nu[N-nn-1], nx[N-nn-1], &sM, (void *) c_ptr[1]);

		dsyrk_ln_libstr(nu[N-nn-1]+nx[N-nn-1]+1, nu[N-nn-1]+nx[N-nn-1], nx[N-nn], 1.0, &sBAbtL, 0, 0, &sBAbtL, 0, 0, 1.0, &hsRSQrq[N-nn-1], 0, 0, &sL, 0, 0);

		// D
//		dgecp_libstr(nu[N-nn-1], nu[N-nn-1], 1.0, &sL, 0, 0, sRSQrq2, nu3[nn], nu3[nn]);
		dtrcp_l_libstr(nu[N-nn-1], 1.0, &sL, 0, 0, sRSQrq2, nu3[nn], nu3[nn]);

		// M
		dgetr_libstr(nx[N-nn-1], nu[N-nn-1], 1.0, &sL, nu[N-nn-1], 0, &sM, 0, 0);

		dgemm_nt_libstr(nu2[N-nn-1]+nx[0]+1, nu[N-nn-1], nx[N-nn-1], 1.0, &hsGamma[N-nn-2], 0, 0, &sM, 0, 0, 0.0, sRSQrq2, nu3[nn+1], nu3[nn], sRSQrq2, nu3[nn+1], nu3[nn]);

		// m
		dgead_libstr(1, nu[N-nn-1], 1.0, &sL, nu[N-nn-1]+nx[N-nn-1], 0, sRSQrq2, nu2[N]+nx[0], nu3[nn]);

		}

	// first stage
	nn = N-1;

	d_create_strmat(nx[N-nn]+1, nx[N-nn], &sLx, (void *) c_ptr[2]);
	d_create_strmat(nu[N-nn-1]+nx[N-nn-1]+1, nx[N-nn], &sBAbtL, (void *) c_ptr[3]);

	dgecp_libstr(nx[N-nn]+1, nx[N-nn], 1.0, &sL, nu[N-nn], nu[N-nn], &sLx, 0, 0);

	dpotrf_l_libstr(nx[N-nn]+1, nx[N-nn], &sLx, 0, 0, &sLx, 0, 0);

	dtrtr_l_libstr(nx[N-nn], 1.0, &sLx, 0, 0, &sLx, 0, 0);

	dtrmm_rutn_libstr(nu[N-nn-1]+nx[N-nn-1]+1, nx[N-nn], 1.0, &hsBAbt[N-nn-1], 0, 0, &sLx, 0, 0, 0.0, &sBAbtL, 0, 0, &sBAbtL, 0, 0);

	dgead_libstr(1, nx[N-nn], 1.0, &sLx, nx[N-nn], 0, &sBAbtL, nu[N-nn-1]+nx[N-nn-1], 0);

	d_create_strmat(nu[N-nn-1]+nx[N-nn-1]+1, nu[N-nn-1]+nx[N-nn-1], &sL, (void *) c_ptr[0]);

	dsyrk_ln_libstr(nu[N-nn-1]+nx[N-nn-1]+1, nu[N-nn-1]+nx[N-nn-1], nx[N-nn], 1.0, &sBAbtL, 0, 0, &sBAbtL, 0, 0, 1.0, &hsRSQrq[N-nn-1], 0, 0, &sL, 0, 0);

	// D, M, m, P, p
//	dgecp_libstr(nu[0]+nx[0]+1, nu[0]+nx[0], 1.0, &sL, 0, 0, sRSQrq2, nu3[N-1], nu3[N-1]); // TODO dtrcp for 'rectangular' matrices
	dtrcp_l_libstr(nu[0]+nx[0], 1.0, &sL, 0, 0, sRSQrq2, nu3[N-1], nu3[N-1]); // TODO dtrcp for 'rectangular' matrices
	dgecp_libstr(1, nu[0]+nx[0], 1.0, &sL, nu[0]+nx[0], 0, sRSQrq2, nu3[N-1]+nu[0]+nx[0], nu3[N-1]); // TODO dtrcp for 'rectangular' matrices

	return;

	}



// TODO general constraints !!!!!
void d_cond_DCtd_libstr(int N, int *nx, int *nu, int *nb, int **hidxb, int *ng, struct d_strmat *hsDCt, struct d_strvec *hsd, struct d_strmat *hsGamma, struct d_strmat *sDCt2, struct d_strvec *sd2, int *idxb2, void *work_space)
	{

	// early return
	if(N<1)
		return;

	double *d2 = sd2->pa;
	double *ptr_d;

	int nu_tmp, ng_tmp;

	int ii, jj;

	// problem size

	int nbb = nb[0]; // box that remain box constraints
	int nbg = 0; // box that becomes general constraints
	for(ii=1; ii<N; ii++)
		for(jj=0; jj<nb[ii]; jj++)
			if(hidxb[ii][jj]<nu[ii])
				nbb++;
			else
				nbg++;

	int nx2 = nx[0];
	int nu2 = nu[0];
	int ngg = ng[0];
	for(ii=1; ii<N; ii++)
		{
		nu2 += nu[ii];
		ngg += ng[ii];
		}
	int ng2 = nbg + ngg;
	int nb2 = nbb;

	// set constraint matrix to zero (it's 2 lower triangular matrices atm)
	dmatse_libstr(nu2+nx2, ng2, 0.0, sDCt2, 0, 0);

	// box constraints

	int idx_gammab = nx[0];
	for(ii=0; ii<N-1; ii++)
		idx_gammab += nu[ii];

	int ib = 0;
	int ig = 0;

	double tmp;
	int idx_g;

	// middle stages
	nu_tmp = 0;
	for(ii=0; ii<N-1; ii++)
		{
		nu_tmp += nu[N-1-ii];
		ptr_d = hsd[N-1-ii].pa;
		for(jj=0; jj<nb[N-1-ii]; jj++)
			{
			if(hidxb[N-1-ii][jj]<nu[N-1-ii]) // input: box constraint
				{
				d2[0*nbb+ib] = ptr_d[0*nb[N-1-ii]+jj];
				d2[1*nbb+ib] = ptr_d[1*nb[N-1-ii]+jj];
				idxb2[ib] = nu_tmp - nu[N-1-ii] + hidxb[N-1-ii][jj];
				ib++;
				}
			else // state: general constraint
				{
				idx_g = hidxb[N-1-ii][jj]-nu[N-1-ii];
				tmp = dmatex1_libstr(&hsGamma[N-2-ii], idx_gammab, idx_g);
				d2[2*nbb+0*ng2+ig] = ptr_d[0*nb[N-1-ii]+jj] - tmp;
				d2[2*nbb+1*ng2+ig] = ptr_d[1*nb[N-1-ii]+jj] - tmp;
				dgecp_libstr(idx_gammab, 1, 1.0, &hsGamma[N-ii-2], 0, idx_g, sDCt2, nu_tmp, ig);
				ig++;
				}
			}
		idx_gammab -= nu[N-2-ii];
		}

	// initial stage: both inputs and states as box constraints
	nu_tmp += nu[0];
	ptr_d = hsd[0].pa;
	for(jj=0; jj<nb[0]; jj++)
		{
		d2[0*nbb+ib] = ptr_d[0*nb[0]+jj];
		d2[1*nbb+ib] = ptr_d[1*nb[0]+jj];
		idxb2[ib] = nu_tmp - nu[0] + hidxb[0][jj];
		ib++;
		}

	// XXX for now, just shift after box-to-general constraints
	// better interleave them, to keep the block lower trianlgular structure !!!

	// general constraints

	struct d_strmat sC;
	struct d_strvec sGammax;
	struct d_strvec sCGammax;

	char *c_ptr;

	nu_tmp = 0;
	ng_tmp = 0;
	for(ii=0; ii<N-1; ii++)
		{

		c_ptr = (char *) work_space;

		dgecp_libstr(nu[N-1-ii], ng[N-1-ii], 1.0, &hsDCt[N-1-ii], 0, 0, sDCt2, nu_tmp, nbg+ng_tmp);

		nu_tmp += nu[N-1-ii];

		d_create_strmat(ng[N-1-ii], nx[N-1-ii], &sC, (void *) c_ptr);
		c_ptr += sC.memory_size;

		dgetr_libstr(nx[N-1-ii], ng[N-1-ii], 1.0, &hsDCt[N-1-ii], nu[N-1-ii], 0, &sC, 0, 0);

		dgemm_nt_libstr(nu2+nx[0]-nu_tmp, ng[N-1-ii], nx[N-1-ii], 1.0, &hsGamma[N-2-ii], 0, 0, &sC, 0, 0, 0.0, sDCt2, nu_tmp, nbg+ng_tmp, sDCt2, nu_tmp, nbg+ng_tmp);

		dveccp_libstr(ng[N-1-ii], 1.0, &hsd[N-1-ii], 2*nb[N-1-ii]+0*ng[N-1-ii], sd2, 2*nb2+0*ng2+nbg+ng_tmp);
		dveccp_libstr(ng[N-1-ii], 1.0, &hsd[N-1-ii], 2*nb[N-1-ii]+1*ng[N-1-ii], sd2, 2*nb2+1*ng2+nbg+ng_tmp);

		d_create_strvec(nx[N-1-ii], &sGammax, (void *) c_ptr);
		c_ptr += sGammax.memory_size;
		d_create_strvec(ng[N-1-ii], &sCGammax, (void *) c_ptr);
		c_ptr += sCGammax.memory_size;

		drowex_libstr(nx[N-1-ii], 1.0, &hsGamma[N-2-ii], nu2+nx[0]-nu_tmp, 0, &sGammax, 0);

		dgemv_n_libstr(ng[N-1-ii], nx[N-1-ii], 1.0, &sC, 0, 0, &sGammax, 0, 0.0, &sCGammax, 0, &sCGammax, 0);

		daxpy_libstr(ng[N-1-ii], -1.0, &sCGammax, 0, sd2, 2*nb2+nbg+ng_tmp);
		daxpy_libstr(ng[N-1-ii], -1.0, &sCGammax, 0, sd2, 2*nb2+ng2+nbg+ng_tmp);

		ng_tmp += ng[N-1-ii];

		}

	ii = N-1;

	dgecp_libstr(nu[0]+nx[0], ng[0], 1.0, &hsDCt[0], 0, 0, sDCt2, nu_tmp, nbg+ng_tmp);

	dveccp_libstr(ng[0], 1.0, &hsd[0], 2*nb[0]+0, sd2, 2*nb2+nbg+ng_tmp);
	dveccp_libstr(ng[0], 1.0, &hsd[0], 2*nb[0]+ng[0], sd2, 2*nb2+ng2+nbg+ng_tmp);

//	ng_tmp += ng[N-1-ii];

	return;

	}



// XXX does not compute hidxb2, since nb2 has to be known to allocate the right space for hidxb2 !!!
void d_part_cond_compute_problem_size_libstr(int N, int *nx, int *nu, int *nb, int **hidxb, int *ng, int N2, int *nx2, int *nu2, int *nb2, int *ng2)
	{

	int ii, jj, kk;

	int N1 = N/N2; // (floor) horizon of small blocks
	int R1 = N - N2*N1; // the first R1 blocks have horizon N1+1
	int M1 = R1>0 ? N1+1 : N1; // (ceil) horizon of large blocks
	int T1; // horizon of current block

	int N_tmp = 0; // temporary sum of horizons
	int nbb; // box constr that remain box constr
	int nbg; // box constr that becomes general constr
	for(ii=0; ii<N2; ii++)
		{
		T1 = ii<R1 ? M1 : N1;
		nx2[ii] = nx[N_tmp+0];
		nu2[ii] = nu[N_tmp+0];
		nb2[ii] = nb[N_tmp+0];
		ng2[ii] = ng[N_tmp+0];
		for(jj=1; jj<T1; jj++)
			{
			nbb = 0;
			nbg = 0;
			for(kk=0; kk<nb[N_tmp+jj]; kk++)
				if(hidxb[N_tmp+jj][kk]<nu[N_tmp+jj])
					nbb++;
				else
					nbg++;
			nx2[ii] += 0;
			nu2[ii] += nu[N_tmp+jj];
			nb2[ii] += nbb;
			ng2[ii] += ng[N_tmp+jj] + nbg;
			}
		N_tmp += T1;
		}
	nx2[N2] = nx[N];
	nu2[N2] = nu[N];
	nb2[N2] = nb[N];
	ng2[N2] = ng[N];

	}



int d_part_cond_work_space_size_bytes_libstr(int N, int *nx, int *nu, int *nb, int **hidxb, int *ng, int N2, int *nx2, int *nu2, int *nb2, int *ng2, int *work_space_sizes)
	{

	int ii, jj;
	int nu_tmp;

	int N1 = N/N2; // (floor) horizon of small blocks
	int R1 = N - N2*N1; // the first R1 blocks have horizon N1+1
	int M1 = R1>0 ? N1+1 : N1; // (ceil) horizon of large blocks
	int T1; // horizon of current block
	int N_tmp; // temporary sum of horizons

	int Gamma_size = 0;
	int pA_size = 0;
	int pC_size = 0;
	work_space_sizes[0] = 0;
	work_space_sizes[1] = 0;
	work_space_sizes[2] = 0;
	work_space_sizes[3] = 0;

	int tmp_size;

	N_tmp = 0;
	for(ii=0; ii<N2; ii++)
		{

		T1 = ii<R1 ? M1 : N1;

		// hpGamma
		nu_tmp = 0;
		tmp_size = 0;
		for(jj=0; jj<T1; jj++)
			{
			nu_tmp += nu[N_tmp+jj];
			tmp_size += d_size_strmat(nx[N_tmp]+nu_tmp+1, nx[N_tmp+jj+1]);
			}
		Gamma_size = tmp_size>Gamma_size ? tmp_size : Gamma_size;

		// sA
		for(jj=0; jj<T1; jj++)
			{
			tmp_size = d_size_strmat(nx[N_tmp+jj+1], nx[N_tmp+jj]);
			pA_size = tmp_size > pA_size ? tmp_size : pA_size;
			}

		// sC
		for(jj=0; jj<T1; jj++)
			{
			tmp_size = d_size_strmat(ng[N_tmp+jj], nx[N_tmp+jj]);
			tmp_size += d_size_strvec(nx[N_tmp+jj]);
			tmp_size += d_size_strvec(ng[N_tmp+jj]);
			pC_size = tmp_size > pC_size ? tmp_size : pC_size;
			}

		// sL : 0 => N-1
		for(jj=0; jj<T1; jj++)
			{
			tmp_size = d_size_strmat(nu[N_tmp+jj]+nx[N_tmp+jj]+1, nu[N_tmp+jj]+nx[N_tmp+jj]);
			work_space_sizes[0] = tmp_size>work_space_sizes[0] ? tmp_size : work_space_sizes[0];
			}

		// sM : 1 => N-1
		for(jj=1; jj<T1; jj++)
			{
			tmp_size = d_size_strmat(nu[N_tmp+jj], nx[N_tmp+jj]);
			work_space_sizes[1] = tmp_size>work_space_sizes[1] ? tmp_size : work_space_sizes[1];
			}

		// sLx : 1 => N-1
		for(jj=1; jj<T1; jj++)
			{
			tmp_size = d_size_strmat(nx[N_tmp+jj]+1, nx[N_tmp+jj]);
			work_space_sizes[2] = tmp_size>work_space_sizes[2] ? tmp_size : work_space_sizes[2];
			}

		// sBAbtL : 0 => N-2
		for(jj=0; jj<T1-1; jj++)
			{
			tmp_size = d_size_strmat(nu[N_tmp+jj]+nx[N_tmp+jj]+1, nx[N_tmp+1+jj]);
			work_space_sizes[3] = tmp_size>work_space_sizes[3] ? tmp_size : work_space_sizes[3];
			}

		N_tmp += T1;

		}

	tmp_size = work_space_sizes[0] + work_space_sizes[1] + work_space_sizes[2] + work_space_sizes[3];
	tmp_size = tmp_size>pA_size ? tmp_size : pA_size;
	tmp_size = tmp_size>pC_size ? tmp_size : pC_size;
	int size = Gamma_size+tmp_size;

	size = (size + 63) / 64 * 64; // make work space multiple of (typical) cache line size

	return size;

	}





int d_part_cond_memory_space_size_bytes_libstr(int N, int *nx, int *nu, int *nb, int **hidxb, int *ng, int N2, int *nx2, int *nu2, int *nb2, int *ng2)
	{

	// early return
	if(N2==N)
		{
		return 0;
		}

	int ii, jj, kk;

	// data matrices
	int size = 0;
	int i_size = 0;
	for(ii=0; ii<N2; ii++)
		{
		// hpBAbt2
		size += d_size_strmat(nu2[ii]+nx2[ii]+1, nx2[ii+1]);
		// hpRSQrq2
		size += d_size_strmat(nu2[ii]+nx2[ii]+1, nu2[ii]+nx2[ii]);
		// hDCt2
		size += d_size_strmat(nu2[ii]+nx2[ii], ng2[ii]);
		// hd2
		size += d_size_strvec(2*nb2[ii]+2*ng2[ii]);
		// hidxb2
		i_size += nb2[ii];
		}

	size = size + i_size*sizeof(int);

	// make memory space multiple of (typical) cache line size
	size = (size + 63) / 64 * 64;

	return size;

	}





void d_part_cond_libstr(int N, int *nx, int *nu, int *nb, int **hidxb, int *ng, struct d_strmat *hsBAbt, struct d_strmat *hsRSQrq, struct d_strmat *hsDCt, struct d_strvec *hsd, int N2, int *nx2, int *nu2, int *nb2, int **hidxb2, int *ng2, struct d_strmat *hsBAbt2, struct d_strmat *hsRSQrq2, struct d_strmat *hsDCt2, struct d_strvec *hsd2, void *memory, void *work, int *work_space_sizes)
	{

	int ii, jj, kk;
	int nu_tmp;

	// early return
	if(N2==N)
		{
		for(ii=0; ii<N; ii++)
			{
			nx2[ii] = nx[ii];
			nu2[ii] = nu[ii];
			nb2[ii] = nb[ii];
			hidxb2[ii] = hidxb[ii];
			ng2[ii] = ng[ii];
			hsBAbt2[ii] = hsBAbt[ii];
			hsRSQrq2[ii] = hsRSQrq[ii];
			hsDCt2[ii] = hsDCt[ii];
			hsd2[ii] = hsd[ii];
			}
		ii = N;
		nx2[ii] = nx[ii];
		nu2[ii] = nu[ii];
		nb2[ii] = nb[ii];
		hidxb2[ii] = hidxb[ii];
		ng2[ii] = ng[ii];
		hsRSQrq2[ii] = hsRSQrq[ii];
		hsDCt2[ii] = hsDCt[ii];
		hsd2[ii] = hsd[ii];
		return;
		}

	// sequential update not implemented
	if(N2>N)
		{
		printf("\nError: it must be N2<=N, sequential update not implemented\n\n");
		exit(1);
		}

	// general constraints not implemented (can be only ng[N]>0)
//	for(ii=0; ii<N; ii++)
//		if(ng[ii]>0)
//			{
//			printf("\nError: it must be ng>0, general constraints case not implemented\n\n");
//			exit(1);
//			}

	int N1 = N/N2; // (floor) horizon of small blocks
	int R1 = N - N2*N1; // the first R1 blocks have horizon N1+1
	int M1 = R1>0 ? N1+1 : N1; // (ceil) horizon of large blocks
	int T1; // horizon of current block
	int N_tmp = 0; // temporary sum of horizons

	// data matrices (memory space) no last stage !!!!!
	char *c_ptr = (char *) memory;
	for(ii=0; ii<N2; ii++)
		{
		d_create_strmat(nu2[ii]+nx2[ii]+1, nx2[ii+1], &hsBAbt2[ii], (void *) c_ptr);
		c_ptr += hsBAbt2[ii].memory_size;
		}
	for(ii=0; ii<N2; ii++)
		{
		d_create_strmat(nu2[ii]+nx2[ii]+1, nu2[ii]+nx2[ii], &hsRSQrq2[ii], (void *) c_ptr);
		c_ptr += hsRSQrq2[ii].memory_size;
		}
	for(ii=0; ii<N2; ii++)
		{
		d_create_strmat(nu2[ii]+nx2[ii], ng2[ii], &hsDCt2[ii], (void *) c_ptr);
		c_ptr += hsDCt2[ii].memory_size;
		}
	for(ii=0; ii<N2; ii++)
		{
		d_create_strvec(2*nb2[ii]+2*ng2[ii], &hsd2[ii], (void *) c_ptr);
		c_ptr += hsd2[ii].memory_size;
		}
	int *i_ptr = (int *) c_ptr;
	for(ii=0; ii<N2; ii++)
		{
		hidxb2[ii] = i_ptr;
		i_ptr += nb2[ii];
		}

	// work space
	struct d_strmat hsGamma[M1];
	int tmp_i, tmp_size;

	// other stages
	N_tmp = 0;
	for(ii=0; ii<N2; ii++)
		{

		T1 = ii<R1 ? M1 : N1;

		c_ptr = (char *) work;

		// sGamma
		nu_tmp = nu[N_tmp+0];
		for(jj=0; jj<T1; jj++)
			{
			d_create_strmat(nx[N_tmp+0]+nu_tmp+1, nx[N_tmp+jj+1], &hsGamma[jj], (void *) c_ptr);
			c_ptr += hsGamma[jj].memory_size;
			nu_tmp += nu[N_tmp+jj+1];
			}
		// sA
		// no c_ptr += ... : overwrite sA
		d_cond_BAbt_libstr(T1, &nx[N_tmp], &nu[N_tmp], &hsBAbt[N_tmp], hsGamma, &hsBAbt2[ii], (void *) c_ptr);

		d_cond_RSQrq_libstr(T1, &nx[N_tmp], &nu[N_tmp], &hsBAbt[N_tmp], &hsRSQrq[N_tmp], hsGamma, &hsRSQrq2[ii], (void *) c_ptr, work_space_sizes);

		d_cond_DCtd_libstr(T1, &nx[N_tmp], &nu[N_tmp], &nb[N_tmp], &hidxb[N_tmp], &ng[N_tmp], &hsDCt[N_tmp], &hsd[N_tmp], hsGamma, &hsDCt2[ii], &hsd2[ii], hidxb2[ii], (void *) c_ptr);
		N_tmp += T1;
//exit(1);
		}

	// last stage
	hsRSQrq2[N2] = hsRSQrq[N];
	hsDCt2[N2] = hsDCt[N];
	hsd2[N2] = hsd[N];
	hidxb2[N2] = hidxb[N];

	return;

	}



int d_part_expand_work_space_size_bytes_libstr(int N, int *nx, int *nu, int *nb, int *ng, int *work_space_sizes)
	{

	int ii, i_tmp;

	int nzM = d_size_strvec(nu[0]+nx[0]);
	int ngM = d_size_strvec(ng[0]);

	for(ii=1; ii<=N; ii++)
		{
		i_tmp = d_size_strvec(nu[ii]+nx[ii]);
		nzM = i_tmp>nzM ? i_tmp : nzM;
		i_tmp = d_size_strvec(ng[ii]);
		ngM = i_tmp>ngM ? i_tmp : ngM;
		}

	work_space_sizes[0] = nzM;
	work_space_sizes[1] = ngM;

	int size = work_space_sizes[0] + work_space_sizes[1];

	size = (size + 63) / 64 * 64; // make multiple of (typical) cache line size

	return size;

	}



void d_part_expand_solution_libstr(int N, int *nx, int *nu, int *nb, int **hidxb, int *ng, struct d_strmat *hsBAbt, struct d_strvec *hsb, struct d_strmat *hsRSQrq, struct d_strvec *hsrq, struct d_strmat *hsDCt, struct d_strvec *hsux, struct d_strvec *hspi, struct d_strvec *hslam, struct d_strvec *hst, int N2, int *nx2, int *nu2, int *nb2, int **hidxb2, int *ng2, struct d_strvec *hsux2, struct d_strvec *hspi2, struct d_strvec *hslam2, struct d_strvec *hst2, void *work_space, int *work_space_sizes)
	{

	int ii, jj, ll;

	struct d_strvec workvec[2];

	double *ptr_work0, *ptr_work1, *ptr_lam, *ptr_t, *ptr_lam2, *ptr_t2;

	char *c_ptr[2];
	c_ptr[0] = (char *) work_space;
	c_ptr[1] = c_ptr[0] + work_space_sizes[0];

	int N1 = N/N2; // (floor) horizon of small blocks
	int R1 = N - N2*N1; // the first R1 blocks have horizion N1+1
	int M1 = R1>0 ? N1+1 : N1; // (ceil) horizon of large blocks
	int T1; // horizon of current block
	int N_tmp, nu_tmp;
	int nbb2, nbg2, ngg2;
	int stg;

	// inputs & initial states
	N_tmp = 0;
	for(ii=0; ii<N2; ii++)
		{
		T1 = ii<R1 ? M1 : N1;
		nu_tmp = 0;
		// final stages: copy only input
		for(jj=0; jj<T1-1; jj++)
			{
			dveccp_libstr(nu[N_tmp+T1-1-jj], 1.0, &hsux2[ii], nu_tmp, &hsux[N_tmp+T1-1-jj], 0);
			nu_tmp += nu[N_tmp+T1-1-jj];
			}
		// first stage: copy input and state
		dveccp_libstr(nu[N_tmp+0]+nx[N_tmp+0], 1.0, &hsux2[ii], nu_tmp, &hsux[N_tmp+0], 0);
		//
		N_tmp += T1;
		}

	// copy final state
	dveccp_libstr(nx[N], 1.0, &hsux2[N2], 0, &hsux[N], 0);

	// compute missing states by simulation within each block
	N_tmp = 0;
	for(ii=0; ii<N2; ii++)
		{
		T1 = ii<R1 ? M1 : N1;
		for(jj=0; jj<T1-1; jj++) // last stage is already there !!!
			{
			dgemv_t_libstr(nu[N_tmp+jj]+nx[N_tmp+jj], nx[N_tmp+jj+1], 1.0, &hsBAbt[N_tmp+jj], 0, 0, &hsux[N_tmp+jj], 0, 1.0, &hsb[N_tmp+jj], 0, &hsux[N_tmp+jj+1], nu[N_tmp+jj+1]);
			}
		//
		N_tmp += T1;
		}

	// slack variables and ineq lagrange multipliers
	N_tmp = 0;
	for(ii=0; ii<N2; ii++)
		{
		ptr_lam2 = hslam2[ii].pa;
		ptr_t2 = hst2[ii].pa;
		nbb2 = 0;
		nbg2 = 0;
		ngg2 = 0;
		T1 = ii<R1 ? M1 : N1;
		// final stages
		for(jj=0; jj<T1-1; jj++)
			{
			stg = N_tmp+T1-1-jj;
			ptr_lam = hslam[stg].pa;
			ptr_t = hst[stg].pa;
			for(ll=0; ll<nb[stg]; ll++)
				{
				if(hidxb[stg][ll]<nu[stg])
					{
					// box as box
					ptr_lam[0*nb[stg]+ll] = ptr_lam2[0*nb2[ii]+nbb2];
					ptr_lam[1*nb[stg]+ll] = ptr_lam2[1*nb2[ii]+nbb2];
					ptr_t[0*nb[stg]+ll] = ptr_t2[0*nb2[ii]+nbb2];
					ptr_t[1*nb[stg]+ll] = ptr_t2[1*nb2[ii]+nbb2];
					nbb2++;
					}
				else
					{
					// box as general XXX change when decide where nbg are placed wrt ng
					ptr_lam[0*nb[stg]+ll] = ptr_lam2[2*nb2[ii]+0*ng2[ii]+nbg2];
					ptr_lam[1*nb[stg]+ll] = ptr_lam2[2*nb2[ii]+1*ng2[ii]+nbg2];
					ptr_t[0*nb[stg]+ll] = ptr_t2[2*nb2[ii]+0*ng2[ii]+nbg2];
					ptr_t[1*nb[stg]+ll] = ptr_t2[2*nb2[ii]+1*ng2[ii]+nbg2];
					nbg2++;
					}
				}
			}
		for(jj=0; jj<T1-1; jj++)
			{
			stg = N_tmp+T1-1-jj;
			ptr_lam = hslam[stg].pa;
			ptr_t = hst[stg].pa;
			for(ll=0; ll<ng[stg]; ll++)
				{
				// general as general
				ptr_lam[2*nb[stg]+0*ng[stg]+ll] = ptr_lam2[2*nb2[ii]+0*ng2[ii]+nbg2+ngg2];
				ptr_lam[2*nb[stg]+1*ng[stg]+ll] = ptr_lam2[2*nb2[ii]+1*ng2[ii]+nbg2+ngg2];
				ptr_t[2*nb[stg]+0*ng[stg]+ll] = ptr_t2[2*nb2[ii]+0*ng2[ii]+nbg2+ngg2];
				ptr_t[2*nb[stg]+1*ng[stg]+ll] = ptr_t2[2*nb2[ii]+1*ng2[ii]+nbg2+ngg2];
				ngg2++;
				}
			}
		// first stage
		stg = N_tmp+T1-1-jj;
		// all box as box
		dveccp_libstr(nb[N_tmp+0], 1.0, &hslam2[ii], 0*nb2[ii]+nbb2, &hslam[stg], 0*nb[stg]);
		dveccp_libstr(nb[N_tmp+0], 1.0, &hslam2[ii], 1*nb2[ii]+nbb2, &hslam[stg], 1*nb[stg]);
		dveccp_libstr(nb[N_tmp+0], 1.0, &hst2[ii], 0*nb2[ii]+nbb2, &hst[stg], 0*nb[stg]);
		dveccp_libstr(nb[N_tmp+0], 1.0, &hst2[ii], 1*nb2[ii]+nbb2, &hst[stg], 1*nb[stg]);
		// first stage: general
		dveccp_libstr(ng[N_tmp+0], 1.0, &hslam2[ii], 2*nb2[ii]+0*ng2[ii]+nbg2+ngg2, &hslam[stg], 2*nb[stg]+0*ng[stg]);
		dveccp_libstr(ng[N_tmp+0], 1.0, &hslam2[ii], 2*nb2[ii]+1*ng2[ii]+nbg2+ngg2, &hslam[stg], 2*nb[stg]+1*ng[stg]);
		dveccp_libstr(ng[N_tmp+0], 1.0, &hst2[ii], 2*nb2[ii]+0*ng2[ii]+nbg2+ngg2, &hst[stg], 2*nb[stg]+0*ng[stg]);
		dveccp_libstr(ng[N_tmp+0], 1.0, &hst2[ii], 2*nb2[ii]+1*ng2[ii]+nbg2+ngg2, &hst[stg], 2*nb[stg]+1*ng[stg]);
		//
		N_tmp += T1;
		}
	// last stage: just copy
	// box
	dveccp_libstr(nb[N], 1.0, &hslam2[N2], 0*nb2[N2], &hslam[N], 0*nb[N]);
	dveccp_libstr(nb[N], 1.0, &hslam2[N2], 1*nb2[N2], &hslam[N], 1*nb[N]);
	dveccp_libstr(nb[N], 1.0, &hst2[N2], 0*nb2[N2], &hst[N], 0*nb[N]);
	dveccp_libstr(nb[N], 1.0, &hst2[N2], 1*nb2[N2], &hst[N], 1*nb[N]);
	// general
	dveccp_libstr(ng[N], 1.0, &hslam2[N2], 2*nb2[N2]+0*ng2[N2], &hslam[N], 2*nb[N]+0*ng[N]);
	dveccp_libstr(ng[N], 1.0, &hslam2[N2], 2*nb2[N2]+1*ng2[N2], &hslam[N], 2*nb[N]+1*ng[N]);
	dveccp_libstr(ng[N], 1.0, &hst2[N2], 2*nb2[N2]+0*ng2[N2], &hst[N], 2*nb[N]+0*ng[N]);
	dveccp_libstr(ng[N], 1.0, &hst2[N2], 2*nb2[N2]+1*ng2[N2], &hst[N], 2*nb[N]+1*ng[N]);

	// lagrange multipliers of equality constraints
	N_tmp = 0;
	for(ii=0; ii<N2; ii++)
		{
		T1 = ii<R1 ? M1 : N1;
		// last stage: just copy
		dveccp_libstr(nx[N_tmp+T1], 1.0, &hspi2[ii+1], 0, &hspi[N_tmp+T1], 0);
		// middle stages: backward simulation
		for(jj=0; jj<T1-1; jj++)
			{
			stg = N_tmp+T1-1-jj;
			d_create_strvec(nu[stg]+nx[stg], &workvec[0], (void *) c_ptr[0]);
			dveccp_libstr(nu[stg]+nx[stg], 1.0, &hsrq[stg], 0, &workvec[0], 0);
			ptr_work0 = workvec[0].pa;
			ptr_lam = hslam[stg].pa;
			for(ll=0; ll<nb[stg]; ll++)
				ptr_work0[hidxb[stg][ll]] += - ptr_lam[0*nb[stg]+ll] + ptr_lam[1*nb[stg]+ll];
			dsymv_l_libstr(nu[stg]+nx[stg], nu[stg]+nx[stg], 1.0, &hsRSQrq[stg], 0, 0, &hsux[stg], 0, 1.0, &workvec[0], 0, &workvec[0], 0);
			dgemv_n_libstr(nu[stg]+nx[stg], nx[stg+1], 1.0, &hsBAbt[stg], 0, 0, &hspi[stg+1], 0, 1.0, &workvec[0], 0, &workvec[0], 0);
			d_create_strvec(ng[stg], &workvec[1], (void *) c_ptr[1]);
			ptr_work1 = workvec[1].pa;
			for(ll=0; ll<ng[stg]; ll++)
				ptr_work1[ll] = ptr_lam[2*nb[stg]+1*ng[stg]+ll] - ptr_lam[2*nb[stg]+0*ng[stg]+ll];
			dgemv_n_libstr(nu[stg]+nx[stg], ng[stg], 1.0, &hsDCt[stg], 0, 0, &workvec[1], 0, 1.0, &workvec[0], 0, &workvec[0], 0);
			dveccp_libstr(nx[stg], 1.0, &workvec[0], nu[stg],  &hspi[stg], 0);
			}
		N_tmp += T1;
		}

	return;

	}

#endif
