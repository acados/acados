/**************************************************************************************************
* acados/external/hpmpc/lqcp_solvers/d_back_ric_rec_libstr.c                                                                                                *
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

#include <stdlib.h>
#include <stdio.h>
// 2017.03.13 Dang add target.h
#include "target.h"

#ifdef BLASFEO

#include <blasfeo_target.h>
#include <blasfeo_common.h>
#include <blasfeo_i_aux.h>
#include <blasfeo_d_aux.h>
#include <blasfeo_d_kernel.h>
#include <blasfeo_d_blas.h>



int d_back_ric_rec_work_space_size_bytes_libstr(int N, int *nx, int *nu, int *nb, int *ng)
	{

	int ii;

	// max sizes
	int nxM  = 0;
	int ngM = 0;
	int nuxM  = 0;
	int nxgM = ng[N];
	for(ii=0; ii<N; ii++)
		{
		nxM = nx[ii]>nxM ? nx[ii] : nxM;
		ngM = ng[ii]>ngM ? ng[ii] : ngM;
		nuxM = nu[ii]+nx[ii]>nuxM ? nu[ii]+nx[ii]+1 : nuxM;
		nxgM = nx[ii+1]+ng[ii]>nxgM ? nx[ii+1]+ng[ii] : nxgM;
		}
	ii = N;
	nxM = nx[ii]>nxM ? nx[ii] : nxM;
	ngM = ng[ii]>ngM ? ng[ii] : ngM;
	nuxM = nu[ii]+nx[ii]>nuxM ? nu[ii]+nx[ii]+1 : nuxM;

	int size = 0;

	size += d_size_strmat(nuxM+1, nxgM); // ric_work_mat[0]
	if(ngM>0)
		size += d_size_strmat(nuxM, nxgM); // ric_work_mat[1]
	size += d_size_strvec(nxM); // ric_work_vec[0]

	// make multiple of (typical) cache line size
	size = (size+63)/64*64;

	return size;
	}



void d_back_ric_rec_sv_libstr(int N, int *nx, int *nu, int *nb, int **hidxb, int *ng, int update_b, struct d_strmat *hsBAbt, struct d_strvec *hsb, int update_q, struct d_strmat *hsRSQrq, struct d_strvec *hsrq, struct d_strmat *hsDCt, struct d_strvec *hsQx, struct d_strvec *hsqx, struct d_strvec *hsux, int compute_pi, struct d_strvec *hspi, int compute_Pb, struct d_strvec *hsPb, struct d_strmat *hsL, struct d_strmat *hsLxt, void *work)
	{

	char *c_ptr;

	struct d_strmat hswork_mat_0, hswork_mat_1;
	struct d_strvec hswork_vec_0;

	int nn;

	// factorization and backward substitution

	// last stage
	dtrcp_l_libstr(nu[N]+nx[N], 1.0, &hsRSQrq[N], 0, 0, &hsL[N], 0, 0); // TODO dtrcp_l_libstr with m and n, for m>=n
	if(update_q)
		{
		drowin_libstr(nu[N]+nx[N], 1.0, &hsrq[N], 0, &hsL[N], nu[N]+nx[N], 0);
		}
	else
		{
		dgecp_libstr(1, nu[N]+nx[N], 1.0, &hsRSQrq[N], nu[N]+nx[N], 0, &hsL[N], nu[N]+nx[N], 0);
		}
	if(nb[N]>0)
		{
		ddiaad_libspstr(nb[N], hidxb[N], 1.0, &hsQx[N], 0, &hsL[N], 0, 0);
		drowad_libspstr(nb[N], hidxb[N], 1.0, &hsqx[N], 0, &hsL[N], nu[N]+nx[N], 0);
		}
	if(ng[N]>0)
		{
		c_ptr = (char *) work;
		d_create_strmat(nu[N]+nx[N]+1, ng[N], &hswork_mat_0, (void *) c_ptr);
		c_ptr += hswork_mat_0.memory_size;
		dgemm_r_diag_libstr(nu[N]+nx[N], ng[N], 1.0, &hsDCt[N], 0, 0, &hsQx[N], nb[N], 0.0, &hswork_mat_0, 0, 0, &hswork_mat_0, 0, 0);
		drowin_libstr(ng[N], 1.0, &hsqx[N], nb[N], &hswork_mat_0, nu[N]+nx[N], 0);
		dsyrk_dpotrf_ln_libstr(nu[N]+nx[N]+1, nu[N]+nx[N], ng[N], &hswork_mat_0, 0, 0, &hsDCt[N], 0, 0, &hsL[N], 0, 0, &hsL[N], 0, 0);
		}
	else
		{
		dpotrf_l_libstr(nu[N]+nx[N]+1, nu[N]+nx[N], &hsL[N], 0, 0, &hsL[N], 0, 0);
		}
	dtrtr_l_libstr(nx[N], 1.0, &hsL[N], nu[N], nu[N], &hsLxt[N], 0, 0);

	// middle stages
	for(nn=0; nn<N; nn++)
		{
		c_ptr = (char *) work;
		d_create_strmat(nu[N-nn-1]+nx[N-nn-1]+1, nx[N-nn]+ng[N-nn-1], &hswork_mat_0, (void *) c_ptr);
		c_ptr += hswork_mat_0.memory_size;
		if(ng[N-nn-1]>0)
			{
			d_create_strmat(nu[N-nn-1]+nx[N-nn-1], nx[N-nn]+ng[N-nn-1], &hswork_mat_1, (void *) c_ptr);
			c_ptr += hswork_mat_1.memory_size;
			}
		d_create_strvec(nx[N-nn], &hswork_vec_0, (void *) c_ptr);
		c_ptr += hswork_vec_0.memory_size;
		if(update_b)
			{
			drowin_libstr(nx[N-nn], 1.0, &hsb[N-nn-1], 0, &hsBAbt[N-nn-1], nu[N-nn-1]+nx[N-nn-1], 0);
			}
		dtrmm_rutn_libstr(nu[N-nn-1]+nx[N-nn-1]+1, nx[N-nn], 1.0, &hsBAbt[N-nn-1], 0, 0, &hsLxt[N-nn], 0, 0, 0.0, &hswork_mat_0, 0, 0, &hswork_mat_0, 0, 0);
		if(compute_Pb)
			{
			drowex_libstr(nx[N-nn], 1.0, &hswork_mat_0, nu[N-nn-1]+nx[N-nn-1], 0, &hswork_vec_0, 0);
			dtrmv_utn_libstr(nx[N-nn], &hsLxt[N-nn], 0, 0, &hswork_vec_0, 0, &hsPb[N-nn], 0);
			}
		dgead_libstr(1, nx[N-nn], 1.0, &hsL[N-nn], nu[N-nn]+nx[N-nn], nu[N-nn], &hswork_mat_0, nu[N-nn-1]+nx[N-nn-1], 0);
		dtrcp_l_libstr(nu[N-nn-1]+nx[N-nn-1], 1.0, &hsRSQrq[N-nn-1], 0, 0, &hsL[N-nn-1], 0, 0);
		if(update_q)
			{
			drowin_libstr(nu[N-nn-1]+nx[N-nn-1], 1.0, &hsrq[N-nn-1], 0, &hsL[N-nn-1], nu[N-nn-1]+nx[N-nn-1], 0);
			}
		else
			{
			dgecp_libstr(1, nu[N-nn-1]+nx[N-nn-1], 1.0, &hsRSQrq[N-nn-1], nu[N-nn-1]+nx[N-nn-1], 0, &hsL[N-nn-1], nu[N-nn-1]+nx[N-nn-1], 0);
			}
		if(nb[N-nn-1]>0)
			{
			ddiaad_libspstr(nb[N-nn-1], hidxb[N-nn-1], 1.0, &hsQx[N-nn-1], 0, &hsL[N-nn-1], 0, 0);
			drowad_libspstr(nb[N-nn-1], hidxb[N-nn-1], 1.0, &hsqx[N-nn-1], 0, &hsL[N-nn-1], nu[N-nn-1]+nx[N-nn-1], 0);
			}
		if(ng[N-nn-1]>0)
			{
			dgemm_r_diag_libstr(nu[N-nn-1]+nx[N-nn-1], ng[N-nn-1], 1.0, &hsDCt[N-nn-1], 0, 0, &hsQx[N-nn-1], nb[N-nn-1], 0.0, &hswork_mat_0, 0, nx[N-nn], &hswork_mat_0, 0, nx[N-nn]);
			drowin_libstr(ng[N-nn-1], 1.0, &hsqx[N-nn-1], nb[N-nn-1], &hswork_mat_0, nu[N-nn-1]+nx[N-nn-1], nx[N-nn]);
			dgecp_libstr(nu[N-nn-1]+nx[N-nn-1], nx[N-nn], 1.0, &hswork_mat_0, 0, 0, &hswork_mat_1, 0, 0);
			dgecp_libstr(nu[N-nn-1]+nx[N-nn-1], ng[N-nn-1], 1.0, &hsDCt[N-nn-1], 0, 0, &hswork_mat_1, 0, nx[N-nn]);
			dsyrk_dpotrf_ln_libstr(nu[N-nn-1]+nx[N-nn-1]+1, nu[N-nn-1]+nx[N-nn-1], nx[N-nn]+ng[N-nn-1], &hswork_mat_0, 0, 0, &hswork_mat_1, 0, 0, &hsL[N-nn-1], 0, 0, &hsL[N-nn-1], 0, 0);
			}
		else
			{
			dsyrk_dpotrf_ln_libstr(nu[N-nn-1]+nx[N-nn-1]+1, nu[N-nn-1]+nx[N-nn-1], nx[N-nn], &hswork_mat_0, 0, 0, &hswork_mat_0, 0, 0, &hsL[N-nn-1], 0, 0, &hsL[N-nn-1], 0, 0);
			}
		dtrtr_l_libstr(nx[N-nn-1], 1.0, &hsL[N-nn-1], nu[N-nn-1], nu[N-nn-1], &hsLxt[N-nn-1], 0, 0);
		}

	// forward substitution

	// first stage
	nn = 0;
	drowex_libstr(nu[nn]+nx[nn], -1.0, &hsL[nn], nu[nn]+nx[nn], 0, &hsux[nn], 0);
	dtrsv_ltn_libstr(nu[nn]+nx[nn], nu[nn]+nx[nn], &hsL[nn], 0, 0, &hsux[nn], 0, &hsux[nn], 0);
	drowex_libstr(nx[nn+1], 1.0, &hsBAbt[nn], nu[nn]+nx[nn], 0, &hsux[nn+1], nu[nn+1]);
	dgemv_t_libstr(nu[nn]+nx[nn], nx[nn+1], 1.0, &hsBAbt[nn], 0, 0, &hsux[nn], 0, 1.0, &hsux[nn+1], nu[nn+1], &hsux[nn+1], nu[nn+1]);
	if(compute_pi)
		{
		c_ptr = (char *) work;
		d_create_strvec(nx[nn+1], &hswork_vec_0, (void *) c_ptr);
		c_ptr += hswork_vec_0.memory_size;
		dveccp_libstr(nx[nn+1], 1.0, &hsux[nn+1], nu[nn+1], &hspi[nn+1], 0);
		drowex_libstr(nx[nn+1], 1.0, &hsL[nn+1], nu[nn+1]+nx[nn+1], nu[nn+1], &hswork_vec_0, 0);
		dtrmv_unn_libstr(nx[nn+1], &hsLxt[nn+1], 0, 0, &hspi[nn+1], 0, &hspi[nn+1], 0);
		daxpy_libstr(nx[nn+1], 1.0, &hswork_vec_0, 0, &hspi[nn+1], 0);
		dtrmv_utn_libstr(nx[nn+1], &hsLxt[nn+1], 0, 0, &hspi[nn+1], 0, &hspi[nn+1], 0);
		}

	// middle stages
	for(nn=1; nn<N; nn++)
		{
		drowex_libstr(nu[nn], -1.0, &hsL[nn], nu[nn]+nx[nn], 0, &hsux[nn], 0);
		dtrsv_ltn_libstr(nu[nn]+nx[nn], nu[nn], &hsL[nn], 0, 0, &hsux[nn], 0, &hsux[nn], 0);
		drowex_libstr(nx[nn+1], 1.0, &hsBAbt[nn], nu[nn]+nx[nn], 0, &hsux[nn+1], nu[nn+1]);
		dgemv_t_libstr(nu[nn]+nx[nn], nx[nn+1], 1.0, &hsBAbt[nn], 0, 0, &hsux[nn], 0, 1.0, &hsux[nn+1], nu[nn+1], &hsux[nn+1], nu[nn+1]);
		if(compute_pi)
			{
			c_ptr = (char *) work;
			d_create_strvec(nx[nn+1], &hswork_vec_0, (void *) c_ptr);
			c_ptr += hswork_vec_0.memory_size;
			dveccp_libstr(nx[nn+1], 1.0, &hsux[nn+1], nu[nn+1], &hspi[nn+1], 0);
			drowex_libstr(nx[nn+1], 1.0, &hsL[nn+1], nu[nn+1]+nx[nn+1], nu[nn+1], &hswork_vec_0, 0);
			dtrmv_unn_libstr(nx[nn+1], &hsLxt[nn+1], 0, 0, &hspi[nn+1], 0, &hspi[nn+1], 0);
			daxpy_libstr(nx[nn+1], 1.0, &hswork_vec_0, 0, &hspi[nn+1], 0);
			dtrmv_utn_libstr(nx[nn+1], &hsLxt[nn+1], 0, 0, &hspi[nn+1], 0, &hspi[nn+1], 0);
			}
		}

	return;

	}



void d_back_ric_rec_trf_libstr(int N, int *nx, int *nu, int *nb, int **hidxb, int *ng, struct d_strmat *hsBAbt, struct d_strmat *hsRSQrq, struct d_strmat *hsDCt, struct d_strvec *hsQx, struct d_strmat *hsL, struct d_strmat *hsLxt, void *work)
	{

	char *c_ptr;

	struct d_strmat hswork_mat_0, hswork_mat_1;

	int nn;

	// factorization

	// last stage
	dtrcp_l_libstr(nu[N]+nx[N], 1.0, &hsRSQrq[N], 0, 0, &hsL[N], 0, 0);
	if(nb[N]>0)
		{
		ddiaad_libspstr(nb[N], hidxb[N], 1.0, &hsQx[N], 0, &hsL[N], 0, 0);
		}
	if(ng[N]>0)
		{
		c_ptr = (char *) work;
		d_create_strmat(nu[N]+nx[N], ng[N], &hswork_mat_0, (void *) c_ptr);
		c_ptr += hswork_mat_0.memory_size;
		dgemm_r_diag_libstr(nu[N]+nx[N], ng[N], 1.0, &hsDCt[N], 0, 0, &hsQx[N], nb[N], 0.0, &hswork_mat_0, 0, 0, &hswork_mat_0, 0, 0);
		dsyrk_dpotrf_ln_libstr(nu[N]+nx[N], nu[N]+nx[N], ng[N], &hswork_mat_0, 0, 0, &hsDCt[N], 0, 0, &hsL[N], 0, 0, &hsL[N], 0, 0);
		}
	else
		{
		dpotrf_l_libstr(nu[N]+nx[N], nu[N]+nx[N], &hsL[N], 0, 0, &hsL[N], 0, 0);
		}
	dtrtr_l_libstr(nx[N], 1.0, &hsL[N], nu[N], nu[N], &hsLxt[N], 0, 0);

	// middle stages
	for(nn=0; nn<N; nn++)
		{
		c_ptr = (char *) work;
		d_create_strmat(nu[N-nn-1]+nx[N-nn-1], nx[N-nn]+ng[N-nn-1], &hswork_mat_0, (void *) c_ptr);
		c_ptr += hswork_mat_0.memory_size;
		if(ng[N-nn-1]>0)
			{
			d_create_strmat(nu[N-nn-1]+nx[N-nn-1], nx[N-nn]+ng[N-nn-1], &hswork_mat_1, (void *) c_ptr);
			c_ptr += hswork_mat_1.memory_size;
			}
		dtrmm_rutn_libstr(nu[N-nn-1]+nx[N-nn-1], nx[N-nn], 1.0, &hsBAbt[N-nn-1], 0, 0, &hsLxt[N-nn], 0, 0, 0.0, &hswork_mat_0, 0, 0, &hswork_mat_0, 0, 0);
		dtrcp_l_libstr(nu[N-nn-1]+nx[N-nn-1], 1.0, &hsRSQrq[N-nn-1], 0, 0, &hsL[N-nn-1], 0, 0);
		if(nb[N-nn-1]>0)
			{
			ddiaad_libspstr(nb[N-nn-1], hidxb[N-nn-1], 1.0, &hsQx[N-nn-1], 0, &hsL[N-nn-1], 0, 0);
			}
		if(ng[N-nn-1]>0)
			{
			dgemm_r_diag_libstr(nu[N-nn-1]+nx[N-nn-1], ng[N-nn-1], 1.0, &hsDCt[N-nn-1], 0, 0, &hsQx[N-nn-1], nb[N], 0.0, &hswork_mat_0, 0, nx[N-nn], &hswork_mat_0, 0, nx[N-nn]);
			dgecp_libstr(nu[N-nn-1]+nx[N-nn-1], nx[N-nn], 1.0, &hswork_mat_0, 0, 0, &hswork_mat_1, 0, 0);
			dgecp_libstr(nu[N-nn-1]+nx[N-nn-1], ng[N-nn-1], 1.0, &hsDCt[N-nn-1], 0, 0, &hswork_mat_1, 0, nx[N-nn]);
			dsyrk_dpotrf_ln_libstr(nu[N-nn-1]+nx[N-nn-1], nu[N-nn-1]+nx[N-nn-1], nx[N-nn]+ng[N-nn-1], &hswork_mat_0, 0, 0, &hswork_mat_1, 0, 0, &hsL[N-nn-1], 0, 0, &hsL[N-nn-1], 0, 0);
			}
		else
			{
			dsyrk_dpotrf_ln_libstr(nu[N-nn-1]+nx[N-nn-1], nu[N-nn-1]+nx[N-nn-1], nx[N-nn], &hswork_mat_0, 0, 0, &hswork_mat_0, 0, 0, &hsL[N-nn-1], 0, 0, &hsL[N-nn-1], 0, 0);
			}
		dtrtr_l_libstr(nx[N-nn-1], 1.0, &hsL[N-nn-1], nu[N-nn-1], nu[N-nn-1], &hsLxt[N-nn-1], 0, 0);
		}

	return;

	}



void d_back_ric_rec_trs_libstr(int N, int *nx, int *nu, int *nb, int **idxb, int *ng, struct d_strmat *hsBAbt, struct d_strvec *hsb, struct d_strvec *hsrq, struct d_strmat *hsDCt, struct d_strvec *hsqx, struct d_strvec *hsux, int compute_pi, struct d_strvec *hspi, int compute_Pb, struct d_strvec *hsPb, struct d_strmat *hsL, struct d_strmat *hsLxt, void *work)
	{

	char *c_ptr;

	struct d_strvec hswork_vec_0;

	int nn;

	// backward substitution

	// last stage
	dveccp_libstr(nu[N]+nx[N], 1.0, &hsrq[N], 0, &hsux[N], 0);
	if(nb[N]>0)
		{
		dvecad_libspstr(nb[N], idxb[N], 1.0, &hsqx[N], 0, &hsux[N], 0);
		}
	// general constraints
	if(ng[N]>0)
		{
		dgemv_n_libstr(nx[N], ng[N], 1.0, &hsDCt[N], 0, 0, &hsqx[N], nb[N], 1.0, &hsux[N], 0, &hsux[N], 0);
		}
//d_print_tran_strvec(nu[N]+nx[N], &hsux[N], 0);

	// middle stages
	for(nn=0; nn<N-1; nn++)
		{
		c_ptr = (char *) work;
		d_create_strvec(nx[N-nn], &hswork_vec_0, (void *) c_ptr);
		c_ptr += hswork_vec_0.memory_size;
		if(compute_Pb)
			{
			dtrmv_unn_libstr(nx[N-nn], &hsLxt[N-nn], 0, 0, &hsb[N-nn-1], 0, &hsPb[N-nn], 0);
			dtrmv_utn_libstr(nx[N-nn], &hsLxt[N-nn], 0, 0, &hsPb[N-nn], 0, &hsPb[N-nn], 0);
			}
		dveccp_libstr(nu[N-nn-1]+nx[N-nn-1], 1.0, &hsrq[N-nn-1], 0, &hsux[N-nn-1], 0);
		if(nb[N-nn-1]>0)
			{
			dvecad_libspstr(nb[N-nn-1], idxb[N-nn-1], 1.0, &hsqx[N-nn-1], 0, &hsux[N-nn-1], 0);
			}
		if(ng[N-nn-1]>0)
			{
			dgemv_n_libstr(nu[N-nn-1]+nx[N-nn-1], ng[N-nn-1], 1.0, &hsDCt[N-nn-1], 0, 0, &hsqx[N-nn-1], nb[N-nn-1], 1.0, &hsux[N-nn-1], 0, &hsux[N-nn-1], 0);
			}
		dveccp_libstr(nx[N-nn], 1.0, &hsPb[N-nn], 0, &hswork_vec_0, 0);
		daxpy_libstr(nx[N-nn], 1.0, &hsux[N-nn], nu[N-nn], &hswork_vec_0, 0);
		dgemv_n_libstr(nu[N-nn-1]+nx[N-nn-1], nx[N-nn], 1.0, &hsBAbt[N-nn-1], 0, 0, &hswork_vec_0, 0, 1.0, &hsux[N-nn-1], 0, &hsux[N-nn-1], 0);
		dtrsv_lnn_libstr(nu[N-nn-1]+nx[N-nn-1], nu[N-nn-1], &hsL[N-nn-1], 0, 0, &hsux[N-nn-1], 0, &hsux[N-nn-1], 0);
//d_print_tran_strvec(nu[N-nn-1]+nx[N-nn-1], &hsux[N-nn-1], 0);
		}

	// first stage
	nn = N-1;
	c_ptr = (char *) work;
	d_create_strvec(nx[N-nn], &hswork_vec_0, (void *) c_ptr);
	c_ptr += hswork_vec_0.memory_size;
	if(compute_Pb)
		{
		dtrmv_unn_libstr(nx[N-nn], &hsLxt[N-nn], 0, 0, &hsb[N-nn-1], 0, &hsPb[N-nn], 0);
		dtrmv_utn_libstr(nx[N-nn], &hsLxt[N-nn], 0, 0, &hsPb[N-nn], 0, &hsPb[N-nn], 0);
		}
	dveccp_libstr(nu[N-nn-1]+nx[N-nn-1], 1.0, &hsrq[N-nn-1], 0, &hsux[N-nn-1], 0);
	if(nb[N-nn-1]>0)
		{
		dvecad_libspstr(nb[N-nn-1], idxb[N-nn-1], 1.0, &hsqx[N-nn-1], 0, &hsux[N-nn-1], 0);
		}
	if(ng[N-nn-1]>0)
		{
		dgemv_n_libstr(nu[N-nn-1]+nx[N-nn-1], ng[N-nn-1], 1.0, &hsDCt[N-nn-1], 0, 0, &hsqx[N-nn-1], nb[N-nn-1], 1.0, &hsux[N-nn-1], 0, &hsux[N-nn-1], 0);
		}
	dveccp_libstr(nx[N-nn], 1.0, &hsPb[N-nn], 0, &hswork_vec_0, 0);
	daxpy_libstr(nx[N-nn], 1.0, &hsux[N-nn], nu[N-nn], &hswork_vec_0, 0);
	dgemv_n_libstr(nu[N-nn-1]+nx[N-nn-1], nx[N-nn], 1.0, &hsBAbt[N-nn-1], 0, 0, &hswork_vec_0, 0, 1.0, &hsux[N-nn-1], 0, &hsux[N-nn-1], 0);
	dtrsv_lnn_libstr(nu[N-nn-1]+nx[N-nn-1], nu[N-nn-1]+nx[N-nn-1], &hsL[N-nn-1], 0, 0, &hsux[N-nn-1], 0, &hsux[N-nn-1], 0);

//d_print_tran_strvec(nu[N-nn-1]+nx[N-nn-1], &hsux[N-nn-1], 0);
	// forward substitution

	// first stage
	nn = 0;
	if(compute_pi)
		{
		dveccp_libstr(nx[nn+1], 1.0, &hsux[nn+1], nu[nn+1], &hspi[nn+1], 0);
		}
	dveccp_libstr(nu[nn]+nx[nn], -1.0, &hsux[nn], 0, &hsux[nn], 0);
	dtrsv_ltn_libstr(nu[nn]+nx[nn], nu[nn]+nx[nn], &hsL[nn], 0, 0, &hsux[nn], 0, &hsux[nn], 0);
	dgemv_t_libstr(nu[nn]+nx[nn], nx[nn+1], 1.0, &hsBAbt[nn], 0, 0, &hsux[nn], 0, 1.0, &hsb[nn], 0, &hsux[nn+1], nu[nn+1]);
	if(compute_pi)
		{
		c_ptr = (char *) work;
		d_create_strvec(nx[nn+1], &hswork_vec_0, (void *) c_ptr);
		c_ptr += hswork_vec_0.memory_size;
		dveccp_libstr(nx[nn+1], 1.0, &hsux[nn+1], nu[nn+1], &hswork_vec_0, 0);
		dtrmv_unn_libstr(nx[nn+1], &hsLxt[nn+1], 0, 0, &hswork_vec_0, 0, &hswork_vec_0, 0);
		dtrmv_utn_libstr(nx[nn+1], &hsLxt[nn+1], 0, 0, &hswork_vec_0, 0, &hswork_vec_0, 0);
		daxpy_libstr(nx[nn+1], 1.0, &hswork_vec_0, 0, &hspi[nn+1], 0);
		}

	// middle stages
	for(nn=1; nn<N; nn++)
		{
		if(compute_pi)
			{
			dveccp_libstr(nx[nn+1], 1.0, &hsux[nn+1], nu[nn+1], &hspi[nn+1], 0);
			}
		dveccp_libstr(nu[nn], -1.0, &hsux[nn], 0, &hsux[nn], 0);
		dtrsv_ltn_libstr(nu[nn]+nx[nn], nu[nn], &hsL[nn], 0, 0, &hsux[nn], 0, &hsux[nn], 0);
		dgemv_t_libstr(nu[nn]+nx[nn], nx[nn+1], 1.0, &hsBAbt[nn], 0, 0, &hsux[nn], 0, 1.0, &hsb[nn], 0, &hsux[nn+1], nu[nn+1]);
		if(compute_pi)
			{
			c_ptr = (char *) work;
			d_create_strvec(nx[nn+1], &hswork_vec_0, (void *) c_ptr);
			c_ptr += hswork_vec_0.memory_size;
			dveccp_libstr(nx[nn+1], 1.0, &hsux[nn+1], nu[nn+1], &hswork_vec_0, 0);
			dtrmv_unn_libstr(nx[nn+1], &hsLxt[nn+1], 0, 0, &hswork_vec_0, 0, &hswork_vec_0, 0);
			dtrmv_utn_libstr(nx[nn+1], &hsLxt[nn+1], 0, 0, &hswork_vec_0, 0, &hswork_vec_0, 0);
			daxpy_libstr(nx[nn+1], 1.0, &hswork_vec_0, 0, &hspi[nn+1], 0);
			}
		}

	return;

	}



void d_back_ric_rec_sv_back_libstr(int N, int *nx, int *nu, int *nb, int **hidxb, int *ng, int update_b, struct d_strmat *hsBAbt, struct d_strvec *hsb, int update_q, struct d_strmat *hsRSQrq, struct d_strvec *hsrq, struct d_strmat *hsDCt, struct d_strvec *hsQx, struct d_strvec *hsqx, struct d_strvec *hsux, int compute_pi, struct d_strvec *hspi, int compute_Pb, struct d_strvec *hsPb, struct d_strmat *hsL, struct d_strmat *hsLxt, void *work)
	{

	char *c_ptr;

	struct d_strmat hswork_mat_0, hswork_mat_1;
	struct d_strvec hswork_vec_0;

	int nn;

	// factorization and backward substitution

	// last stage
	dtrcp_l_libstr(nu[N]+nx[N], 1.0, &hsRSQrq[N], 0, 0, &hsL[N], 0, 0); // TODO dtrcp_l_libstr with m and n, for m>=n
	if(update_q)
		{
		drowin_libstr(nu[N]+nx[N], 1.0, &hsrq[N], 0, &hsL[N], nu[N]+nx[N], 0);
		}
	else
		{
		dgecp_libstr(1, nu[N]+nx[N], 1.0, &hsRSQrq[N], nu[N]+nx[N], 0, &hsL[N], nu[N]+nx[N], 0);
		}
	if(nb[N]>0)
		{
		ddiaad_libspstr(nb[N], hidxb[N], 1.0, &hsQx[N], 0, &hsL[N], 0, 0);
		drowad_libspstr(nb[N], hidxb[N], 1.0, &hsqx[N], 0, &hsL[N], nu[N]+nx[N], 0);
		}
	if(ng[N]>0)
		{
		c_ptr = (char *) work;
		d_create_strmat(nu[N]+nx[N]+1, ng[N], &hswork_mat_0, (void *) c_ptr);
		c_ptr += hswork_mat_0.memory_size;
		dgemm_r_diag_libstr(nu[N]+nx[N], ng[N], 1.0, &hsDCt[N], 0, 0, &hsQx[N], nb[N], 0.0, &hswork_mat_0, 0, 0, &hswork_mat_0, 0, 0);
		drowin_libstr(ng[N], 1.0, &hsqx[N], nb[N], &hswork_mat_0, nu[N]+nx[N], 0);
		dsyrk_dpotrf_ln_libstr(nu[N]+nx[N]+1, nu[N]+nx[N], ng[N], &hswork_mat_0, 0, 0, &hsDCt[N], 0, 0, &hsL[N], 0, 0, &hsL[N], 0, 0);
		}
	else
		{
		dpotrf_l_libstr(nu[N]+nx[N]+1, nu[N]+nx[N], &hsL[N], 0, 0, &hsL[N], 0, 0);
		}
	dtrtr_l_libstr(nx[N], 1.0, &hsL[N], nu[N], nu[N], &hsLxt[N], 0, 0);

	// middle stages
	for(nn=0; nn<N; nn++)
		{
		c_ptr = (char *) work;
		d_create_strmat(nu[N-nn-1]+nx[N-nn-1]+1, nx[N-nn]+ng[N-nn-1], &hswork_mat_0, (void *) c_ptr);
		c_ptr += hswork_mat_0.memory_size;
		if(ng[N-nn-1]>0)
			{
			d_create_strmat(nu[N-nn-1]+nx[N-nn-1], nx[N-nn]+ng[N-nn-1], &hswork_mat_1, (void *) c_ptr);
			c_ptr += hswork_mat_1.memory_size;
			}
		d_create_strvec(nx[N-nn], &hswork_vec_0, (void *) c_ptr);
		c_ptr += hswork_vec_0.memory_size;
		if(update_b)
			{
			drowin_libstr(nx[N-nn], 1.0, &hsb[N-nn-1], 0, &hsBAbt[N-nn-1], nu[N-nn-1]+nu[N-nn-1], 0);
			}
		dtrmm_rutn_libstr(nu[N-nn-1]+nx[N-nn-1]+1, nx[N-nn], 1.0, &hsBAbt[N-nn-1], 0, 0, &hsLxt[N-nn], 0, 0, 0.0, &hswork_mat_0, 0, 0, &hswork_mat_0, 0, 0);
		if(compute_Pb)
			{
			drowex_libstr(nx[N-nn], 1.0, &hswork_mat_0, nu[N-nn-1]+nx[N-nn-1], 0, &hswork_vec_0, 0);
			dtrmv_utn_libstr(nx[N-nn], &hsLxt[N-nn], 0, 0, &hswork_vec_0, 0, &hsPb[N-nn], 0);
			}
		dgead_libstr(1, nx[N-nn], 1.0, &hsL[N-nn], nu[N-nn]+nx[N-nn], nu[N-nn], &hswork_mat_0, nu[N-nn-1]+nx[N-nn-1], 0);
		dtrcp_l_libstr(nu[N-nn-1]+nx[N-nn-1], 1.0, &hsRSQrq[N-nn-1], 0, 0, &hsL[N-nn-1], 0, 0);
		if(update_q)
			{
			drowin_libstr(nu[N-nn-1]+nx[N-nn-1], 1.0, &hsrq[N-nn-1], 0, &hsL[N-nn-1], nu[N-nn-1]+nx[N-nn-1], 0);
			}
		else
			{
			dgecp_libstr(1, nu[N-nn-1]+nx[N-nn-1], 1.0, &hsRSQrq[N-nn-1], nu[N-nn-1]+nx[N-nn-1], 0, &hsL[N-nn-1], nu[N-nn-1]+nx[N-nn-1], 0);
			}
		if(nb[N-nn-1]>0)
			{
			ddiaad_libspstr(nb[N-nn-1], hidxb[N-nn-1], 1.0, &hsQx[N-nn-1], 0, &hsL[N-nn-1], 0, 0);
			drowad_libspstr(nb[N-nn-1], hidxb[N-nn-1], 1.0, &hsqx[N-nn-1], 0, &hsL[N-nn-1], nu[N-nn-1]+nx[N-nn-1], 0);
			}
		if(ng[N-nn-1]>0)
			{
			dgemm_r_diag_libstr(nu[N-nn-1]+nx[N-nn-1], ng[N-nn-1], 1.0, &hsDCt[N-nn-1], 0, 0, &hsQx[N-nn-1], nb[N], 0.0, &hswork_mat_0, 0, nx[N-nn], &hswork_mat_0, 0, nx[N-nn]);
			drowin_libstr(ng[N-nn-1], 1.0, &hsqx[N-nn-1], nb[N-nn-1], &hswork_mat_0, nu[N-nn-1]+nx[N-nn-1], nx[N-nn]);
			dgecp_libstr(nu[N-nn-1]+nx[N-nn-1], nx[N-nn], 1.0, &hswork_mat_0, 0, 0, &hswork_mat_1, 0, 0);
			dgecp_libstr(nu[N-nn-1]+nx[N-nn-1], ng[N-nn-1], 1.0, &hsDCt[N-nn-1], 0, 0, &hswork_mat_1, 0, nx[N-nn]);
			dsyrk_dpotrf_ln_libstr(nu[N-nn-1]+nx[N-nn-1]+1, nu[N-nn-1]+nx[N-nn-1], nx[N-nn]+ng[N-nn-1], &hswork_mat_0, 0, 0, &hswork_mat_1, 0, 0, &hsL[N-nn-1], 0, 0, &hsL[N-nn-1], 0, 0);
			}
		else
			{
			dsyrk_dpotrf_ln_libstr(nu[N-nn-1]+nx[N-nn-1]+1, nu[N-nn-1]+nx[N-nn-1], nx[N-nn], &hswork_mat_0, 0, 0, &hswork_mat_0, 0, 0, &hsL[N-nn-1], 0, 0, &hsL[N-nn-1], 0, 0);
			}
		dtrtr_l_libstr(nx[N-nn-1], 1.0, &hsL[N-nn-1], nu[N-nn-1], nu[N-nn-1], &hsLxt[N-nn-1], 0, 0);
		}


	return;

	}



void d_back_ric_rec_sv_forw_libstr(int N, int *nx, int *nu, int *nb, int **hidxb, int *ng, int update_b, struct d_strmat *hsBAbt, struct d_strvec *hsb, int update_q, struct d_strmat *hsRSQrq, struct d_strvec *hsrq, struct d_strmat *hsDCt, struct d_strvec *hsQx, struct d_strvec *hsqx, struct d_strvec *hsux, int compute_pi, struct d_strvec *hspi, int compute_Pb, struct d_strvec *hsPb, struct d_strmat *hsL, struct d_strmat *hsLxt, void *work)
	{

	char *c_ptr;

	struct d_strvec hswork_vec_0;

	int nn;

	// forward substitution (without first stage)

	// middle stages
	for(nn=0; nn<N; nn++)
		{
		drowex_libstr(nu[nn], -1.0, &hsL[nn], nu[nn]+nx[nn], 0, &hsux[nn], 0);
		dtrsv_ltn_libstr(nu[nn]+nx[nn], nu[nn], &hsL[nn], 0, 0, &hsux[nn], 0, &hsux[nn], 0);
		drowex_libstr(nx[nn+1], 1.0, &hsBAbt[nn], nu[nn]+nx[nn], 0, &hsux[nn+1], nu[nn+1]);
		dgemv_t_libstr(nu[nn]+nx[nn], nx[nn+1], 1.0, &hsBAbt[nn], 0, 0, &hsux[nn], 0, 1.0, &hsux[nn+1], nu[nn+1], &hsux[nn+1], nu[nn+1]);
		if(compute_pi)
			{
			c_ptr = (char *) work;
			d_create_strvec(nx[nn+1], &hswork_vec_0, (void *) c_ptr);
			c_ptr += hswork_vec_0.memory_size;
			dveccp_libstr(nx[nn+1], 1.0, &hsux[nn+1], nu[nn+1], &hspi[nn+1], 0);
			drowex_libstr(nx[nn+1], 1.0, &hsL[nn+1], nu[nn+1]+nx[nn+1], nu[nn+1], &hswork_vec_0, 0);
			dtrmv_unn_libstr(nx[nn+1], &hsLxt[nn+1], 0, 0, &hspi[nn+1], 0, &hspi[nn+1], 0);
			daxpy_libstr(nx[nn+1], 1.0, &hswork_vec_0, 0, &hspi[nn+1], 0);
			dtrmv_utn_libstr(nx[nn+1], &hsLxt[nn+1], 0, 0, &hspi[nn+1], 0, &hspi[nn+1], 0);
			}
		}

	return;

	}



#endif
