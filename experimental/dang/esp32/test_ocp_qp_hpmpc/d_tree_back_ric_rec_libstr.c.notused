/**************************************************************************************************
* acados/experimental/dang/esp32/test_ocp_qp_hpmpc/                                                                                                *
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

#include "tree.h"



// work space
int d_tree_back_ric_rec_work_space_size_bytes_libstr(int Nn, struct node *tree, int *nx, int *nu, int *nb, int *ng)
	{

	int ii;

	int size = 0;

	// max sizes
	int nxM  = 0;
	int ngM = 0;
	int nuxM  = 0;
	for(ii=0; ii<Nn; ii++)
		{
		nxM = nx[ii]>nxM ? nx[ii] : nxM;
		ngM = ng[ii]>ngM ? ng[ii] : ngM;
		nuxM = nu[ii]+nx[ii]>nuxM ? nu[ii]+nx[ii]+1 : nuxM;
		}
	int nxgM = nxM>ngM ? nxM : ngM;

	size += d_size_strmat(nuxM+1, nxgM); // ric_work_mat[0]
	size += d_size_strvec(nxM); // ric_work_vec[0]

	return size;
	}



// help routines

void d_back_ric_trf_funnel1_libstr(int nkids, int nx0, int nx1, int nu0, int nb0, int *hidxb0, int ng0, struct d_strmat *hsBAbt, struct d_strmat *hsRSQrq, struct d_strmat *hsDCt, struct d_strvec *hsQx, struct d_strmat *hsL, struct d_strmat *hsLxt0, struct d_strmat *hsLxt1, void *work)
	{

	char *c_ptr;

	struct d_strmat hswork_mat_0, hswork_mat_1;

	int ii;

	c_ptr = (char *) work;
	d_create_strmat(nu0+nx0, nx1, &hswork_mat_0, (void *) c_ptr);
	c_ptr += hswork_mat_0.memory_size;

	// initialize with hessian
	dtrcp_l_libstr(nu0+nx0, 1.0, &hsRSQrq[0], 0, 0, &hsL[0], 0, 0);

	// all kids: update
	for(ii=0; ii<nkids; ii++)
		{
		dtrmm_rutn_libstr(nu0+nx0, nx1, 1.0, &hsBAbt[ii], 0, 0, &hsLxt1[ii], 0, 0, 0.0, &hswork_mat_0, 0, 0, &hswork_mat_0, 0, 0);
		dsyrk_ln_libstr(nu0+nx0, nu0+nx0, nx1, 1.0, &hswork_mat_0, 0, 0, &hswork_mat_0, 0, 0, 1.0, &hsL[0], 0, 0, &hsL[0], 0, 0);
		}

	// update with box constraints
	if(nb0>0)
		{
		ddiaad_libspstr(nb0, hidxb0, 1.0, &hsQx[0], 0, &hsL[0], 0, 0);
		}
	// update with general constraints and factorize at the end
	if(ng0>0)
		{
		c_ptr = (char *) work;
		d_create_strmat(nu0+nx0, ng0, &hswork_mat_0, (void *) c_ptr);
		c_ptr += hswork_mat_0.memory_size;
		dgemm_r_diag_libstr(nu0+nx0, ng0, 1.0, &hsDCt[0], 0, 0, &hsQx[0], nb0, 0.0, &hswork_mat_0, 0, 0, &hswork_mat_0, 0, 0);
		dsyrk_dpotrf_ln_libstr(nu0+nx0, nu0+nx0, ng0, &hsDCt[0], 0, 0, &hswork_mat_0, 0, 0, &hsL[0], 0, 0, &hsL[0], 0, 0);
//		printf("\nd_back_ric_trf_funnel1_libstr: ng0>0: feature not implemented yet\n\n");
//		exit(1);
		}
	else
		{
		dpotrf_l_libstr(nu0+nx0, nu0+nx0, &hsL[0], 0, 0, &hsL[0], 0, 0);
		}
	dtrtr_l_libstr(nx0, 1.0, &hsL[0], nu0, nu0, &hsLxt0[0], 0, 0);

	return;

	}



void d_back_ric_trf_legN_libstr(int nx0, int nb0, int *hidxb0, int ng0, struct d_strmat *hsRSQrq, struct d_strmat *hsDCt, struct d_strvec *hsQx, struct d_strmat *hsL, struct d_strmat *hsLxt, void *work)
	{

	char *c_ptr;

	struct d_strmat hswork_mat_0;

	dtrcp_l_libstr(nx0, 1.0, &hsRSQrq[0], 0, 0, &hsL[0], 0, 0);
	if(nb0>0)
		{
		ddiaad_libspstr(nb0, hidxb0, 1.0, &hsQx[0], 0, &hsL[0], 0, 0);
		}
	if(ng0>0)
		{
		c_ptr = (char *) work;
		d_create_strmat(nx0, ng0, &hswork_mat_0, (void *) c_ptr);
		c_ptr += hswork_mat_0.memory_size;
		dgemm_r_diag_libstr(nx0, ng0, 1.0, &hsDCt[0], 0, 0, &hsQx[0], nb0, 0.0, &hswork_mat_0, 0, 0, &hswork_mat_0, 0, 0);
		dsyrk_dpotrf_ln_libstr(nx0, nx0, ng0, &hswork_mat_0, 0, 0, &hsDCt[0], 0, 0, &hsL[0], 0, 0, &hsL[0], 0, 0);
//		printf("\nd_back_ric_trf_legN_libstr: ng>0: feature not implemented yet\n\n");
//		exit(1);
		}
	else
		{
		dpotrf_l_libstr(nx0, nx0, &hsL[0], 0, 0, &hsL[0], 0, 0);
		}
	dtrtr_l_libstr(nx0, 1.0, &hsL[0], 0, 0, &hsLxt[0], 0, 0);


	return;

	}



void d_back_ric_trs_back_leg0_libstr(int nx0, int nx1, int nu0, int nu1, int nb0, int *hidxb0, int ng0, struct d_strmat *hsBAbt, struct d_strvec *hsb, struct d_strvec *hsrq, struct d_strmat *hsDCt, struct d_strvec *hsqx, struct d_strmat *hsL, struct d_strmat *hsLxt, int compute_Pb, struct d_strvec *hsPb, struct d_strvec *hsux0, struct d_strvec *hsux1, void *work)
	{

	char *c_ptr;

	struct d_strvec hswork_vec_0;

	c_ptr = (char *) work;
	d_create_strvec(nx1, &hswork_vec_0, (void *) c_ptr);
	c_ptr += hswork_vec_0.memory_size;
	if(compute_Pb)
		{
		dtrmv_unn_libstr(nx1, &hsLxt[0], 0, 0, &hsb[0], 0, &hsPb[0], 0);
		dtrmv_utn_libstr(nx1, &hsLxt[0], 0, 0, &hsPb[0], 0, &hsPb[0], 0);
		}
	dveccp_libstr(nu0+nx0, 1.0, &hsrq[0], 0, &hsux0[0], 0);
	if(nb0>0)
		{
		dvecad_libspstr(nb0, hidxb0, 1.0, &hsqx[0], 0, &hsux0[0], 0);
		}
	if(ng0>0)
		{
		dgemv_n_libstr(nu0+nx0, ng0, 1.0, &hsDCt[0], 0, 0, &hsqx[0], nb0, 1.0, &hsux0[0], 0, &hsux0[0], 0);
		}
	dveccp_libstr(nx1, 1.0, &hsPb[0], 0, &hswork_vec_0, 0);
	daxpy_libstr(nx1, 1.0, &hsux1[0], nu1, &hswork_vec_0, 0);
	dgemv_n_libstr(nu0+nx0, nx1, 1.0, &hsBAbt[0], 0, 0, &hswork_vec_0, 0, 1.0, &hsux0[0], 0, &hsux0[0], 0);
	dtrsv_lnn_libstr(nu0+nx0, nu0+nx0, &hsL[0], 0, 0, &hsux0[0], 0, &hsux0[0], 0);

	return;

	}



void d_back_ric_trs_back_leg1_libstr(int nx0, int nx1, int nu0, int nu1, int nb0, int *hidxb0, int ng0, struct d_strmat *hsBAbt, struct d_strvec *hsb, struct d_strvec *hsrq, struct d_strmat *hsDCt, struct d_strvec *hsqx, struct d_strmat *hsL, struct d_strmat *hsLxt, int compute_Pb, struct d_strvec *hsPb, struct d_strvec *hsux0, struct d_strvec *hsux1, void *work)
	{

	char *c_ptr;

	struct d_strvec hswork_vec_0;

	c_ptr = (char *) work;
	d_create_strvec(nx1, &hswork_vec_0, (void *) c_ptr);
	c_ptr += hswork_vec_0.memory_size;
	if(compute_Pb)
		{
		dtrmv_unn_libstr(nx1, &hsLxt[0], 0, 0, &hsb[0], 0, &hsPb[0], 0);
		dtrmv_utn_libstr(nx1, &hsLxt[0], 0, 0, &hsPb[0], 0, &hsPb[0], 0);
		}
	dveccp_libstr(nu0+nx0, 1.0, &hsrq[0], 0, &hsux0[0], 0);
	if(nb0>0)
		{
		dvecad_libspstr(nb0, hidxb0, 1.0, &hsqx[0], 0, &hsux0[0], 0);
		}
	if(ng0>0)
		{
		dgemv_n_libstr(nu0+nx0, ng0, 1.0, &hsDCt[0], 0, 0, &hsqx[0], nb0, 1.0, &hsux0[0], 0, &hsux0[0], 0);
		}
	dveccp_libstr(nx1, 1.0, &hsPb[0], 0, &hswork_vec_0, 0);
	daxpy_libstr(nx1, 1.0, &hsux1[0], nu1, &hswork_vec_0, 0);
	dgemv_n_libstr(nu0+nx0, nx1, 1.0, &hsBAbt[0], 0, 0, &hswork_vec_0, 0, 1.0, &hsux0[0], 0, &hsux0[0], 0);
	dtrsv_lnn_libstr(nu0+nx0, nu0, &hsL[0], 0, 0, &hsux0[0], 0, &hsux0[0], 0);

	return;

	}



void d_back_ric_trs_back_legN_libstr(int nx0, int nu0, int nb0, int *hidxb0, int ng0, struct d_strvec *hsrq, struct d_strmat *hsDCt, struct d_strvec *hsqx, struct d_strvec *hsux)
	{

	dveccp_libstr(nu0+nx0, 1.0, &hsrq[0], 0, &hsux[0], 0);
	if(nb0>0)
		{
		dvecad_libspstr(nb0, hidxb0, 1.0, &hsqx[0], 0, &hsux[0], 0);
		}
	if(ng0>0)
		{
		dgemv_n_libstr(nx0, ng0, 1.0, &hsDCt[0], 0, 0, &hsqx[0], nb0, 1.0, &hsux[0], 0, &hsux[0], 0);
		}

	return;

	}



void d_back_ric_trs_back_funnel0_libstr(int nkids, int nx0, int nx1, int nu0, int nu1, int nb0, int *hidxb0, int ng0, struct d_strmat *hsBAbt, struct d_strvec *hsb, struct d_strvec *hsrq, struct d_strmat *hsDCt, struct d_strvec *hsqx, struct d_strmat *hsL, struct d_strmat *hsLxt, int compute_Pb, struct d_strvec *hsPb, struct d_strvec *hsux0, struct d_strvec *hsux1, void *work)
	{

	int ii;

	char *c_ptr;

	struct d_strvec hswork_vec_0;

	c_ptr = (char *) work;
	d_create_strvec(nx1, &hswork_vec_0, (void *) c_ptr);
	c_ptr += hswork_vec_0.memory_size;

	// initialize with gradient
	dveccp_libstr(nu0+nx0, 1.0, &hsrq[0], 0, &hsux0[0], 0);

	// all kids: update
	for(ii=0; ii<nkids; ii++)
		{
		if(compute_Pb)
			{
			dtrmv_unn_libstr(nx1, &hsLxt[ii], 0, 0, &hsb[ii], 0, &hsPb[ii], 0);
			dtrmv_utn_libstr(nx1, &hsLxt[ii], 0, 0, &hsPb[ii], 0, &hsPb[ii], 0);
			}
		dveccp_libstr(nx1, 1.0, &hsPb[ii], 0, &hswork_vec_0, 0);
		daxpy_libstr(nx1, 1.0, &hsux1[ii], nu1, &hswork_vec_0, 0);
		dgemv_n_libstr(nu0+nx0, nx1, 1.0, &hsBAbt[ii], 0, 0, &hswork_vec_0, 0, 1.0, &hsux0[0], 0, &hsux0[0], 0);
		}

	// update with box constraints
	if(nb0>0)
		{
		dvecad_libspstr(nb0, hidxb0, 1.0, &hsqx[0], 0, &hsux0[0], 0);
		}
	// update with general constraints
	if(ng0>0)
		{
		dgemv_n_libstr(nu0+nx0, ng0, 1.0, &hsDCt[0], 0, 0, &hsqx[0], nb0, 1.0, &hsux0[0], 0, &hsux0[0], 0);
		}
	// solve at the end
	dtrsv_lnn_libstr(nu0+nx0, nu0+nx0, &hsL[0], 0, 0, &hsux0[0], 0, &hsux0[0], 0);

	return;

	}



void d_back_ric_trs_back_funnel1_libstr(int nkids, int nx0, int nx1, int nu0, int nu1, int nb0, int *hidxb0, int ng0, struct d_strmat *hsBAbt, struct d_strvec *hsb, struct d_strvec *hsrq, struct d_strmat *hsDCt, struct d_strvec *hsqx, struct d_strmat *hsL, struct d_strmat *hsLxt, int compute_Pb, struct d_strvec *hsPb, struct d_strvec *hsux0, struct d_strvec *hsux1, void *work)
	{

	int ii;

	char *c_ptr;

	struct d_strvec hswork_vec_0;

	c_ptr = (char *) work;
	d_create_strvec(nx1, &hswork_vec_0, (void *) c_ptr);
	c_ptr += hswork_vec_0.memory_size;

	// initialize with gradient
	dveccp_libstr(nu0+nx0, 1.0, &hsrq[0], 0, &hsux0[0], 0);

	// all kids: update
	for(ii=0; ii<nkids; ii++)
		{
		if(compute_Pb)
			{
			dtrmv_unn_libstr(nx1, &hsLxt[ii], 0, 0, &hsb[ii], 0, &hsPb[ii], 0);
			dtrmv_utn_libstr(nx1, &hsLxt[ii], 0, 0, &hsPb[ii], 0, &hsPb[ii], 0);
			}
		dveccp_libstr(nx1, 1.0, &hsPb[ii], 0, &hswork_vec_0, 0);
		daxpy_libstr(nx1, 1.0, &hsux1[ii], nu1, &hswork_vec_0, 0);
		dgemv_n_libstr(nu0+nx0, nx1, 1.0, &hsBAbt[ii], 0, 0, &hswork_vec_0, 0, 1.0, &hsux0[0], 0, &hsux0[0], 0);
		}

	// update with box constraints
	if(nb0>0)
		{
		dvecad_libspstr(nb0, hidxb0, 1.0, &hsqx[0], 0, &hsux0[0], 0);
		}
	// update with general constraints
	if(ng0>0)
		{
		dgemv_n_libstr(nu0+nx0, ng0, 1.0, &hsDCt[0], 0, 0, &hsqx[0], nb0, 1.0, &hsux0[0], 0, &hsux0[0], 0);
		}
	// solve at the end
	dtrsv_lnn_libstr(nu0+nx0, nu0, &hsL[0], 0, 0, &hsux0[0], 0, &hsux0[0], 0);

	return;

	}



void d_back_ric_trs_forw_leg0_libstr(int nx0, int nx1, int nu0, int nu1, struct d_strmat *hsBAbt, struct d_strvec *hsb, struct d_strmat *hsL, struct d_strmat *hsLxt, struct d_strvec *hsux0, struct d_strvec *hsux1, int compute_pi, struct d_strvec *hspi, void *work)
	{

	char *c_ptr;

	struct d_strvec hswork_vec_0;

	if(compute_pi)
		{
		dveccp_libstr(nx1, 1.0, &hsux1[0], nu1, &hspi[0], 0);
		}
	dveccp_libstr(nu0+nx0, -1.0, &hsux0[0], 0, &hsux0[0], 0);
	dtrsv_ltn_libstr(nu0+nx0, nu0+nx0, &hsL[0], 0, 0, &hsux0[0], 0, &hsux0[0], 0);
	dgemv_t_libstr(nu0+nx0, nx1, 1.0, &hsBAbt[0], 0, 0, &hsux0[0], 0, 1.0, &hsb[0], 0, &hsux1[0], nu1);
	if(compute_pi)
		{
		c_ptr = (char *) work;
		d_create_strvec(nx1, &hswork_vec_0, (void *) c_ptr);
		c_ptr += hswork_vec_0.memory_size;
		dveccp_libstr(nx1, 1.0, &hsux1[0], nu1, &hswork_vec_0, 0);
		dtrmv_unn_libstr(nx1, &hsLxt[0], 0, 0, &hswork_vec_0, 0, &hswork_vec_0, 0);
		dtrmv_utn_libstr(nx1, &hsLxt[0], 0, 0, &hswork_vec_0, 0, &hswork_vec_0, 0);
		daxpy_libstr(nx1, 1.0, &hswork_vec_0, 0, &hspi[0], 0);
		}

	return;

	}



void d_back_ric_trs_forw_leg1_libstr(int nx0, int nx1, int nu0, int nu1, struct d_strmat *hsBAbt, struct d_strvec *hsb, struct d_strmat *hsL, struct d_strmat *hsLxt, struct d_strvec *hsux0, struct d_strvec *hsux1, int compute_pi, struct d_strvec *hspi, void *work)
	{

	char *c_ptr;

	struct d_strvec hswork_vec_0;

	if(compute_pi)
		{
		c_ptr = (char *) work;
		d_create_strvec(nx1, &hswork_vec_0, (void *) c_ptr);
		c_ptr += hswork_vec_0.memory_size;
		dveccp_libstr(nx1, 1.0, &hsux1[0], nu1, &hspi[0], 0);
		}
	dveccp_libstr(nu0, -1.0, &hsux0[0], 0, &hsux0[0], 0);
	dtrsv_ltn_libstr(nu0+nx0, nu0, &hsL[0], 0, 0, &hsux0[0], 0, &hsux0[0], 0);
	dgemv_t_libstr(nu0+nx0, nx1, 1.0, &hsBAbt[0], 0, 0, &hsux0[0], 0, 1.0, &hsb[0], 0, &hsux1[0], nu1);
	if(compute_pi)
		{
		dveccp_libstr(nx1, 1.0, &hsux1[0], nu1, &hswork_vec_0, 0);
		dtrmv_unn_libstr(nx1, &hsLxt[0], 0, 0, &hswork_vec_0, 0, &hswork_vec_0, 0);
		dtrmv_utn_libstr(nx1, &hsLxt[0], 0, 0, &hswork_vec_0, 0, &hswork_vec_0, 0);
		daxpy_libstr(nx1, 1.0, &hswork_vec_0, 0, &hspi[0], 0);
		}

	return;

	}



void d_back_ric_trs_forw_funnel0_libstr(int nkids, int nx0, int nx1, int nu0, int nu1, struct d_strmat *hsBAbt, struct d_strvec *hsb, struct d_strmat *hsL, struct d_strmat *hsLxt, struct d_strvec *hsux0, struct d_strvec *hsux1, int compute_pi, struct d_strvec *hspi, void *work)
	{

	char *c_ptr;

	struct d_strvec hswork_vec_0;

	int ii;

	if(compute_pi)
		{
		c_ptr = (char *) work;
		d_create_strvec(nx1, &hswork_vec_0, (void *) c_ptr);
		c_ptr += hswork_vec_0.memory_size;
		for(ii=0; ii<nkids; ii++)
			{
			dveccp_libstr(nx1, 1.0, &hsux1[ii], nu1, &hspi[ii], 0);
			}
		}
	dveccp_libstr(nu0+nx0, -1.0, &hsux0[0], 0, &hsux0[0], 0);
	dtrsv_ltn_libstr(nu0+nx0, nu0+nx0, &hsL[0], 0, 0, &hsux0[0], 0, &hsux0[0], 0);
	for(ii=0; ii<nkids; ii++)
		{
		dgemv_t_libstr(nu0+nx0, nx1, 1.0, &hsBAbt[ii], 0, 0, &hsux0[0], 0, 1.0, &hsb[ii], 0, &hsux1[ii], nu1);
		if(compute_pi)
			{
			dveccp_libstr(nx1, 1.0, &hsux1[ii], nu1, &hswork_vec_0, 0);
			dtrmv_unn_libstr(nx1, &hsLxt[ii], 0, 0, &hswork_vec_0, 0, &hswork_vec_0, 0);
			dtrmv_utn_libstr(nx1, &hsLxt[ii], 0, 0, &hswork_vec_0, 0, &hswork_vec_0, 0);
			daxpy_libstr(nx1, 1.0, &hswork_vec_0, 0, &hspi[ii], 0);
			}
		}

	return;

	}



void d_back_ric_trs_forw_funnel1_libstr(int nkids, int nx0, int nx1, int nu0, int nu1, struct d_strmat *hsBAbt, struct d_strvec *hsb, struct d_strmat *hsL, struct d_strmat *hsLxt, struct d_strvec *hsux0, struct d_strvec *hsux1, int compute_pi, struct d_strvec *hspi, void *work)
	{

	char *c_ptr;

	struct d_strvec hswork_vec_0;

	int ii;

	if(compute_pi)
		{
		c_ptr = (char *) work;
		d_create_strvec(nx1, &hswork_vec_0, (void *) c_ptr);
		c_ptr += hswork_vec_0.memory_size;
		for(ii=0; ii<nkids; ii++)
			{
			dveccp_libstr(nx1, 1.0, &hsux1[ii], nu1, &hspi[ii], 0);
			}
		}
	dveccp_libstr(nu0, -1.0, &hsux0[0], 0, &hsux0[0], 0);
	dtrsv_ltn_libstr(nu0+nx0, nu0, &hsL[0], 0, 0, &hsux0[0], 0, &hsux0[0], 0);
	for(ii=0; ii<nkids; ii++)
		{
		dgemv_t_libstr(nu0+nx0, nx1, 1.0, &hsBAbt[ii], 0, 0, &hsux0[0], 0, 1.0, &hsb[ii], 0, &hsux1[ii], nu1);
		if(compute_pi)
			{
			dveccp_libstr(nx1, 1.0, &hsux1[ii], nu1, &hswork_vec_0, 0);
			dtrmv_unn_libstr(nx1, &hsLxt[ii], 0, 0, &hswork_vec_0, 0, &hswork_vec_0, 0);
			dtrmv_utn_libstr(nx1, &hsLxt[ii], 0, 0, &hswork_vec_0, 0, &hswork_vec_0, 0);
			daxpy_libstr(nx1, 1.0, &hswork_vec_0, 0, &hspi[ii], 0);
			}
		}

	return;

	}



// Riccati recursion routines

void d_tree_back_ric_rec_trf_libstr(int Nn, struct node *tree, int *nx, int *nu, int *nb, int **hidxb, int *ng, struct d_strmat *hsBAbt, struct d_strmat *hsRSQrq, struct d_strmat *hsDCt, struct d_strvec *hsQx, struct d_strmat *hsL, struct d_strmat *hsLxt, void *work)
	{

	int nn;

	int dad, nkids, idxkid;

	// factorization

	// process one node at the time, starting from the last one
	for(nn=Nn-1; nn>=0; nn--)
		{
		dad = tree[nn].dad;
		nkids = tree[nn].nkids;
		if(nkids==0) // has no kids: last stage
			{
			d_back_ric_trf_legN_libstr(nx[nn], nb[nn], hidxb[nn], ng[nn], &hsRSQrq[nn], &hsDCt[nn], &hsQx[nn], &hsL[nn], &hsLxt[nn], work);
			}
		else // has one kid: middle stages
			{
			idxkid = tree[nn].kids[0];
			d_back_ric_trf_funnel1_libstr(nkids, nx[nn], nx[idxkid], nu[nn], nb[nn], hidxb[nn], ng[nn], &hsBAbt[idxkid], &hsRSQrq[nn], &hsDCt[nn], &hsQx[nn], &hsL[nn], &hsLxt[nn], &hsLxt[idxkid], work);
			}
		}

	return;

	}



void d_tree_back_ric_rec_trs_libstr(int Nn, struct node *tree, int *nx, int *nu, int *nb, int **hidxb, int *ng, struct d_strmat *hsBAbt, struct d_strvec *hsb, struct d_strvec *hsrq, struct d_strmat *hsDCt, struct d_strvec *hsqx, struct d_strvec *hsux, int compute_pi, struct d_strvec *hspi, int compute_Pb, struct d_strvec *hsPb, struct d_strmat *hsL, struct d_strmat *hsLxt, void *work)
	{

	int nn;

	int dad, nkids, idxkid;

	// backward substitution

	// process one node at the time, starting from the last one
	for(nn=Nn-1; nn>=0; nn--)
		{
		dad = tree[nn].dad;
		nkids = tree[nn].nkids;
		if(dad<0) // root
			{
			if(nkids>1) // has many kids => funnel
				{
				idxkid = tree[nn].kids[0];
				d_back_ric_trs_back_funnel0_libstr(nkids, nx[nn], nx[idxkid], nu[nn], nu[idxkid], nb[nn], hidxb[nn], ng[nn], &hsBAbt[idxkid], &hsb[idxkid], &hsrq[nn], &hsDCt[nn], &hsqx[nn], &hsL[nn], &hsLxt[idxkid], compute_Pb, &hsPb[idxkid], &hsux[nn], &hsux[idxkid], work);
				}
			else if(nkids==1) // has one kid => leg
				{
				idxkid = tree[nn].kids[0];
				d_back_ric_trs_back_leg0_libstr(nx[nn], nx[idxkid], nu[nn], nu[idxkid], nb[nn], hidxb[nn], ng[nn], &hsBAbt[idxkid], &hsb[idxkid], &hsrq[nn], &hsDCt[nn], &hsqx[nn], &hsL[nn], &hsLxt[idxkid], compute_Pb, &hsPb[idxkid], &hsux[nn], &hsux[idxkid], work);
				}
			else // has no kids: last stage
				{
				// TODO
				}
			}
		else // kid
			{
			if(nkids>1) // has many kids => funnel
				{
				idxkid = tree[nn].kids[0];
				d_back_ric_trs_back_funnel1_libstr(nkids, nx[nn], nx[idxkid], nu[nn], nu[idxkid], nb[nn], hidxb[nn], ng[nn], &hsBAbt[idxkid], &hsb[idxkid], &hsrq[nn], &hsDCt[nn], &hsqx[nn], &hsL[nn], &hsLxt[idxkid], compute_Pb, &hsPb[idxkid], &hsux[nn], &hsux[idxkid], work);
				}
			else if(nkids==1)// has one kid => leg
				{
				idxkid = tree[nn].kids[0];
				d_back_ric_trs_back_leg1_libstr(nx[nn], nx[idxkid], nu[nn], nu[idxkid], nb[nn], hidxb[nn], ng[nn], &hsBAbt[idxkid], &hsb[idxkid], &hsrq[nn], &hsDCt[nn], &hsqx[nn], &hsL[nn], &hsLxt[idxkid], compute_Pb, &hsPb[idxkid], &hsux[nn], &hsux[idxkid], work);
				}
			else // has no kids: last stage
				{
				d_back_ric_trs_back_legN_libstr(nx[nn], nu[nn], nb[nn], hidxb[nn], ng[nn], &hsrq[nn], &hsDCt[nn], &hsqx[nn], &hsux[nn]);
				}
			}
		}

	// forward substitution

	// process one node at the time, starting from the first one
	for(nn=0; nn<Nn; nn++)
		{
		dad = tree[nn].dad;
		nkids = tree[nn].nkids;
		if(dad<0) // root
			{
			if(nkids>1) // has many kids: funnel
				{
				idxkid = tree[nn].kids[0];
				d_back_ric_trs_forw_funnel0_libstr(nkids, nx[nn], nx[idxkid], nu[nn], nu[idxkid], &hsBAbt[idxkid], &hsb[idxkid], &hsL[nn], &hsLxt[idxkid], &hsux[nn], &hsux[idxkid], compute_pi, &hspi[idxkid], work);
				}
			else if(nkids==1) // has one kid: leg
				{
				idxkid = tree[nn].kids[0];
				d_back_ric_trs_forw_leg0_libstr(nx[nn], nx[idxkid], nu[nn], nu[idxkid], &hsBAbt[idxkid], &hsb[idxkid], &hsL[nn], &hsLxt[idxkid], &hsux[nn], &hsux[idxkid], compute_pi, &hspi[idxkid], work);
				}
			else // no kids
				{
				// TODO
				}
			}
		else // kids
			{
			if(nkids>1) // has many kids: funnel
				{
				idxkid = tree[nn].kids[0];
				d_back_ric_trs_forw_funnel1_libstr(nkids, nx[nn], nx[idxkid], nu[nn], nu[idxkid], &hsBAbt[idxkid], &hsb[idxkid], &hsL[nn], &hsLxt[idxkid], &hsux[nn], &hsux[idxkid], compute_pi, &hspi[idxkid], work);
				}
			else if(nkids==1) // has one kid: leg
				{
				idxkid = tree[nn].kids[0];
				d_back_ric_trs_forw_leg1_libstr(nx[nn], nx[idxkid], nu[nn], nu[idxkid], &hsBAbt[idxkid], &hsb[idxkid], &hsL[nn], &hsLxt[idxkid], &hsux[nn], &hsux[idxkid], compute_pi, &hspi[idxkid], work);
				}
			else // has no kids: last stage
				{
				// nothing to do
				}
			}
		}

	return;

	}



#endif
