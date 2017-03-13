/**************************************************************************************************
* acados/external/hpmpc/mpc_solvers/c99/d_res_ip_res_hard_libstr.c                                                                                                *
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

// 2017.03.13 Dang add target.h
#include "target.h"

#ifdef BLASFEO

#include <blasfeo_target.h>
#include <blasfeo_common.h>
#include <blasfeo_d_blas.h>
#include <blasfeo_d_aux.h>



int d_res_res_mpc_hard_work_space_size_bytes_libstr(int N, int *nx, int *nu, int *nb, int *ng)
	{

	int ii;

	int size = 0;

	int ngM = 0;
	for(ii=0; ii<=N; ii++)
		{
		ngM = ng[ii]>ngM ? ng[ii] : ngM;
		}

	size += 2*d_size_strvec(ngM); // res_work[0], res_work[1]

	// make multiple of (typical) cache line size
	size = (size+63)/64*64;

	return size;

	}



void d_res_res_mpc_hard_libstr(int N, int *nx, int *nu, int *nb, int **idxb, int *ng, struct d_strmat *hsBAbt, struct d_strvec *hsb, struct d_strmat *hsQ, struct d_strvec *hsq, struct d_strvec *hsux, struct d_strmat *hsDCt, struct d_strvec *hsd, struct d_strvec *hspi, struct d_strvec *hslam, struct d_strvec *hst, struct d_strvec *hsrq, struct d_strvec *hsrb, struct d_strvec *hsrd, struct d_strvec *hsrm, double *mu, void *work)
	{

	int ii, jj;

	char *c_ptr;

	struct d_strvec hswork_0, hswork_1;
	double *work0, *work1;


	double
		*ptr_b, *ptr_q, *ptr_d, *ptr_ux, *ptr_pi, *ptr_lam, *ptr_t, *ptr_rb, *ptr_rq, *ptr_rd, *ptr_rm;

	int
		*ptr_idxb;

	int nu0, nu1, nx0, nx1, nxm, nb0, ng0, nb_tot;

	double
		mu2;

	// initialize mu
	nb_tot = 0;
	mu2 = 0;



	// first stage
	ii = 0;
	nu0 = nu[ii];
	nu1 = nu[ii+1];
	nx0 = nx[ii]; // nx1;
	nx1 = nx[ii+1];
	nb0 = nb[ii];
	ng0 = ng[ii];

	ptr_b = hsb[ii].pa;
	ptr_q = hsq[ii].pa;
	ptr_ux = hsux[ii].pa;
	ptr_pi = hspi[ii].pa;
	ptr_rq = hsrq[ii].pa;
	ptr_rb = hsrb[ii].pa;

	if(nb0>0 | ng0>0)
		{
		ptr_d = hsd[ii].pa;
		ptr_lam = hslam[ii].pa;
		ptr_t = hst[ii].pa;
		ptr_rd = hsrd[ii].pa;
		ptr_rm = hsrm[ii].pa;
		}

	for(jj=0; jj<nu0+nx0; jj++)
		ptr_rq[jj] = ptr_q[jj];

	if(nb0>0)
		{

		ptr_idxb = idxb[ii];
		nb_tot += nb0;

		for(jj=0; jj<nb0; jj++)
			{
			ptr_rq[ptr_idxb[jj]] += - ptr_lam[jj] + ptr_lam[nb0+jj];

			ptr_rd[jj]     = ptr_d[jj]     - ptr_ux[ptr_idxb[jj]] + ptr_t[jj];
			ptr_rd[nb0+jj] = ptr_d[nb0+jj] - ptr_ux[ptr_idxb[jj]] - ptr_t[nb0+jj];

			ptr_rm[jj]     = ptr_lam[jj]     * ptr_t[jj];
			ptr_rm[nb0+jj] = ptr_lam[nb0+jj] * ptr_t[nb0+jj];
			mu2 += ptr_rm[jj] + ptr_rm[nb0+jj];
			}
		}

	dsymv_l_libstr(nu0+nx0, nu0+nx0, 1.0, &hsQ[ii], 0, 0, &hsux[ii], 0, 1.0, &hsrq[ii], 0, &hsrq[ii], 0);

	ptr_ux = hsux[ii+1].pa;
	for(jj=0; jj<nx1; jj++)
		ptr_rb[jj] = ptr_b[jj] - ptr_ux[nu1+jj];

	dgemv_nt_libstr(nu0+nx0, nx1, 1.0, 1.0, &hsBAbt[ii], 0, 0, &hspi[ii+1], 0, &hsux[ii], 0, 1.0, 1.0, &hsrq[ii], 0, &hsrb[ii], 0, &hsrq[ii], 0, &hsrb[ii], 0);

	if(ng0>0)
		{

		c_ptr = (char *) work;
		d_create_strvec(ng0, &hswork_0, (void *) c_ptr);
		c_ptr += hswork_0.memory_size;
		d_create_strvec(ng0, &hswork_1, (void *) c_ptr);
		c_ptr += hswork_1.memory_size;
		work0 = hswork_0.pa;
		work1 = hswork_1.pa;

		ptr_d   += 2*nb0;
		ptr_lam += 2*nb0;
		ptr_t   += 2*nb0;
		ptr_rd  += 2*nb0;
		ptr_rm  += 2*nb0;

		nb_tot += ng0;

		for(jj=0; jj<ng0; jj++)
			{
			work0[jj] = ptr_lam[jj+ng0] - ptr_lam[jj+0];

			ptr_rd[jj]     = ptr_d[jj]     + ptr_t[jj];
			ptr_rd[ng0+jj] = ptr_d[ng0+jj] - ptr_t[ng0+jj];

			ptr_rm[jj]     = ptr_lam[jj]     * ptr_t[jj];
			ptr_rm[ng0+jj] = ptr_lam[ng0+jj] * ptr_t[ng0+jj];
			mu2 += ptr_rm[jj] + ptr_rm[ng0+jj];
			}

		dgemv_nt_libstr(nu0+nx0, ng0, 1.0, 1.0, &hsDCt[ii], 0, 0, &hswork_0, 0, &hsux[ii], 0, 1.0, 0.0, &hsrq[ii], 0, &hswork_1, 0, &hsrq[ii], 0, &hswork_1, 0);

		for(jj=0; jj<ng0; jj++)
			{
			ptr_rd[jj]     -= work1[jj];
			ptr_rd[ng0+jj] -= work1[jj];
			}

		}



	// middle stages
	for(ii=1; ii<N; ii++)
		{

		nu0 = nu1;
		nu1 = nu[ii+1];
		nx0 = nx1;
		nx1 = nx[ii+1];
		nb0 = nb[ii];
		ng0 = ng[ii];

		ptr_b = hsb[ii].pa;
		ptr_q = hsq[ii].pa;
		ptr_ux = hsux[ii].pa;
		ptr_pi = hspi[ii].pa;
		ptr_rq = hsrq[ii].pa;
		ptr_rb = hsrb[ii].pa;

		if(nb0>0 | ng0>0)
			{
			ptr_d = hsd[ii].pa;
			ptr_lam = hslam[ii].pa;
			ptr_t = hst[ii].pa;
			ptr_rd = hsrd[ii].pa;
			ptr_rm = hsrm[ii].pa;
			}

		for(jj=0; jj<nu0; jj++)
			ptr_rq[jj] = ptr_q[jj];

		for(jj=0; jj<nx0; jj++)
			ptr_rq[nu0+jj] = ptr_q[nu0+jj] - ptr_pi[jj];

		if(nb0>0)
			{

			ptr_idxb = idxb[ii];
			nb_tot += nb0;

			for(jj=0; jj<nb0; jj++)
				{
				ptr_rq[ptr_idxb[jj]] += - ptr_lam[jj] + ptr_lam[nb0+jj];

				ptr_rd[jj]     = ptr_d[jj]     - ptr_ux[ptr_idxb[jj]] + ptr_t[jj];
				ptr_rd[nb0+jj] = ptr_d[nb0+jj] - ptr_ux[ptr_idxb[jj]] - ptr_t[nb0+jj];

				ptr_rm[jj]     = ptr_lam[jj]     * ptr_t[jj];
				ptr_rm[nb0+jj] = ptr_lam[nb0+jj] * ptr_t[nb0+jj];
				mu2 += ptr_rm[jj] + ptr_rm[nb0+jj];
				}
			}

		dsymv_l_libstr(nu0+nx0, nu0+nx0, 1.0, &hsQ[ii], 0, 0, &hsux[ii], 0, 1.0, &hsrq[ii], 0, &hsrq[ii], 0);

		ptr_ux = hsux[ii+1].pa;
		for(jj=0; jj<nx1; jj++)
			ptr_rb[jj] = ptr_b[jj] - ptr_ux[nu1+jj];

		dgemv_nt_libstr(nu0+nx0, nx1, 1.0, 1.0, &hsBAbt[ii], 0, 0, &hspi[ii+1], 0, &hsux[ii], 0, 1.0, 1.0, &hsrq[ii], 0, &hsrb[ii], 0, &hsrq[ii], 0, &hsrb[ii], 0);

		if(ng0>0)
			{

			c_ptr = (char *) work;
			d_create_strvec(ng0, &hswork_0, (void *) c_ptr);
			c_ptr += hswork_0.memory_size;
			d_create_strvec(ng0, &hswork_1, (void *) c_ptr);
			c_ptr += hswork_1.memory_size;
			work0 = hswork_0.pa;
			work1 = hswork_1.pa;

			ptr_d   += 2*nb0;
			ptr_lam += 2*nb0;
			ptr_t   += 2*nb0;
			ptr_rd  += 2*nb0;
			ptr_rm  += 2*nb0;

			nb_tot += ng0;

			for(jj=0; jj<ng0; jj++)
				{
				work0[jj] = ptr_lam[jj+ng0] - ptr_lam[jj+0];

				ptr_rd[jj]     = ptr_d[jj]     + ptr_t[jj];
				ptr_rd[ng0+jj] = ptr_d[ng0+jj] - ptr_t[ng0+jj];

				ptr_rm[jj]     = ptr_lam[jj]     * ptr_t[jj];
				ptr_rm[ng0+jj] = ptr_lam[ng0+jj] * ptr_t[ng0+jj];
				mu2 += ptr_rm[jj] + ptr_rm[ng0+jj];
				}

			dgemv_nt_libstr(nu0+nx0, ng0, 1.0, 1.0, &hsDCt[ii], 0, 0, &hswork_0, 0, &hsux[ii], 0, 1.0, 0.0, &hsrq[ii], 0, &hswork_1, 0, &hsrq[ii], 0, &hswork_1, 0);

			for(jj=0; jj<ng0; jj++)
				{
				ptr_rd[jj]     -= work1[jj];
				ptr_rd[ng0+jj] -= work1[jj];
				}

			}

		}



	// last stage
	ii = N;
	nu0 = nu1;
	nx0 = nx1;
	nb0 = nb[ii];
	ng0 = ng[ii];

	ptr_q = hsq[ii].pa;
	ptr_ux = hsux[ii].pa;
	ptr_pi = hspi[ii].pa;
	ptr_rq = hsrq[ii].pa;

	if(nb0>0 | ng0>0)
		{
		ptr_d = hsd[ii].pa;
		ptr_lam = hslam[ii].pa;
		ptr_t = hst[ii].pa;
		ptr_rd = hsrd[ii].pa;
		ptr_rm = hsrm[ii].pa;
		}

	for(jj=0; jj<nx0; jj++)
		ptr_rq[nu0+jj] = - ptr_pi[jj] + ptr_q[nu0+jj];

	if(nb0>0)
		{

		ptr_idxb = idxb[ii];
		nb_tot += nb0;

		for(jj=0; jj<nb0; jj++)
			{
			ptr_rq[ptr_idxb[jj]] += - ptr_lam[jj] + ptr_lam[nb0+jj];

			ptr_rd[jj]     = ptr_d[jj]     - ptr_ux[ptr_idxb[jj]] + ptr_t[jj];
			ptr_rd[nb0+jj] = ptr_d[nb0+jj] - ptr_ux[ptr_idxb[jj]] - ptr_t[nb0+jj];

			ptr_rm[jj]     = ptr_lam[jj]     * ptr_t[jj];
			ptr_rm[nb0+jj] = ptr_lam[nb0+jj] * ptr_t[nb0+jj];
			mu2 += ptr_rm[jj] + ptr_rm[nb0+jj];
			}
		}

	dsymv_l_libstr(nx0, nx0, 1.0, &hsQ[ii], 0, 0, &hsux[ii], 0, 1.0, &hsrq[ii], 0, &hsrq[ii], 0);

	if(ng0>0)
		{

		c_ptr = (char *) work;
		d_create_strvec(ng0, &hswork_0, (void *) c_ptr);
		c_ptr += hswork_0.memory_size;
		d_create_strvec(ng0, &hswork_1, (void *) c_ptr);
		c_ptr += hswork_1.memory_size;
		work0 = hswork_0.pa;
		work1 = hswork_1.pa;

		ptr_d   += 2*nb0;
		ptr_lam += 2*nb0;
		ptr_t   += 2*nb0;
		ptr_rd  += 2*nb0;
		ptr_rm  += 2*nb0;

		nb_tot += ng0;

		for(jj=0; jj<ng0; jj++)
			{
			work0[jj] = ptr_lam[jj+ng0] - ptr_lam[jj+0];

			ptr_rd[jj]     = ptr_d[jj]     + ptr_t[jj];
			ptr_rd[ng0+jj] = ptr_d[ng0+jj] - ptr_t[ng0+jj];

			ptr_rm[jj]     = ptr_lam[jj]     * ptr_t[jj];
			ptr_rm[ng0+jj] = ptr_lam[ng0+jj] * ptr_t[ng0+jj];
			mu2 += ptr_rm[jj] + ptr_rm[ng0+jj];
			}

		dgemv_nt_libstr(nu0+nx0, ng0, 1.0, 1.0, &hsDCt[ii], 0, 0, &hswork_0, 0, &hsux[ii], 0, 1.0, 0.0, &hsrq[ii], 0, &hswork_1, 0, &hsrq[ii], 0, &hswork_1, 0);

		for(jj=0; jj<ng0; jj++)
			{
			ptr_rd[jj]     -= work1[jj];
			ptr_rd[ng0+jj] -= work1[jj];
			}

		}


	// normalize mu
	if(nb_tot!=0)
		{
		mu2 /= 2.0*nb_tot;
		mu[0] = mu2;
		}



	return;

	}



#endif
