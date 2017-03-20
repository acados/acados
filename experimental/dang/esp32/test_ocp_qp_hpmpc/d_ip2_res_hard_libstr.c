/**************************************************************************************************
* acados/external/hpmpc/mpc_solvers/d_ip2_res_hard_libstr.c                                                                                                *
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
#include <math.h>
// 2017.03.13 Dang add target.h
#include "target.h"

#ifdef BLASFEO

#include <blasfeo_target.h>
#include <blasfeo_common.h>
#include <blasfeo_d_aux.h>
#include <blasfeo_d_blas.h>

//#else
#include "blas_d.h"
//#endif

#include "aux_d.h"
#include "aux_s.h"
#include "lqcp_solvers.h"
#include "block_size.h"
#include "mpc_aux.h"
#include "d_blas_aux.h"


// use iterative refinement to increase accuracy of the solution of the equality constrained sub-problems
#define ITER_REF 0
#define THR_ITER_REF 1e-5
//#define ITER_REF_REG 0.0
#define CORRECTOR_LOW 1
#define CORRECTOR_HIGH 1



// work space size
int d_ip2_res_mpc_hard_work_space_size_bytes_libstr(int N, int *nx, int *nu, int *nb, int *ng)
	{

	int ii;

	int size = 0;

	for(ii=0; ii<=N; ii++)
		{
		size += d_size_strmat(nu[ii]+nx[ii]+1, nu[ii]+nx[ii]); // L
		size += d_size_strmat(nx[ii], nx[ii]); // Lxt
		size += 5*d_size_strvec(nx[ii]); // b, dpi, Pb, res_b, pi_bkp
		size += 4*d_size_strvec(nu[ii]+nx[ii]); // dux, rq, res_rq, ux_bkp
		size += 8*d_size_strvec(2*nb[ii]+2*ng[ii]); // dlam, dt, tinv, lamt, res_d, res_m, t_bkp, lam_bkp
		size += 2*d_size_strvec(nb[ii]+ng[ii]); // Qx, qx
		}

	// residuals work space size
	size += d_res_res_mpc_hard_work_space_size_bytes_libstr(N, nx, nu, nb, ng);

	// riccati work space size
	size += d_back_ric_rec_work_space_size_bytes_libstr(N, nx, nu, nb, ng);

	// make multiple of (typical) cache line size
	size = (size+63)/64*64;

	return size;
	}



// basic working version

/* primal-dual interior-point method computing residuals at each iteration, hard constraints, time variant matrices, time variant size (mpc version) */
int d_ip2_res_mpc_hard_libstr(int *kk, int k_max, double mu0, double mu_tol, double alpha_min, int warm_start, double *stat, int N, int *nx, int *nu, int *nb, int **idxb, int *ng, struct d_strmat *hsBAbt, struct d_strmat *hsRSQrq, struct d_strmat *hsDCt, struct d_strvec *hsd, struct d_strvec *hsux, int compute_mult, struct d_strvec *hspi, struct d_strvec *hslam, struct d_strvec *hst, void *work)
	{
// printf("d_ip2_res_hard_libstr.c, line 96\n"); // debug
	// indeces
// printf("Free heap size: %d\n",esp_get_free_heap_size()); // for debug
	int jj, ll, ii, it_ref;
// printf("Free heap size: %d\n",esp_get_free_heap_size()); // for debug

	struct d_strmat *hsmatdummy;
	struct d_strvec *hsvecdummy;
// printf("d_ip2_res_hard_libstr.c, line 103\n"); // debug
	// TODO do not use variable size arrays !!!!!
	struct d_strvec hsb[N];
	struct d_strvec hsrq[N+1];
	struct d_strvec hsQx[N+1];
	struct d_strvec hsqx[N+1];
	struct d_strvec hsdux[N+1];
	struct d_strvec hsdpi[N+1];
	struct d_strvec hsdt[N+1];
	struct d_strvec hsdlam[N+1];
	struct d_strvec hstinv[N+1];
	struct d_strvec hslamt[N+1];
	struct d_strvec hsPb[N+1];
	struct d_strmat hsL[N+1];
	struct d_strmat hsLxt[N+1];
	struct d_strvec hsres_rq[N+1];
	struct d_strvec hsres_b[N];
	struct d_strvec hsres_d[N+1];
	struct d_strvec hsres_m[N+1];
	struct d_strmat hsric_work_mat[2];
	struct d_strvec hsric_work_vec[1];
	struct d_strvec hsux_bkp[N+1];
	struct d_strvec hspi_bkp[N+1];
	struct d_strvec hst_bkp[N+1];
	struct d_strvec hslam_bkp[N+1];

	void *d_back_ric_rec_work_space;
	void *d_res_res_mpc_hard_work_space;

	char *c_ptr = work;
// printf("d_ip2_res_hard_libstr.c, line 133\n"); // debug
// printf("Free heap size: %d\n",esp_get_free_heap_size()); // for debug
	// L
	for(ii=0; ii<=N; ii++)
		{
		d_create_strmat(nu[ii]+nx[ii]+1, nu[ii]+nx[ii], &hsL[ii], (void *) c_ptr);
		c_ptr += hsL[ii].memory_size;
		}
// printf("d_ip2_res_hard_libstr.c, line 140\n"); // debug
	// Lxt
	for(ii=0; ii<=N; ii++)
		{
		d_create_strmat(nx[ii], nx[ii], &hsLxt[ii], (void *) c_ptr);
		c_ptr += hsLxt[ii].memory_size;
		}
		// printf("d_ip2_res_hard_libstr.c, line 147\n"); // debug
	// riccati work space
	d_back_ric_rec_work_space = (void *) c_ptr;
	c_ptr += d_back_ric_rec_work_space_size_bytes_libstr(N, nx, nu, nb, ng);
// printf("d_ip2_res_hard_libstr.c, line 151\n"); // debug
	// b as vector
	for(ii=0; ii<N; ii++)
		{
		d_create_strvec(nx[ii+1], &hsb[ii], (void *) c_ptr);
		c_ptr += hsb[ii].memory_size;
		}

	// inputs and states step
	for(ii=0; ii<=N; ii++)
		{
		d_create_strvec(nu[ii]+nx[ii], &hsdux[ii], (void *) c_ptr);
		c_ptr += hsdux[ii].memory_size;
		}

	// equality constr multipliers step
	for(ii=0; ii<=N; ii++)
		{
		d_create_strvec(nx[ii], &hsdpi[ii], (void *) c_ptr);
		c_ptr += hsdpi[ii].memory_size;
		}

	// backup of P*b
	for(ii=0; ii<=N; ii++)
		{
		d_create_strvec(nx[ii], &hsPb[ii], (void *) c_ptr);
		c_ptr += hsPb[ii].memory_size;
		}

	// linear part of cost function
	for(ii=0; ii<=N; ii++)
		{
		d_create_strvec(nu[ii]+nx[ii], &hsrq[ii], (void *) c_ptr);
		c_ptr += hsrq[ii].memory_size;
		}

	// slack variables, Lagrangian multipliers for inequality constraints and work space
	for(ii=0; ii<=N; ii++)
		{
		d_create_strvec(2*nb[ii]+2*ng[ii], &hsdlam[ii], (void *) c_ptr);
		c_ptr += hsdlam[ii].memory_size;
		d_create_strvec(2*nb[ii]+2*ng[ii], &hsdt[ii], (void *) c_ptr);
		c_ptr += hsdt[ii].memory_size;
		}

	for(ii=0; ii<=N; ii++)
		{
		d_create_strvec(2*nb[ii]+2*ng[ii], &hstinv[ii], (void *) c_ptr);
		c_ptr += hstinv[ii].memory_size;
		}

	for(ii=0; ii<=N; ii++)
		{
		d_create_strvec(2*nb[ii]+2*ng[ii], &hslamt[ii], (void *) c_ptr);
		c_ptr += hslamt[ii].memory_size;
		}

	for(ii=0; ii<=N; ii++)
		{
		d_create_strvec(nb[ii]+ng[ii], &hsQx[ii], (void *) c_ptr);
		c_ptr += hsQx[ii].memory_size;
		d_create_strvec(nb[ii]+ng[ii], &hsqx[ii], (void *) c_ptr);
		c_ptr += hsqx[ii].memory_size;
		}
		// printf("d_ip2_res_hard_libstr.c, line 215\n"); // debug
	// residuals work space
	d_res_res_mpc_hard_work_space = (void *) c_ptr;
	c_ptr += d_res_res_mpc_hard_work_space_size_bytes_libstr(N, nx, nu, nb, ng);


	for(ii=0; ii<=N; ii++)
		{
		d_create_strvec(nu[ii]+nx[ii], &hsres_rq[ii], (void *) c_ptr);
		c_ptr += hsres_rq[ii].memory_size;
		}

	for(ii=0; ii<N; ii++)
		{
		d_create_strvec(nx[ii+1], &hsres_b[ii], (void *) c_ptr);
		c_ptr += hsres_b[ii].memory_size;
		}

	for(ii=0; ii<=N; ii++)
		{
		d_create_strvec(2*nb[ii]+2*ng[ii], &hsres_d[ii], (void *) c_ptr);
		c_ptr += hsres_d[ii].memory_size;
		}

	for(ii=0; ii<=N; ii++)
		{
		d_create_strvec(2*nb[ii]+2*ng[ii], &hsres_m[ii], (void *) c_ptr);
		c_ptr += hsres_m[ii].memory_size;
		}

	// backup solution
	for(ii=0; ii<=N; ii++)
		{
		d_create_strvec(nu[ii]+nx[ii], &hsux_bkp[ii], (void *) c_ptr);
		c_ptr += hsux_bkp[ii].memory_size;
		d_create_strvec(nx[ii], &hspi_bkp[ii], (void *) c_ptr);
		c_ptr += hspi_bkp[ii].memory_size;
		d_create_strvec(2*nb[ii]+2*ng[ii], &hslam_bkp[ii], (void *) c_ptr);
		c_ptr += hslam_bkp[ii].memory_size;
		d_create_strvec(2*nb[ii]+2*ng[ii], &hst_bkp[ii], (void *) c_ptr);
		c_ptr += hst_bkp[ii].memory_size;
		}



	// extract b
	for(jj=0; jj<N; jj++)
		{
		drowex_libstr(nx[jj+1], 1.0, &hsBAbt[jj], nu[jj]+nx[jj], 0, &hsb[jj], 0);
		}

	// extract q
	for(jj=0; jj<=N; jj++)
		{
		drowex_libstr(nu[jj]+nx[jj], 1.0, &hsRSQrq[jj], nu[jj]+nx[jj], 0, &hsrq[jj], 0);
		}



	double temp0, temp1;
	double alpha, mu, mu_aff;

	// check if there are inequality constraints
	double mu_scal = 0.0;
	for(jj=0; jj<=N; jj++) mu_scal += 2*nb[jj] + 2*ng[jj];
	if(mu_scal!=0.0) // there are some constraints
		{
		mu_scal = 1.0 / mu_scal;
		}
	else // call the riccati solver and return
		{
		d_back_ric_rec_sv_libstr(N, nx, nu, nb, idxb, ng, 0, hsBAbt, hsb, 0, hsRSQrq, hsrq, hsmatdummy, hsvecdummy, hsvecdummy, hsux, compute_mult, hspi, 1, hsPb, hsL, hsLxt, d_back_ric_rec_work_space);
		// no IPM iterations
		*kk = 0;
		// return success
		return 0;
		}
		// printf("d_ip2_res_hard_libstr.c, line 292\n"); // debug
	//printf("\nmu_scal = %f\n", mu_scal);
	double sigma = 0.0;
	//for(ii=0; ii<=N; ii++)
	//	printf("\n%d %d\n", nb[ii], ng[ii]);
	//exit(1);



	// initialize ux & pi & t>0 & lam>0
	d_init_var_mpc_hard_libstr(N, nx, nu, nb, idxb, ng, hsux, hspi, hsDCt, hsd, hst, hslam, mu0, warm_start);
		// printf("d_ip2_res_hard_libstr.c, line 303\n"); // debug
#if 0
printf("\nux\n");
for(ii=0; ii<=N; ii++)
	d_print_tran_strvec(nu[ii]+nx[ii], &hsux[ii], 0);
printf("\npi\n");
for(ii=0; ii<N; ii++)
	d_print_tran_strvec(nx[ii+1], &hspi[ii+1], 0);
printf("\nlam\n");
for(ii=0; ii<=N; ii++)
	d_print_tran_strvec(2*nb[ii]+2*ng[ii], &hslam[ii], 0);
printf("\nt\n");
for(ii=0; ii<=N; ii++)
	d_print_tran_strvec(2*nb[ii]+2*ng[ii], &hst[ii], 0);
exit(1);
#endif



	// compute the duality gap
	mu = mu0;

	// set to zero iteration count
	*kk = 0;

	// larger than minimum accepted step size
	alpha = 1.0;





	//
	// loop without residuals compuation at early iterations
	//

//	double mu_tol_low = mu_tol;
	double mu_tol_low = mu_tol<THR_ITER_REF ? THR_ITER_REF : mu_tol ;

#if 0
	if(0)
#else
	while( *kk<k_max && mu>mu_tol_low && alpha>=alpha_min )
#endif
		{

//		printf("\nkk = %d (no res)\n", *kk);

		// printf("d_ip2_res_hard_libstr.c, line 351\n"); // debug

		//update cost function matrices and vectors (box constraints)
		d_update_hessian_mpc_hard_libstr(N, nx, nu, nb, ng, hsd, 0.0, hst, hstinv, hslam, hslamt, hsdlam, hsQx, hsqx);

#if 0
for(ii=0; ii<=N; ii++)
	d_print_tran_strvec(2*nb[ii]+2*ng[ii], &hslam[ii], 0);
for(ii=0; ii<=N; ii++)
	d_print_tran_strvec(2*nb[ii]+2*ng[ii], &hst[ii], 0);
for(ii=0; ii<=N; ii++)
	d_print_tran_strvec(nb[ii]+ng[ii], &hsQx[ii], 0);
for(ii=0; ii<=N; ii++)
	d_print_tran_strvec(nb[ii]+ng[ii], &hsqx[ii], 0);
//if(*kk==1)
exit(1);
#endif

		// printf("d_ip2_res_hard_libstr.c, line 369\n"); // debug
		// compute the search direction: factorize and solve the KKT system
#if 1
		d_back_ric_rec_sv_libstr(N, nx, nu, nb, idxb, ng, 0, hsBAbt, hsb, 1, hsRSQrq, hsrq, hsDCt, hsQx, hsqx, hsdux, compute_mult, hsdpi, 1, hsPb, hsL, hsLxt, d_back_ric_rec_work_space);
#else
		d_back_ric_rec_trf_tv_res(N, nx, nu, pBAbt, pQ, pL, dL, work, nb, idxb, ng, pDCt, Qx, bd);
		d_back_ric_rec_trs_tv_res(N, nx, nu, pBAbt, b, pL, dL, q, l, dux, work, 1, Pb, compute_mult, dpi, nb, idxb, ng, pDCt, qx);
#endif


#if 0
for(ii=0; ii<=N; ii++)
	d_print_strmat(nu[ii]+nx[ii]+1, nu[ii]+nx[ii], &hsRSQrq[ii], 0, 0);
//exit(1);
#endif
#if 0
for(ii=0; ii<=N; ii++)
	d_print_pmat(pnz[ii], cnux[ii], hpL[ii], cnux[ii]);
if(*kk==2)
exit(1);
#endif
#if 0
printf("\ndux\n");
for(ii=0; ii<=N; ii++)
	d_print_tran_strvec(nu[ii]+nx[ii], &hsdux[ii], 0);
printf("\ndpi\n");
for(ii=1; ii<=N; ii++)
	d_print_tran_strvec(nx[ii], &hsdpi[ii], 0);
//if(*kk==2)
exit(1);
#endif


#if CORRECTOR_LOW==1 // IPM1

		// compute t_aff & dlam_aff & dt_aff & alpha
		alpha = 1.0;
		d_compute_alpha_mpc_hard_libstr(N, nx, nu, nb, idxb, ng, &alpha, hst, hsdt, hslam, hsdlam, hslamt, hsdux, hsDCt, hsd);



		stat[5*(*kk)] = sigma;
		stat[5*(*kk)+1] = alpha;

		alpha *= 0.995;

#if 0
printf("\nalpha = %f\n", alpha);
exit(1);
#endif

		// printf("d_ip2_res_hard_libstr.c, line 420\n"); // debug
		// compute the affine duality gap
		d_compute_mu_mpc_hard_libstr(N, nx, nu, nb, ng, &mu_aff, mu_scal, alpha, hslam, hsdlam, hst, hsdt);

		stat[5*(*kk)+2] = mu_aff;

#if 0
printf("\nmu = %f\n", mu_aff);
exit(1);
#endif



		// compute sigma
		sigma = mu_aff/mu;
		sigma = sigma*sigma*sigma;
//		if(sigma<sigma_min)
//			sigma = sigma_min;
//printf("\n%f %f %f %f\n", mu_aff, mu, sigma, mu_scal);
//exit(1);

#if 0
for(ii=0; ii<=N; ii++)
	d_print_mat(1, 2*pnb[ii], dt[ii], 1);
for(ii=0; ii<=N; ii++)
	d_print_mat(1, 2*pnb[ii], dlam[ii], 1);
for(ii=0; ii<=N; ii++)
	d_print_mat(1, 2*pnb[ii], t_inv[ii], 1);
for(ii=0; ii<=N; ii++)
	d_print_mat(1, nb[ii], pl[ii], 1);
//exit(1);
#endif

		// printf("d_ip2_res_hard_libstr.c, line 453\n"); // debug
//		d_update_gradient_mpc_hard_tv(N, nx, nu, nb, ng, sigma*mu, dt, dlam, t_inv, pl, qx);
		d_update_gradient_mpc_hard_libstr(N, nx, nu, nb, ng, sigma*mu, hsdt, hsdlam, hstinv, hsqx);

#if 0
for(ii=0; ii<=N; ii++)
	d_print_mat(1, nb[ii]+ng[ii], qx[ii], 1);
//if(*kk==1)
exit(1);
#endif



//		// copy b into x
//		for(ii=0; ii<N; ii++)
//			for(jj=0; jj<nx[ii+1]; jj++)
//				dux[ii+1][nu[ii+1]+jj] = pBAbt[ii][(nu[ii]+nx[ii])/bs*bs*cnx[ii+1]+(nu[ii]+nx[ii])%bs+bs*jj]; // copy b


		// printf("d_ip2_res_hard_libstr.c, line 472\n"); // debug
		// solve the system
		d_back_ric_rec_trs_libstr(N, nx, nu, nb, idxb, ng, hsBAbt, hsb, hsrq, hsDCt, hsqx, hsdux, compute_mult, hsdpi, 0, hsPb, hsL, hsLxt, d_back_ric_rec_work_space);

#if 0
printf("\ndux\n");
for(ii=0; ii<=N; ii++)
	d_print_tran_strvec(nu[ii]+nx[ii], &hsdux[ii], 0);
//if(*kk==1)
exit(1);
#endif



#endif // end of IPM1

		// printf("d_ip2_res_hard_libstr.c, line 488\n"); // debug
		// compute t & dlam & dt & alpha
		alpha = 1.0;
		d_compute_alpha_mpc_hard_libstr(N, nx, nu, nb, idxb, ng, &alpha, hst, hsdt, hslam, hsdlam, hslamt, hsdux, hsDCt, hsd);

#if 0
printf("\nalpha = %f\n", alpha);
printf("\ndt\n");
for(ii=0; ii<=N; ii++)
	d_print_tran_strvec(2*nb[ii]+2*ng[ii], &hsdt[ii], 0);
printf("\ndlam\n");
for(ii=0; ii<=N; ii++)
	d_print_tran_strvec(2*nb[ii]+2*ng[ii], &hsdlam[ii], 0);
exit(2);
#endif

		stat[5*(*kk)] = sigma;
		stat[5*(*kk)+3] = alpha;

		alpha *= 0.995;


#if 0
printf("\ndux\n");
for(ii=0; ii<=N; ii++)
	d_print_tran_strvec(nu[ii]+nx[ii], &hsdux[ii], 0);
printf("\ndpi\n");
for(ii=1; ii<=N; ii++)
	d_print_tran_strvec(nx[ii], &hsdpi[ii], 0);
printf("\ndlam\n");
for(ii=0; ii<=N; ii++)
	d_print_tran_strvec(2*nb[ii]+2*ng[ii], &hsdlam[ii], 0);
printf("\ndt\n");
for(ii=0; ii<=N; ii++)
	d_print_tran_strvec(2*nb[ii]+2*ng[ii], &hsdt[ii], 0);
//if(*kk==1)
exit(1);
#endif

		// printf("d_ip2_res_hard_libstr.c, line 527\n"); // debug
		// compute step & update x, u, lam, t & compute the duality gap mu
		d_update_var_mpc_hard_libstr(N, nx, nu, nb, ng, &mu, mu_scal, alpha, hsux, hsdux, hst, hsdt, hslam, hsdlam, hspi, hsdpi);

		stat[5*(*kk)+4] = mu;


#if 0
printf("\nux\n");
for(ii=0; ii<=N; ii++)
	d_print_tran_strvec(nu[ii]+nx[ii], &hsux[ii], 0);
printf("\npi\n");
for(ii=1; ii<=N; ii++)
	d_print_tran_strvec(nx[ii], &hspi[ii], 0);
printf("\nlam\n");
for(ii=0; ii<=N; ii++)
	d_print_tran_strvec(2*nb[ii]+2*ng[ii], &hslam[ii], 0);
printf("\nt\n");
for(ii=0; ii<=N; ii++)
	d_print_tran_strvec(2*nb[ii]+2*ng[ii], &hst[ii], 0);
//if(*kk==1)
//exit(1);
#endif


		// increment loop index
		(*kk)++;


		} // end of IP loop



#if 0
printf("\nux\n");
for(jj=0; jj<=N; jj++)
	d_print_tran_strvec(nu[jj]+nx[jj], &hsux[jj], 0);
printf("\npi\n");
for(jj=1; jj<=N; jj++)
	d_print_tran_strvec(nx[jj], &hspi[jj], 0);
printf("\nlam\n");
for(ii=0; ii<=N; ii++)
	d_print_tran_strvec(2*nb[ii]+2*ng[ii], &hslam[ii], 0);
printf("\nt\n");
for(ii=0; ii<=N; ii++)
	d_print_tran_strvec(2*nb[ii]+2*ng[ii], &hst[ii], 0);
exit(2);
#endif


	//
	// loop with residuals computation and iterative refinement for high-accuracy result
	//
		// printf("d_ip2_res_hard_libstr.c, line 580\n"); // debug
	// compute residuals
	d_res_res_mpc_hard_libstr(N, nx, nu, nb, idxb, ng, hsBAbt, hsb, hsRSQrq, hsrq, hsux, hsDCt, hsd, hspi, hslam, hst, hsres_rq, hsres_b, hsres_d, hsres_m, &mu, d_res_res_mpc_hard_work_space);

#if 0
	printf("kk = %d\n", *kk);
	printf("\nres_q\n");
	for(jj=0; jj<=N; jj++)
		d_print_e_tran_strvec(nu[jj]+nx[jj], &hsres_rq[jj], 0);
	printf("\nres_b\n");
	for(jj=0; jj<N; jj++)
		d_print_e_tran_strvec(nx[jj+1], &hsres_b[jj], 0);
	printf("\nres_d\n");
	for(jj=0; jj<=N; jj++)
		d_print_e_tran_strvec(2*nb[jj]+2*ng[jj], &hsres_d[jj], 0);
	printf("\nres_m\n");
	for(jj=0; jj<=N; jj++)
		d_print_e_tran_strvec(2*nb[jj]+2*ng[jj], &hsres_m[jj], 0);
	printf("\nmu\n");
	d_print_e_mat(1, 1, &mu, 1);
	exit(2);
#endif





	// IP loop
#if 0
	int ipm_it;
	for(ipm_it=0; ipm_it<3; ipm_it++)
#else
	while( *kk<k_max && mu>mu_tol && alpha>=alpha_min ) // XXX exit conditions on residuals???
#endif
		{

//		printf("\nkk = %d (res)\n", *kk);


#if 0
printf("\nIPM it %d\n", *kk);
#endif


		// printf("d_ip2_res_hard_libstr.c, line 624\n"); // debug
		// compute the update of Hessian and gradient from box and general constraints
		d_update_hessian_gradient_res_mpc_hard_libstr(N, nx, nu, nb, ng, hsres_d, hsres_m, hst, hslam, hstinv, hsQx, hsqx);

#if 0
for(ii=0; ii<=N; ii++)
	d_print_tran_strvec(nb[ii]+ng[ii], &hsQx[ii], 0);
for(ii=0; ii<=N; ii++)
	d_print_tran_strvec(nb[ii]+ng[ii], &hsqx[ii], 0);
//if(*kk==1)
exit(1);
#endif


		// printf("d_ip2_res_hard_libstr.c, line 638\n"); // debug
		// compute the search direction: factorize and solve the KKT system
#if ITER_REF>0
// TODO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		// update Hessian and gradient
		for(ii=0; ii<=N; ii++)
			{

			// box constraints
			// gradient
			for(jj=0; jj<nu[ii]+nx[ii]; jj++)
				q2[ii][jj] = res_q[ii][jj];
			dvecad_libsp(nb[ii], idxb[ii], 1.0, qx[ii], q2[ii]);
//			d_print_mat_e(1, nu[ii]+nx[ii], q2[ii], 1);
			// hessian
			for(jj=0; jj<pnz[ii]*cnux[ii]; jj++)
				pQ2[ii][jj] = pQ[ii][jj];
			ddiaad_libsp(nb[ii], idxb[ii], 1.0, Qx[ii], pQ2[ii], cnux[ii]);
#ifdef BLASFEO
			drowin_lib(cnux[ii], 1.0, q2[ii], pQ2[ii]+(nu[ii]+nx[ii])/bs*bs*cnux[ii]+(nu[ii]+nx[ii])%bs);
#else
			drowin_lib(cnux[ii], q2[ii], pQ2[ii]+(nu[ii]+nx[ii])/bs*bs*cnux[ii]+(nu[ii]+nx[ii])%bs);
#endif
//			drowin_lib(cnux[ii], res_q[ii], pQ2[ii]+(nu[ii]+nx[ii])/bs*bs*cnux[ii]+(nu[ii]+nx[ii])%bs);
//			drowad_libsp(nb[ii], idxb[ii], 1.0, qx[ii], pQ2[ii]+(nu[ii]+nx[ii])/bs*bs*cnux[ii]+(nu[ii]+nx[ii])%bs);
//			d_print_pmat_e(nu[ii]+nx[ii]+1, nu[ii]+nx[ii], bs, pQ2[ii], cnux[ii]);

			// general constraints
			if(ng[ii]>0) // TODO unsymmetric update not requiring sqrt & div ???
				{
				work2 = work + pnz[ii]*cng[ii];
//				for(jj=0; jj<ng[ii]; jj++)
//					Qx[ii][pnb[ii]+jj] = sqrt(Qx[ii][pnb[ii]+jj]); // XXX
				dgemm_diag_right_lib(nu[ii]+nx[ii], ng[ii], pDCt[ii], cng[ii], Qx[ii]+pnb[ii], 0, work, cng[ii], work, cng[ii]);
#ifdef BLASFEO
				drowin_lib(ng[ii], 1.0, qx[ii]+pnb[ii], work+(nu[ii]+nx[ii])/bs*cng[ii]*bs+(nu[ii]+nx[ii])%bs);
#else
				drowin_lib(ng[ii], qx[ii]+pnb[ii], work+(nu[ii]+nx[ii])/bs*cng[ii]*bs+(nu[ii]+nx[ii])%bs);
#endif
//				for(jj=0; jj<ng[ii]; jj++)
//					work[(nu[ii]+nx[ii])/bs*cng[ii]*bs+(nu[ii]+nx[ii])%bs+jj*bs] /= Qx[ii][pnb[ii]+jj];
#ifdef BLASFEO
				dgecp_lib(nu[ii]+nx[ii], 1.0, ng[ii], 0, pDCt[ii], cng[ii], 0, work2, png[ii]);
				dsyrk_nt_l_lib(nu[ii]+nx[ii]+1, nu[ii]+nx[ii], ng[ii], 1.0, work, cng[ii], work2, cng[ii], 1.0, pQ2[ii], cnux[ii], pQ2[ii], cnux[ii]);
				drowex_lib(cnux[ii], 1.0, pQ2[ii]+(nu[ii]+nx[ii])/bs*bs*cnux[ii]+(nu[ii]+nx[ii])%bs, q2[ii]);
#else
				dgecp_lib(nu[ii]+nx[ii], ng[ii], 0, pDCt[ii], cng[ii], 0, work2, png[ii]);
				dsyrk_nt_lib(nu[ii]+nx[ii]+1, nu[ii]+nx[ii], ng[ii], work, cng[ii], work2, cng[ii], 1, pQ2[ii], cnux[ii], pQ2[ii], cnux[ii]);
				drowex_lib(cnux[ii], pQ2[ii]+(nu[ii]+nx[ii])/bs*bs*cnux[ii]+(nu[ii]+nx[ii])%bs, q2[ii]);
#endif
				}

//			d_print_pmat_e(nu[ii]+nx[ii]+1, nu[ii]+nx[ii], bs, pQ2[ii], cnux[ii]);

//			// regularization
//			ddiareg_lib(nu[ii]+nx[ii], ITER_REF_REG, 0, pQ2[ii], cnux[ii]);

			}
//		exit(2);
		// printf("d_ip2_res_hard_libstr.c, line 697\n"); // debug
		// factorize & solve KKT system
		d_back_ric_rec_sv_libstr(N, nx, nu, nb2, idxb, ng2, 1, hsBAbt, hsres_b, 0, hsRSQrq2, hsres_q, hsmatdummy, hsvecdummy, hsvecdummy, hsdux, compute_mult, hsdpi, 1, hsPb, hsL, hsLxt, d_back_ric_rec_work_space);
		// printf("d_ip2_res_hard_libstr.c, line 700\n"); // debug
#if CORRECTOR_HIGH==1
		if(0)
#else
		for(it_ref=0; it_ref<ITER_REF; it_ref++)
#endif
			{

//			// remove regularization
//			for(ii=0; ii<=N; ii++)
//				ddiareg_lib(nu[ii]+nx[ii], -ITER_REF_REG, 0, pQ2[ii], cnux[ii]);
		// printf("d_ip2_res_hard_libstr.c, line 711\n"); // debug
			// compute residuals
			d_res_res_mpc_hard_libstr(N, nx, nu, nb2, idxb, ng2, hsBAbt, hsres_b, hsRSQrq2, hsrq2, hsdux, hsmatdummy, hsvecdummy, hsdpi, hsvecdummy, hsvecdummy, hsres_q2, hsres_b2, hsvecdummy, hsvecdummy, pdummy, d_res_res_mpc_hard_work_space);

#if 0
			printf("\niterative refinemet %d\n", it_ref);
			printf("\nres_q2\n");
			for(ii=0; ii<=N; ii++)
				d_print_mat_e(1, nu[ii]+nx[ii], res_q2[ii], 1);
			printf("\nres_b2\n");
			for(ii=1; ii<=N; ii++)
				d_print_mat_e(1, nx[ii], res_b2[ii], 1);
//			exit(2);
#endif
		// printf("d_ip2_res_hard_libstr.c, line 725\n"); // debug
			// solve for residuals
//			d_back_ric_rec_trs_tv_res(N, nx, nu, nb2, idxb, ng2, pBAbt, res_b2, res_q2, ppdummy, ppdummy, dux2, compute_mult, dpi2, 1, Pb2, memory, work);
			d_back_ric_rec_trs_libstr(N, nx, nu, nb2, idxb, ng2, hsBAbt, hsres_b2, hsres_q2, hsmatdummy, hsvecdummy, hsdux2, compute_mult, hsdpi, 1, hsPb, hsL, hsLxt, d_back_ric_rec_work_space);

	//		printf("\nux2\n");
	//		for(ii=0; ii<=N; ii++)
	//			d_print_mat_e(1, nu[ii]+nx[ii], dux2[ii], 1);
	//		printf("\npi2\n");
	//		for(ii=0; ii<N; ii++)
	//			d_print_mat_e(1, nx[ii+1], dpi2[ii], 1);
	//		exit(2);

			// update solution
			for(ii=0; ii<=N; ii++)
				for(jj=0; jj<nu[ii]+nx[ii]; jj++)
					dux[ii][jj] += dux2[ii][jj];
			for(ii=1; ii<=N; ii++)
				for(jj=0; jj<nx[ii]; jj++)
					dpi[ii][jj] += dpi2[ii][jj];

			}

#if 0
		// compute residuals again
		d_res_res_mpc_hard_libstr(N, nx, nu, nb2, idxb, ng2, hsBAbt, hsres_b, hsRSQrq2, hsrq2, hsdux, hsmatdummy, hsvecdummy, hsdpi, hsvecdummy, hsvecdummy, hsres_q2, hsres_b2, hsvecdummy, hsvecdummy, pdummy, d_res_res_mpc_hard_work_space);

		printf("\nres_q2\n");
		for(ii=0; ii<=N; ii++)
			d_print_mat_e(1, nu[ii]+nx[ii], res_q2[ii], 1);
		printf("\nres_b2\n");
		for(ii=0; ii<N; ii++)
			d_print_mat_e(1, nx[ii+1], res_b2[ii], 1);
//		exit(2);
#endif


#else // no iterative refinement
#if 0
for(ii=0; ii<=N; ii++)
	d_print_e_tran_strvec(nu[ii]+nx[ii], &hsres_rq[ii], 0);
for(ii=0; ii<N; ii++)
	d_print_e_tran_strvec(nx[ii+1], &hsres_b[ii], 0);
for(ii=0; ii<=N; ii++)
	d_print_tran_strvec(nb[ii]+ng[ii], &hsQx[ii], 0);
for(ii=0; ii<=N; ii++)
	d_print_tran_strvec(nb[ii]+ng[ii], &hsqx[ii], 0);
for(ii=0; ii<=N; ii++)
	d_print_strmat(nu[ii]+nx[ii], nu[ii]+nx[ii], &hsRSQrq[ii], 0, 0);
for(ii=0; ii<N; ii++)
	d_print_strmat(nu[ii]+nx[ii], nx[ii+1], &hsBAbt[ii+1], 0, 0);
exit(1);
#endif
#if 1
		d_back_ric_rec_sv_libstr(N, nx, nu, nb, idxb, ng, 1, hsBAbt, hsres_b, 1, hsRSQrq, hsres_rq, hsDCt, hsQx, hsqx, hsdux, compute_mult, hsdpi, 1, hsPb, hsL, hsLxt, d_back_ric_rec_work_space);
#else
		d_back_ric_rec_trf_tv_res(N, nx, nu, pBAbt, pQ, pL, dL, work, nb, idxb, ng, pDCt, Qx, bd);
		d_back_ric_rec_trs_tv_res(N, nx, nu, pBAbt, res_b, pL, dL, res_q, l, dux, work, 1, Pb, compute_mult, dpi, nb, idxb, ng, pDCt, qx);
#endif
#endif


#if 0
//printf("\npL\n");
//for(ii=0; ii<=N; ii++)
//	d_print_pmat(nu[ii]+nx[ii]+1, nu[ii]+nx[ii]+1, bs, pL[ii], cnux[ii]);
printf("\ndL\n");
for(ii=0; ii<=N; ii++)
	d_print_mat_e(1, nu[ii]+nx[ii], dL[ii], 1);
//exit(1);
#endif
#if 0
printf("\ndux\n");
for(ii=0; ii<=N; ii++)
	d_print_tran_strvec(nu[ii]+nx[ii], &hsdux[ii], 0);
printf("\ndpi\n");
for(ii=0; ii<=N; ii++)
	d_print_tran_strvec(nx[ii], &hsdpi[ii], 0);
//if(*kk==1)
exit(1);
#endif



#if 0
for(ii=0; ii<=N; ii++)
	d_print_pmat(nu[ii]+nx[ii]+1, nu[ii]+nx[ii]+1, bs, pQ[ii], cnux[ii]);
//exit(1);
#endif
#if 0
for(ii=0; ii<=N; ii++)
	d_print_pmat(nu[ii]+nx[ii]+1, nu[ii]+nx[ii], bs, pL[ii], cnux[ii]);
//exit(1);
#endif
#if 0
printf("\nux_aff\n");
for(ii=0; ii<=N; ii++)
	d_print_mat(1, nu[ii]+nx[ii], dux[ii], 1);
printf("\npi_aff\n");
for(ii=1; ii<=N; ii++)
	d_print_mat(1, nx[ii], dpi[ii], 1);
//if(*kk==1)
exit(1);
#endif


#if CORRECTOR_HIGH==1 // IPM1
		// printf("d_ip2_res_hard_libstr.c, line 832\n"); // debug
		// compute t_aff & dlam_aff & dt_aff & alpha
		alpha = 1.0;
		d_compute_alpha_res_mpc_hard_libstr(N, nx, nu, nb, idxb, ng, hsdux, hst, hstinv, hslam, hsDCt, hsres_d, hsres_m, hsdt, hsdlam, &alpha);



		stat[5*(*kk)] = sigma;
		stat[5*(*kk)+1] = alpha;

		alpha *= 0.995;

#if 0
printf("\nalpha = %f\n", alpha);
exit(1);
#endif

		// printf("d_ip2_res_hard_libstr.c, line 849\n"); // debug
		// compute the affine duality gap
		d_compute_mu_res_mpc_hard_libstr(N, nx, nu, nb, ng, alpha, hslam, hsdlam, hst, hsdt, &mu_aff, mu_scal);

		stat[5*(*kk)+2] = mu_aff;

#if 0
printf("\nmu = %f\n", mu_aff);
exit(1);
#endif



		// compute sigma
		sigma = mu_aff/mu;
		sigma = sigma*sigma*sigma;
//		if(sigma<sigma_min)
//			sigma = sigma_min;
//printf("\n%f %f %f %f\n", mu_aff, mu, sigma, mu_scal);
//exit(1);

#if 0
for(ii=0; ii<=N; ii++)
	d_print_mat(1, 2*pnb[ii], dt[ii], 1);
for(ii=0; ii<=N; ii++)
	d_print_mat(1, 2*pnb[ii], dlam[ii], 1);
for(ii=0; ii<=N; ii++)
	d_print_mat(1, 2*pnb[ii], t_inv[ii], 1);
for(ii=0; ii<=N; ii++)
	d_print_mat(1, nb[ii], pl[ii], 1);
//exit(1);
#endif

		// printf("d_ip2_res_hard_libstr.c, line 882\n"); // debug
		// update res_m
		d_compute_centering_correction_res_mpc_hard_libstr(N, nb, ng, sigma*mu, hsdt, hsdlam, hsres_m);

		// printf("d_ip2_res_hard_libstr.c, line 886\n"); // debug

		// update gradient
		d_update_gradient_res_mpc_hard_libstr(N, nx, nu, nb, ng, hsres_d, hsres_m, hslam, hstinv, hsqx);

#if 0
for(ii=0; ii<=N; ii++)
	d_print_mat(1, nb[ii]+ng[ii], qx[ii], 1);
if(*kk==1)
exit(1);
#endif



#if ITER_REF>0

// TODO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		// update gradient
		for(ii=0; ii<=N; ii++)
			{
			// copy gradient
			for(jj=0; jj<nu[ii]+nx[ii]; jj++)
				q2[ii][jj] = res_q[ii][jj];
			// box constraints
			if(nb[ii]>0)
				dvecad_libsp(nb[ii], idxb[ii], 1.0, qx[ii], q2[ii]);
			// general constraints
			if(ng[ii]>0)
				dgemv_n_libstr(nu[ii]+nx[ii], ng[ii], 1.0, &hsDCt[ii], 0, 0, &hsqx[ii], nb[ii], 1.0, &hsrq2[ii], 0, &hsrq2[ii], 0);
			}
		// printf("d_ip2_res_hard_libstr.c, line 916\n"); // debug
		// solve the KKT system
//		d_back_ric_rec_trs_libstr(N, nx, nu, nb2, idxb, ng2, hsBAbt, hsres_b, hsq2, hsvecdummy, hsvecdummy, hsdux, compute_mult, hsdpi, 0, hsPb, memory, work);
		d_back_ric_rec_trs_libstr(N, nx, nu, nb2, idxb, ng2, hsBAbt, hsres_b, hsrq2, hsmatdummy, hsvecdummy, hsdux, compute_mult, hsdpi, 0, hsPb, hsL, hsLxt, d_back_ric_rec_work_space);

#if 0
printf("\ndux\n");
for(ii=0; ii<=N; ii++)
	d_print_mat(1, nu[ii]+nx[ii], dux[ii], 1);
//if(*kk==1)
exit(1);
#endif



#if 0
		if(0)
#else
		for(it_ref=0; it_ref<ITER_REF; it_ref++)
#endif
			{

//			// remove regularization
//			for(ii=0; ii<=N; ii++)
//				ddiareg_lib(nu[ii]+nx[ii], -ITER_REF_REG, 0, pQ2[ii], cnux[ii]);

			// compute residuals
			d_res_res_mpc_hard_libstr(N, nx, nu, nb2, idxb, ng2, hsBAbt, hsres_b, hsRSQrq2, hsrq2, hsdux, hsmatdummy, hsvecdummy, hsdpi, hsvecdummy, hsvecdummy, hsres_q2, hsres_b2, hsvecdummy, hsvecdummy, pdummy, d_res_res_mpc_hard_work_space);

#if 0
			printf("\niterative refinemet %d\n", it_ref);
			printf("\nres_q2\n");
			for(ii=0; ii<=N; ii++)
				d_print_mat_e(1, nu[ii]+nx[ii], res_q2[ii], 1);
			printf("\nres_b2\n");
			for(ii=0; ii<N; ii++)
				d_print_mat_e(1, nx[ii+1], res_b2[ii], 1);
//			exit(2);
#endif

			// solve for residuals
//			d_back_ric_rec_trs_tv_res(N, nx, nu, nb2, idxb, ng2, pBAbt, res_b2, res_q2, ppdummy, ppdummy, dux2, compute_mult, dpi2, 1, Pb2, memory, work);
			d_back_ric_rec_trs_libstr(N, nx, nu, nb2, idxb, ng2, hsBAbt, hsres_b, hsrq2, hsmatdummy, hsvecdummy, hsdux2, compute_mult, hsdpi2, 0, hsPb2, hsL, hsLxt, d_back_ric_rec_work_space);

	//		printf("\nux2\n");
	//		for(ii=0; ii<=N; ii++)
	//			d_print_mat_e(1, nu[ii]+nx[ii], dux2[ii], 1);
	//		printf("\npi2\n");
	//		for(ii=0; ii<N; ii++)
	//			d_print_mat_e(1, nx[ii+1], dpi2[ii], 1);
	//		exit(2);

			// update solution
			for(ii=0; ii<=N; ii++)
				for(jj=0; jj<nu[ii]+nx[ii]; jj++)
					dux[ii][jj] += dux2[ii][jj];
			for(ii=1; ii<=N; ii++)
				for(jj=0; jj<nx[ii]; jj++)
					dpi[ii][jj] += dpi2[ii][jj];

			}

#if 0
		// compute residuals again
		d_res_res_mpc_hard_tv(N, nx, nu, nb2, idxb, ng2, pBAbt, res_b, pQ2, q2, dux, ppdummy, ppdummy, dpi, ppdummy, ppdummy, res_work, res_q2, res_b2, ppdummy, ppdummy, pdummy);

		printf("\nres_q2\n");
		for(ii=0; ii<=N; ii++)
			d_print_mat_e(1, nu[ii]+nx[ii], res_q2[ii], 1);
		printf("\nres_b2\n");
		for(ii=0; ii<N; ii++)
			d_print_mat_e(1, nx[ii+1], res_b2[ii], 1);
//		exit(2);
#endif


#else // no iter res

		// solve the KKT system
//		d_back_ric_rec_trs_tv_res(N, nx, nu, pBAbt, res_b, pL, dL, res_q, l, dux, work, 0, Pb, compute_mult, dpi, nb, idxb, ng, pDCt, qx);
		d_back_ric_rec_trs_libstr(N, nx, nu, nb, idxb, ng, hsBAbt, hsres_b, hsres_rq, hsDCt, hsqx, hsdux, compute_mult, hsdpi, 0, hsPb, hsL, hsLxt, d_back_ric_rec_work_space);


#endif



#endif // end of IPM1

		// printf("d_ip2_res_hard_libstr.c, line 1005\n"); // debug
		// compute t & dlam & dt & alpha
		alpha = 1.0;
		d_compute_alpha_res_mpc_hard_libstr(N, nx, nu, nb, idxb, ng, hsdux, hst, hstinv, hslam, hsDCt, hsres_d, hsres_m, hsdt, hsdlam, &alpha);

#if 0
printf("\nalpha = %f\n", alpha);
printf("\nd\n");
for(ii=0; ii<=N; ii++)
	d_print_tran_strvec(2*nb[ii]+2*ng[ii], &hsd[ii], 0);
printf("\nres_d\n");
for(ii=0; ii<=N; ii++)
	d_print_tran_strvec(2*nb[ii]+2*ng[ii], &hsres_d[ii], 0);
printf("\ndt\n");
for(ii=0; ii<=N; ii++)
	d_print_tran_strvec(2*nb[ii]+2*ng[ii], &hsdt[ii], 0);
printf("\ndlam\n");
for(ii=0; ii<=N; ii++)
	d_print_tran_strvec(2*nb[ii]+2*ng[ii], &hsdlam[ii], 0);
exit(2);
#endif

		stat[5*(*kk)] = sigma;
		stat[5*(*kk)+3] = alpha;

		alpha *= 0.995;


		// printf("d_ip2_res_hard_libstr.c, line 1033\n"); // debug
		// backup & update x, u, pi, lam, t
		d_backup_update_var_res_mpc_hard_libstr(N, nx, nu, nb, ng, alpha, hsux_bkp, hsux, hsdux, hspi_bkp, hspi, hsdpi, hst_bkp, hst, hsdt, hslam_bkp, hslam, hsdlam);
//		d_update_var_res_mpc_hard_tv(N, nx, nu, nb, ng, alpha, ux, dux, pi, dpi, t, dt, lam, dlam);


#if 0
printf("\nux\n");
for(ii=0; ii<=N; ii++)
	d_print_tran_strvec(nu[ii]+nx[ii], &hsux[ii], 0);
printf("\npi\n");
for(ii=1; ii<=N; ii++)
	d_print_tran_strvec(nx[ii], &hspi[ii], 0);
printf("\nlam\n");
for(ii=0; ii<=N; ii++)
	d_print_tran_strvec(2*nb[ii]+2*ng[ii], &hslam[ii], 0);
printf("\nt\n");
for(ii=0; ii<=N; ii++)
	d_print_tran_strvec(2*nb[ii]+2*ng[ii], &hst[ii], 0);
//if(*kk==1)
//exit(1);
#endif


		// restore dynamics
		for(jj=0; jj<N; jj++)
			drowin_libstr(nx[jj+1], 1.0, &hsb[jj], 0, &hsBAbt[jj], nu[jj]+nx[jj], 0);
		// printf("d_ip2_res_hard_libstr.c, line 1060\n"); // debug
		d_res_res_mpc_hard_libstr(N, nx, nu, nb, idxb, ng, hsBAbt, hsb, hsRSQrq, hsrq, hsux, hsDCt, hsd, hspi, hslam, hst, hsres_rq, hsres_b, hsres_d, hsres_m, &mu, d_res_res_mpc_hard_work_space);

#if 0
	printf("\nres_q\n");
	for(jj=0; jj<=N; jj++)
		d_print_e_tran_strvec(nu[jj]+nx[jj], &hsres_rq[jj], 0);
	printf("\nres_b\n");
	for(jj=0; jj<N; jj++)
		d_print_e_tran_strvec(nx[jj+1], &hsres_b[jj], 0);
	printf("\nres_d\n");
	for(jj=0; jj<=N; jj++)
		d_print_e_tran_strvec(2*nb[jj]+2*ng[jj], &hsres_d[jj], 0);
	printf("\nres_m\n");
	for(jj=0; jj<=N; jj++)
		d_print_e_tran_strvec(2*nb[jj]+2*ng[jj], &hsres_m[jj], 0);
	printf("\nmu\n");
	d_print_e_mat(1, 1, &mu, 1);
	exit(2);
#endif

		stat[5*(*kk)+4] = mu;




		// increment loop index
		(*kk)++;


		} // end of IP loop



	// restore Hessian XXX and gradient
//	for(jj=0; jj<=N; jj++)
//		{
//		for(ll=0; ll<nb[jj]; ll++)
//			{
//			idx = idxb[jj][ll];
//			pQ[jj][idx/bs*bs*cnux[jj]+idx%bs+idx*bs] = bd[jj][ll];
//			pQ[jj][(nu[jj]+nx[jj])/bs*bs*cnux[jj]+(nu[jj]+nx[jj])%bs+idx*bs] = bl[jj][ll];
//			}
//		}

#if 0
printf("\nABb\n");
for(jj=0; jj<N; jj++)
	d_print_pmat(nu[jj]+nx[jj]+1, nx[jj+1], bs, pBAbt[jj], cnx[jj+1]);
printf("\nQ\n");
for(jj=0; jj<=N; jj++)
	d_print_pmat(nu[jj]+nx[jj]+1, nu[jj]+nx[jj], bs, pQ[jj], cnux[jj]);
//printf("\nux\n");
//for(jj=0; jj<=N; jj++)
//	d_print_mat(1, nu[jj]+nx[jj], ux[jj], 1);
//exit(2);
#endif

#if 0
printf("\nux\n");
for(jj=0; jj<=N; jj++)
	d_print_mat(1, nu[jj]+nx[jj], ux[jj], 1);
printf("\npi\n");
for(jj=0; jj<N; jj++)
	d_print_mat(1, nx[jj+1], pi[jj], 1);
printf("\nlam\n");
for(ii=0; ii<=N; ii++)
	d_print_mat(1, 2*pnb[ii]+2*png[ii], lam[ii], 1);
printf("\nt\n");
for(ii=0; ii<=N; ii++)
	d_print_mat(1, 2*pnb[ii]+2*png[ii], t[ii], 1);
exit(2);
#endif



	// TODO if mu is nan, recover solution !!!
//	if(mu==1.0/0.0 || mu==-1.0/0.0)
//		{
//		printf("\nnan!!!\n");
//		exit(3);
//		}








//exit(2);

	// successful exit
	if(mu<=mu_tol)
		return 0;

	// max number of iterations reached
	if(*kk>=k_max)
		return 1;

	// no improvement
	if(alpha<alpha_min)
		return 2;

	// impossible
	return -1;

	} // end of ipsolver



#endif
