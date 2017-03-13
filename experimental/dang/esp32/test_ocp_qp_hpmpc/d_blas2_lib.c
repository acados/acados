/**************************************************************************************************
* acados/external/blasfeo/blas/d_blas2_lib.c                                                                                                *
* This file is part of BLASFEO.                                                                   *
*                                                                                                 *
* BLASFEO -- BLAS For Embedded Optimization.                                                      *
* Copyright (C) 2016-2017 by Gianluca Frison.                                                     *
* Developed at IMTEK (University of Freiburg) under the supervision of Moritz Diehl.              *
* All rights reserved.                                                                            *
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
*                          gianluca.frison (at) imtek.uni-freiburg.de                             *
*                                                                                                 *
**************************************************************************************************/

#include <stdlib.h>
#include <stdio.h>

#if defined(LA_BLAS)
#if defined(LA_BLAS_OPENBLAS)
#include <f77blas.h>
#elif defined(LA_BLAS_BLIS)
#elif defined(LA_BLAS_NETLIB)
#include "d_blas.h"
#elif defined(LA_BLAS_MKL)
#include <mkl_blas.h>
#endif
#endif

#include "blasfeo_common.h"
#include "blasfeo_d_aux.h"


#if defined(LA_REFERENCE)



void dgemv_n_libstr(int m, int n, double alpha, struct d_strmat *sA, int ai, int aj, struct d_strvec *sx, int xi, double beta, struct d_strvec *sy, int yi, struct d_strvec *sz, int zi)
	{
	int ii, jj;
	double
		y_0, y_1, y_2, y_3,
		x_0, x_1;
	int lda = sA->m;
	double *pA = sA->pA + ai + aj*lda;
	double *x = sx->pa + xi;
	double *y = sy->pa + yi;
	double *z = sz->pa + zi;
#if 1 // y reg version
	ii = 0;
	for(; ii<m-1; ii+=2)
		{
		y_0 = 0.0;
		y_1 = 0.0;
		jj = 0;
		for(; jj<n-1; jj+=2)
			{
			y_0 += pA[ii+0+lda*(jj+0)] * x[jj+0] + pA[ii+0+lda*(jj+1)] * x[jj+1];
			y_1 += pA[ii+1+lda*(jj+0)] * x[jj+0] + pA[ii+1+lda*(jj+1)] * x[jj+1];
			}
		if(jj<n)
			{
			y_0 += pA[ii+0+lda*jj] * x[jj];
			y_1 += pA[ii+1+lda*jj] * x[jj];
			}
		z[ii+0] = beta * y[ii+0] + alpha * y_0;
		z[ii+1] = beta * y[ii+1] + alpha * y_1;
		}
	for(; ii<m; ii++)
		{
		y_0 = 0.0;
		for(jj=0; jj<n; jj++)
			{
			y_0 += pA[ii+lda*jj] * x[jj];
			}
		z[ii] = beta * y[ii] + alpha * y_0;
		}
#else // x reg version
	for(ii=0; ii<n; ii++)
		{
		z[ii] = beta * y[ii];
		}
	jj = 0;
	for(; jj<n-1; jj+=2)
		{
		x_0 = alpha * x[jj+0];
		x_1 = alpha * x[jj+1];
		ii = 0;
		for(; ii<m-1; ii+=2)
			{
			z[ii+0] += pA[ii+0+lda*(jj+0)] * x_0 + pA[ii+0+lda*(jj+1)] * x_1;
			z[ii+1] += pA[ii+1+lda*(jj+0)] * x_0 + pA[ii+1+lda*(jj+1)] * x_1;
			}
		for(; ii<m; ii++)
			{
			z[ii] += pA[ii+lda*(jj+0)] * x_0;
			z[ii] += pA[ii+lda*(jj+1)] * x_1;
			}
		}
	for(; jj<n; jj++)
		{
		x_0 = alpha * x[jj+0];
		for(ii=0; ii<m; ii++)
			{
			z[ii] += pA[ii+lda*(jj+0)] * x_0;
			}
		}
#endif
	return;
	}



void dgemv_t_libstr(int m, int n, double alpha, struct d_strmat *sA, int ai, int aj, struct d_strvec *sx, int xi, double beta, struct d_strvec *sy, int yi, struct d_strvec *sz, int zi)
	{
	int ii, jj;
	double
		y_0, y_1;
	int lda = sA->m;
	double *pA = sA->pA + ai + aj*lda;
	double *x = sx->pa + xi;
	double *y = sy->pa + yi;
	double *z = sz->pa + zi;
	jj = 0;
	for(; jj<n-1; jj+=2)
		{
		y_0 = 0.0;
		y_1 = 0.0;
		ii = 0;
		for(; ii<m-1; ii+=2)
			{
			y_0 += pA[ii+0+lda*(jj+0)] * x[ii+0] + pA[ii+1+lda*(jj+0)] * x[ii+1];
			y_1 += pA[ii+0+lda*(jj+1)] * x[ii+0] + pA[ii+1+lda*(jj+1)] * x[ii+1];
			}
		if(ii<m)
			{
			y_0 += pA[ii+lda*(jj+0)] * x[ii];
			y_1 += pA[ii+lda*(jj+1)] * x[ii];
			}
		z[jj+0] = beta * y[jj+0] + alpha * y_0;
		z[jj+1] = beta * y[jj+1] + alpha * y_1;
		}
	for(; jj<n; jj++)
		{
		y_0 = 0.0;
		for(ii=0; ii<m; ii++)
			{
			y_0 += pA[ii+lda*(jj+0)] * x[ii];
			}
		z[jj+0] = beta * y[jj+0] + alpha * y_0;
		}
	return;
	}



// TODO optimize !!!!!
void dgemv_nt_libstr(int m, int n, double alpha_n, double alpha_t, struct d_strmat *sA, int ai, int aj, struct d_strvec *sx_n, int xi_n, struct d_strvec *sx_t, int xi_t, double beta_n, double beta_t, struct d_strvec *sy_n, int yi_n, struct d_strvec *sy_t, int yi_t, struct d_strvec *sz_n, int zi_n, struct d_strvec *sz_t, int zi_t)
	{
	int ii, jj;
	double
		a_00,
		x_n_0,
		y_t_0;
	int lda = sA->m;
	double *pA = sA->pA + ai + aj*lda;
	double *x_n = sx_n->pa + xi_n;
	double *x_t = sx_t->pa + xi_t;
	double *y_n = sy_n->pa + yi_n;
	double *y_t = sy_t->pa + yi_t;
	double *z_n = sz_n->pa + zi_n;
	double *z_t = sz_t->pa + zi_t;
	for(ii=0; ii<m; ii++)
		{
		z_n[ii] = beta_n * y_n[ii];
		}
	for(jj=0; jj<n; jj++)
		{
		y_t_0 = 0.0;
		x_n_0 = alpha_n * x_n[jj];
		for(ii=0; ii<m; ii++)
			{
			a_00 = pA[ii+lda*jj];
			z_n[ii] += a_00 * x_n_0;
			y_t_0 += a_00 * x_t[ii];
			}
		z_t[jj] = beta_t * y_t[jj] + alpha_t * y_t_0;
		}
	return;
	}



// TODO optimize !!!!!
void dsymv_l_libstr(int m, int n, double alpha, struct d_strmat *sA, int ai, int aj, struct d_strvec *sx, int xi, double beta, struct d_strvec *sy, int yi, struct d_strvec *sz, int zi)
	{
	int ii, jj;
	double
		y_0;
	int lda = sA->m;
	double *pA = sA->pA + ai + aj*lda;
	double *x = sx->pa + xi;
	double *y = sy->pa + yi;
	double *z = sz->pa + zi;
	for(ii=0; ii<n; ii++)
		{
		y_0 = 0.0;
		jj = 0;
		for(; jj<=ii; jj++)
			{
			y_0 += pA[ii+lda*jj] * x[jj];
			}
		for( ; jj<m; jj++)
			{
			y_0 += pA[jj+lda*ii] * x[jj];
			}
		z[ii] = beta * y[ii] + alpha * y_0;
		}
	return;
	}



void dtrmv_lnn_libstr(int m, struct d_strmat *sA, int ai, int aj, struct d_strvec *sx, int xi, struct d_strvec *sz, int zi)
	{
	int ii, jj;
	double
		y_0, y_1;
	int lda = sA->m;
	double *pA = sA->pA + ai + aj*lda;
	double *x = sx->pa + xi;
	double *z = sz->pa + zi;
	if(m%2!=0)
		{
		ii = m-1;
		y_0 = x[ii];
		y_0 *= pA[ii+lda*ii];
		for(jj=0; jj<ii; jj++)
			{
			y_0 += pA[ii+lda*jj] * x[jj];
			}
		z[ii] = y_0;
		m -= 1;
		}
	for(ii=m-2; ii>=0; ii-=2)
		{
		y_0 = x[ii+0];
		y_1 = x[ii+1];
		y_1 *= pA[ii+1+lda*(ii+1)];
		y_1 += pA[ii+1+lda*(ii+0)] * y_0;
		y_0 *= pA[ii+0+lda*(ii+0)];
		jj = 0;
		for(; jj<ii-1; jj+=2)
			{
			y_0 += pA[ii+0+lda*(jj+0)] * x[jj+0] + pA[ii+0+lda*(jj+1)] * x[jj+1];
			y_1 += pA[ii+1+lda*(jj+0)] * x[jj+0] + pA[ii+1+lda*(jj+1)] * x[jj+1];
			}
//	XXX there is no clean up loop !!!!!
//		for(; jj<ii; jj++)
//			{
//			y_0 += pA[ii+0+lda*jj] * x[jj];
//			y_1 += pA[ii+1+lda*jj] * x[jj];
//			}
		z[ii+0] = y_0;
		z[ii+1] = y_1;
		}
	return;
	}



void dtrmv_ltn_libstr(int m, struct d_strmat *sA, int ai, int aj, struct d_strvec *sx, int xi, struct d_strvec *sz, int zi)
	{
	int ii, jj;
	double
		y_0, y_1;
	int lda = sA->m;
	double *pA = sA->pA + ai + aj*lda;
	double *x = sx->pa + xi;
	double *z = sz->pa + zi;
	jj = 0;
	for(; jj<m-1; jj+=2)
		{
		y_0 = x[jj+0];
		y_1 = x[jj+1];
		y_0 *= pA[jj+0+lda*(jj+0)];
		y_0 += pA[jj+1+lda*(jj+0)] * y_1;
		y_1 *= pA[jj+1+lda*(jj+1)];
		ii = jj+2;
		for(; ii<m-1; ii+=2)
			{
			y_0 += pA[ii+0+lda*(jj+0)] * x[ii+0] + pA[ii+1+lda*(jj+0)] * x[ii+1];
			y_1 += pA[ii+0+lda*(jj+1)] * x[ii+0] + pA[ii+1+lda*(jj+1)] * x[ii+1];
			}
		for(; ii<m; ii++)
			{
			y_0 += pA[ii+lda*(jj+0)] * x[ii];
			y_1 += pA[ii+lda*(jj+1)] * x[ii];
			}
		z[jj+0] = y_0;
		z[jj+1] = y_1;
		}
	for(; jj<m; jj++)
		{
		y_0 = x[jj];
		y_0 *= pA[jj+lda*jj];
		for(ii=jj+1; ii<m; ii++)
			{
			y_0 += pA[ii+lda*jj] * x[ii];
			}
		z[jj] = y_0;
		}
	return;
	}



void dtrmv_unn_libstr(int m, struct d_strmat *sA, int ai, int aj, struct d_strvec *sx, int xi, struct d_strvec *sz, int zi)
	{
	int ii, jj;
	double
		y_0, y_1,
		x_0, x_1;
	int lda = sA->m;
	double *pA = sA->pA + ai + aj*lda;
	double *x = sx->pa + xi;
	double *z = sz->pa + zi;
#if 1 // y reg version
	jj = 0;
	for(; jj<m-1; jj+=2)
		{
		y_0 = x[jj+0];
		y_1 = x[jj+1];
		y_0 = pA[jj+0+lda*(jj+0)] * y_0;
		y_0 += pA[jj+0+lda*(jj+1)] * y_1;
		y_1 = pA[jj+1+lda*(jj+1)] * y_1;
		ii = jj+2;
		for(; ii<m-1; ii+=2)
			{
			y_0 += pA[jj+0+lda*(ii+0)] * x[ii+0] + pA[jj+0+lda*(ii+1)] * x[ii+1];
			y_1 += pA[jj+1+lda*(ii+0)] * x[ii+0] + pA[jj+1+lda*(ii+1)] * x[ii+1];
			}
		if(ii<m)
			{
			y_0 += pA[jj+0+lda*(ii+0)] * x[ii+0];
			y_1 += pA[jj+1+lda*(ii+0)] * x[ii+0];
			}
		z[jj+0] = y_0;
		z[jj+1] = y_1;
		}
	for(; jj<m; jj++)
		{
		y_0 = pA[jj+lda*jj] * x[jj];
		for(ii=jj+1; ii<m; ii++)
			{
			y_0 += pA[jj+lda*ii] * x[ii];
			}
		z[jj] = y_0;
		}
#else // x reg version
	if(x != z)
		{
		for(ii=0; ii<m; ii++)
			z[ii] = x[ii];
		}
	jj = 0;
	for(; jj<m-1; jj+=2)
		{
		x_0 = z[jj+0];
		x_1 = z[jj+1];
		ii = 0;
		for(; ii<jj-1; ii+=2)
			{
			z[ii+0] += pA[ii+0+lda*(jj+0)] * x_0 + pA[ii+0+lda*(jj+1)] * x_1;
			z[ii+1] += pA[ii+1+lda*(jj+0)] * x_0 + pA[ii+1+lda*(jj+1)] * x_1;
			}
//	XXX there is no clean-up loop, since jj+=2 !!!!!
//		for(; ii<jj; ii++)
//			{
//			z[ii+0] += pA[ii+0+lda*(jj+0)] * x_0 + pA[ii+0+lda*(jj+1)] * x_1;
//			}
		x_0 *= pA[jj+0+lda*(jj+0)];
		x_0 += pA[jj+0+lda*(jj+1)] * x_1;
		x_1 *= pA[jj+1+lda*(jj+1)];
		z[jj+0] = x_0;
		z[jj+1] = x_1;
		}
	for(; jj<m; jj++)
		{
		x_0 = z[jj];
		for(ii=0; ii<jj; ii++)
			{
			z[ii] += pA[ii+lda*jj] * x_0;
			}
		x_0 *= pA[jj+lda*jj];
		z[jj] = x_0;
		}
#endif
	return;
	}



void dtrmv_utn_libstr(int m, struct d_strmat *sA, int ai, int aj, struct d_strvec *sx, int xi, struct d_strvec *sz, int zi)
	{
	int ii, jj;
	double
		y_0, y_1;
	int lda = sA->m;
	double *pA = sA->pA + ai + aj*lda;
	double *x = sx->pa + xi;
	double *z = sz->pa + zi;
	if(m%2!=0)
		{
		jj = m-1;
		y_0 = pA[jj+lda*jj] * x[jj];
		for(ii=0; ii<jj; ii++)
			{
			y_0 += pA[ii+lda*jj] * x[ii];
			}
		z[jj] = y_0;
		m -= 1; // XXX
		}
	for(jj=m-2; jj>=0; jj-=2)
		{
		y_1 = pA[jj+1+lda*(jj+1)] * x[jj+1];
		y_1 += pA[jj+0+lda*(jj+1)] * x[jj+0];
		y_0 = pA[jj+0+lda*(jj+0)] * x[jj+0];
		for(ii=0; ii<jj-1; ii+=2)
			{
			y_0 += pA[ii+0+lda*(jj+0)] * x[ii+0] + pA[ii+1+lda*(jj+0)] * x[ii+1];
			y_1 += pA[ii+0+lda*(jj+1)] * x[ii+0] + pA[ii+1+lda*(jj+1)] * x[ii+1];
			}
//	XXX there is no clean-up loop !!!!!
//		if(ii<jj)
//			{
//			y_0 += pA[ii+lda*(jj+0)] * x[ii];
//			y_1 += pA[ii+lda*(jj+1)] * x[ii];
//			}
		z[jj+0] = y_0;
		z[jj+1] = y_1;
		}
	return;
	}



void dtrsv_lnn_libstr(int m, int n, struct d_strmat *sA, int ai, int aj, struct d_strvec *sx, int xi, struct d_strvec *sz, int zi)
	{
	int ii, jj, j1;
	double
		y_0, y_1,
		x_0, x_1;
	int lda = sA->m;
	double *pA = sA->pA + ai + aj*lda;
	double *dA = sA->dA;
	double *x = sx->pa + xi;
	double *z = sz->pa + zi;
	if(ai==0 & aj==0)
		{
		if(sA->use_dA!=1)
			{
			for(ii=0; ii<n; ii++)
				dA[ii] = 1.0 / pA[ii+lda*ii];
			sA->use_dA = 1;
			}
		}
	else
		{
		for(ii=0; ii<n; ii++)
			dA[ii] = 1.0 / pA[ii+lda*ii];
		sA->use_dA = 0;
		}
#if 1 // y reg version
	ii = 0;
	for(; ii<n-1; ii+=2)
		{
		y_0 = x[ii+0];
		y_1 = x[ii+1];
		jj = 0;
		for(; jj<ii-1; jj+=2)
			{
			y_0 -= pA[ii+0+lda*(jj+0)] * z[jj+0] + pA[ii+0+lda*(jj+1)] * z[jj+1];
			y_1 -= pA[ii+1+lda*(jj+0)] * z[jj+0] + pA[ii+1+lda*(jj+1)] * z[jj+1];
			}
//	XXX there is no clean-up loop !!!!!
//		if(jj<ii)
//			{
//			y_0 -= pA[ii+0+lda*(jj+0)] * z[jj+0];
//			y_1 -= pA[ii+1+lda*(jj+0)] * z[jj+0];
//			}
		y_0 *= dA[ii+0];
		y_1 -= pA[ii+1+lda*(jj+0)] * y_0;
		y_1 *= dA[ii+1];
		z[ii+0] = y_0;
		z[ii+1] = y_1;
		}
	for(; ii<n; ii++)
		{
		y_0 = x[ii];
		for(jj=0; jj<ii; jj++)
			{
			y_0 -= pA[ii+lda*jj] * z[jj];
			}
		y_0 *= dA[ii];
		z[ii] = y_0;
		}
	for(; ii<m-1; ii+=2)
		{
		y_0 = x[ii+0];
		y_1 = x[ii+1];
		jj = 0;
		for(; jj<n-1; jj+=2)
			{
			y_0 -= pA[ii+0+lda*(jj+0)] * z[jj+0] + pA[ii+0+lda*(jj+1)] * z[jj+1];
			y_1 -= pA[ii+1+lda*(jj+0)] * z[jj+0] + pA[ii+1+lda*(jj+1)] * z[jj+1];
			}
		if(jj<n)
			{
			y_0 -= pA[ii+0+lda*(jj+0)] * z[jj+0];
			y_1 -= pA[ii+1+lda*(jj+0)] * z[jj+0];
			}
		z[ii+0] = y_0;
		z[ii+1] = y_1;
		}
	for(; ii<m; ii++)
		{
		y_0 = x[ii];
		for(jj=0; jj<n; jj++)
			{
			y_0 -= pA[ii+lda*jj] * z[jj];
			}
		z[ii] = y_0;
		}
#else // x reg version
	if(x != z)
		{
		for(ii=0; ii<m; ii++)
			z[ii] = x[ii];
		}
	jj = 0;
	for(; jj<n-1; jj+=2)
		{
		x_0 = dA[jj+0] * z[jj+0];
		x_1 = z[jj+1] - pA[jj+1+lda*(jj+0)] * x_0;
		x_1 = dA[jj+1] * x_1;
		z[jj+0] = x_0;
		z[jj+1] = x_1;
		ii = jj+2;
		for(; ii<m-1; ii+=2)
			{
			z[ii+0] -= pA[ii+0+lda*(jj+0)] * x_0 + pA[ii+0+lda*(jj+1)] * x_1;
			z[ii+1] -= pA[ii+1+lda*(jj+0)] * x_0 + pA[ii+1+lda*(jj+1)] * x_1;
			}
		for(; ii<m; ii++)
			{
			z[ii] -= pA[ii+lda*(jj+0)] * x_0 + pA[ii+lda*(jj+1)] * x_1;
			}
		}
	for(; jj<n; jj++)
		{
		x_0 = dA[jj] * z[jj];
		z[jj] = x_0;
		for(ii=jj+1; ii<m; ii++)
			{
			z[ii] -= pA[ii+lda*jj] * x_0;
			}
		}
#endif
	return;
	}



void dtrsv_ltn_libstr(int m, int n, struct d_strmat *sA, int ai, int aj, struct d_strvec *sx, int xi, struct d_strvec *sz, int zi)
	{
	int ii, jj;
	double
		y_0, y_1;
	int lda = sA->m;
	double *pA = sA->pA + ai + aj*lda;
	double *dA = sA->dA;
	double *x = sx->pa + xi;
	double *z = sz->pa + zi;
	if(ai==0 & aj==0)
		{
		if(sA->use_dA!=1)
			{
			for(ii=0; ii<n; ii++)
				dA[ii] = 1.0 / pA[ii+lda*ii];
			sA->use_dA = 1;
			}
		}
	else
		{
		for(ii=0; ii<n; ii++)
			dA[ii] = 1.0 / pA[ii+lda*ii];
		sA->use_dA = 0;
		}
	if(n%2!=0)
		{
		jj = n-1;
		y_0 = x[jj];
		for(ii=jj+1; ii<m; ii++)
			{
			y_0 -= pA[ii+lda*jj] * z[ii];
			}
		y_0 *= dA[jj];
		z[jj] = y_0;
		jj -= 2;
		}
	else
		{
		jj = n-2;
		}
	for(; jj>=0; jj-=2)
		{
		y_0 = x[jj+0];
		y_1 = x[jj+1];
		ii = jj+2;
		for(; ii<m-1; ii+=2)
			{
			y_0 -= pA[ii+0+lda*(jj+0)] * z[ii+0] + pA[ii+1+lda*(jj+0)] * z[ii+1];
			y_1 -= pA[ii+0+lda*(jj+1)] * z[ii+0] + pA[ii+1+lda*(jj+1)] * z[ii+1];
			}
		if(ii<m)
			{
			y_0 -= pA[ii+lda*(jj+0)] * z[ii];
			y_1 -= pA[ii+lda*(jj+1)] * z[ii];
			}
		y_1 *= dA[jj+1];
		y_0 -= pA[jj+1+lda*(jj+0)] * y_1;
		y_0 *= dA[jj+0];
		z[jj+0] = y_0;
		z[jj+1] = y_1;
		}
	return;
	}



#elif defined(LA_BLAS)



void dgemv_n_libstr(int m, int n, double alpha, struct d_strmat *sA, int ai, int aj, struct d_strvec *sx, int xi, double beta, struct d_strvec *sy, int yi, struct d_strvec *sz, int zi)
	{
	char cl = 'l';
	char cn = 'n';
	char cr = 'r';
	char ct = 't';
	char cu = 'u';
	int i1 = 1;
	int lda = sA->m;
	double *pA = sA->pA + ai + aj*lda;
	double *x = sx->pa + xi;
	double *y = sy->pa + yi;
	double *z = sz->pa + zi;
#if defined(LA_BLAS_MKL)
	dcopy(&m, y, &i1, z, &i1);
	dgemv(&cn, &m, &n, &alpha, pA, &lda, x, &i1, &beta, z, &i1);
#else
	dcopy_(&m, y, &i1, z, &i1);
	dgemv_(&cn, &m, &n, &alpha, pA, &lda, x, &i1, &beta, z, &i1);
#endif
	return;
	}



void dgemv_t_libstr(int m, int n, double alpha, struct d_strmat *sA, int ai, int aj, struct d_strvec *sx, int xi, double beta, struct d_strvec *sy, int yi, struct d_strvec *sz, int zi)
	{
	char cl = 'l';
	char cn = 'n';
	char cr = 'r';
	char ct = 't';
	char cu = 'u';
	int i1 = 1;
	int lda = sA->m;
	double *pA = sA->pA + ai + aj*lda;
	double *x = sx->pa + xi;
	double *y = sy->pa + yi;
	double *z = sz->pa + zi;
#if defined(LA_BLAS_MKL)
	dcopy(&n, y, &i1, z, &i1);
	dgemv(&ct, &m, &n, &alpha, pA, &lda, x, &i1, &beta, z, &i1);
#else
	dcopy_(&n, y, &i1, z, &i1);
	dgemv_(&ct, &m, &n, &alpha, pA, &lda, x, &i1, &beta, z, &i1);
#endif
	return;
	}



void dgemv_nt_libstr(int m, int n, double alpha_n, double alpha_t, struct d_strmat *sA, int ai, int aj, struct d_strvec *sx_n, int xi_n, struct d_strvec *sx_t, int xi_t, double beta_n, double beta_t, struct d_strvec *sy_n, int yi_n, struct d_strvec *sy_t, int yi_t, struct d_strvec *sz_n, int zi_n, struct d_strvec *sz_t, int zi_t)
	{
	char cl = 'l';
	char cn = 'n';
	char cr = 'r';
	char ct = 't';
	char cu = 'u';
	int i1 = 1;
	int lda = sA->m;
	double *pA = sA->pA + ai + aj*lda;
	double *x_n = sx_n->pa + xi_n;
	double *x_t = sx_t->pa + xi_t;
	double *y_n = sy_n->pa + yi_n;
	double *y_t = sy_t->pa + yi_t;
	double *z_n = sz_n->pa + zi_n;
	double *z_t = sz_t->pa + zi_t;
#if defined(LA_BLAS_MKL)
	dcopy(&m, y_n, &i1, z_n, &i1);
	dgemv(&cn, &m, &n, &alpha_n, pA, &lda, x_n, &i1, &beta_n, z_n, &i1);
	dcopy(&n, y_t, &i1, z_t, &i1);
	dgemv(&ct, &m, &n, &alpha_t, pA, &lda, x_t, &i1, &beta_t, z_t, &i1);
#else
	dcopy_(&m, y_n, &i1, z_n, &i1);
	dgemv_(&cn, &m, &n, &alpha_n, pA, &lda, x_n, &i1, &beta_n, z_n, &i1);
	dcopy_(&n, y_t, &i1, z_t, &i1);
	dgemv_(&ct, &m, &n, &alpha_t, pA, &lda, x_t, &i1, &beta_t, z_t, &i1);
#endif
	return;
	}



void dsymv_l_libstr(int m, int n, double alpha, struct d_strmat *sA, int ai, int aj, struct d_strvec *sx, int xi, double beta, struct d_strvec *sy, int yi, struct d_strvec *sz, int zi)
	{
	char cl = 'l';
	char cn = 'n';
	char cr = 'r';
	char ct = 't';
	char cu = 'u';
	int i1 = 1;
	double d1 = 1.0;
	int lda = sA->m;
	double *pA = sA->pA + ai + aj*lda;
	double *x = sx->pa + xi;
	double *y = sy->pa + yi;
	double *z = sz->pa + zi;
	int tmp = m-n;
#if defined(LA_BLAS_MKL)
	dcopy(&m, y, &i1, z, &i1);
	dsymv(&cl, &n, &alpha, pA, &lda, x, &i1, &beta, z, &i1);
	dgemv(&cn, &tmp, &n, &alpha, pA+n, &lda, x, &i1, &beta, z+n, &i1);
	dgemv(&ct, &tmp, &n, &alpha, pA+n, &lda, x+n, &i1, &d1, z, &i1);
#else
	dcopy_(&m, y, &i1, z, &i1);
	dsymv_(&cl, &n, &alpha, pA, &lda, x, &i1, &beta, z, &i1);
	dgemv_(&cn, &tmp, &n, &alpha, pA+n, &lda, x, &i1, &beta, z+n, &i1);
	dgemv_(&ct, &tmp, &n, &alpha, pA+n, &lda, x+n, &i1, &d1, z, &i1);
#endif
	return;
	}



void dtrmv_lnn_libstr(int m, struct d_strmat *sA, int ai, int aj, struct d_strvec *sx, int xi, struct d_strvec *sz, int zi)
	{
	char cl = 'l';
	char cn = 'n';
	char cr = 'r';
	char ct = 't';
	char cu = 'u';
	int i1 = 1;
	double d1 = 1.0;
	double dm1 = -1.0;
	int lda = sA->m;
	double *pA = sA->pA + ai + aj*lda;
	double *x = sx->pa + xi;
	double *z = sz->pa + zi;
#if defined(LA_BLAS_MKL)
	dcopy(&m, x, &i1, z, &i1);
	dtrmv(&cl, &cn, &cn, &m, pA, &lda, z, &i1);
#else
	dcopy_(&m, x, &i1, z, &i1);
	dtrmv_(&cl, &cn, &cn, &m, pA, &lda, z, &i1);
#endif
	return;
	}



void dtrmv_ltn_libstr(int m, struct d_strmat *sA, int ai, int aj, struct d_strvec *sx, int xi, struct d_strvec *sz, int zi)
	{
	char cl = 'l';
	char cn = 'n';
	char cr = 'r';
	char ct = 't';
	char cu = 'u';
	int i1 = 1;
	double d1 = 1.0;
	double dm1 = -1.0;
	int lda = sA->m;
	double *pA = sA->pA + ai + aj*lda;
	double *x = sx->pa + xi;
	double *z = sz->pa + zi;
#if defined(LA_BLAS_MKL)
	dcopy(&m, x, &i1, z, &i1);
	dtrmv(&cl, &ct, &cn, &m, pA, &lda, z, &i1);
#else
	dcopy_(&m, x, &i1, z, &i1);
	dtrmv_(&cl, &ct, &cn, &m, pA, &lda, z, &i1);
#endif
	return;
	}



void dtrmv_unn_libstr(int m, struct d_strmat *sA, int ai, int aj, struct d_strvec *sx, int xi, struct d_strvec *sz, int zi)
	{
	char cl = 'l';
	char cn = 'n';
	char cr = 'r';
	char ct = 't';
	char cu = 'u';
	int i1 = 1;
	double d1 = 1.0;
	double dm1 = -1.0;
	int lda = sA->m;
	double *pA = sA->pA + ai + aj*lda;
	double *x = sx->pa + xi;
	double *z = sz->pa + zi;
#if defined(LA_BLAS_MKL)
	dcopy(&m, x, &i1, z, &i1);
	dtrmv(&cu, &cn, &cn, &m, pA, &lda, z, &i1);
#else
	dcopy_(&m, x, &i1, z, &i1);
	dtrmv_(&cu, &cn, &cn, &m, pA, &lda, z, &i1);
#endif
	return;
	}



void dtrmv_utn_libstr(int m, struct d_strmat *sA, int ai, int aj, struct d_strvec *sx, int xi, struct d_strvec *sz, int zi)
	{
	char cl = 'l';
	char cn = 'n';
	char cr = 'r';
	char ct = 't';
	char cu = 'u';
	int i1 = 1;
	double d1 = 1.0;
	double dm1 = -1.0;
	int lda = sA->m;
	double *pA = sA->pA + ai + aj*lda;
	double *x = sx->pa + xi;
	double *z = sz->pa + zi;
#if defined(LA_BLAS_MKL)
	dcopy(&m, x, &i1, z, &i1);
	dtrmv(&cu, &ct, &cn, &m, pA, &lda, z, &i1);
#else
	dcopy_(&m, x, &i1, z, &i1);
	dtrmv_(&cu, &ct, &cn, &m, pA, &lda, z, &i1);
#endif
	return;
	}



void dtrsv_lnn_libstr(int m, int n, struct d_strmat *sA, int ai, int aj, struct d_strvec *sx, int xi, struct d_strvec *sz, int zi)
	{
	char cl = 'l';
	char cn = 'n';
	char cr = 'r';
	char ct = 't';
	char cu = 'u';
	int i1 = 1;
	double d1 = 1.0;
	double dm1 = -1.0;
	int mmn = m-n;
	int lda = sA->m;
	double *pA = sA->pA + ai + aj*lda;
	double *x = sx->pa + xi;
	double *z = sz->pa + zi;
#if defined(LA_BLAS_MKL)
	dcopy(&m, x, &i1, z, &i1);
	dtrsv(&cl, &cn, &cn, &n, pA, &lda, z, &i1);
	dgemv(&cn, &mmn, &n, &dm1, pA+n, &lda, z, &i1, &d1, z+n, &i1);
#else
	dcopy_(&m, x, &i1, z, &i1);
	dtrsv_(&cl, &cn, &cn, &n, pA, &lda, z, &i1);
	dgemv_(&cn, &mmn, &n, &dm1, pA+n, &lda, z, &i1, &d1, z+n, &i1);
#endif
	return;
	}



void dtrsv_ltn_libstr(int m, int n, struct d_strmat *sA, int ai, int aj, struct d_strvec *sx, int xi, struct d_strvec *sz, int zi)
	{
	char cl = 'l';
	char cn = 'n';
	char cr = 'r';
	char ct = 't';
	char cu = 'u';
	int i1 = 1;
	double d1 = 1.0;
	double dm1 = -1.0;
	int mmn = m-n;
	int lda = sA->m;
	double *pA = sA->pA + ai + aj*lda;
	double *x = sx->pa + xi;
	double *z = sz->pa + zi;
#if defined(LA_BLAS_MKL)
	dcopy(&m, x, &i1, z, &i1);
	dgemv(&ct, &mmn, &n, &dm1, pA+n, &lda, z+n, &i1, &d1, z, &i1);
	dtrsv(&cl, &ct, &cn, &n, pA, &lda, z, &i1);
#else
	dcopy_(&m, x, &i1, z, &i1);
	dgemv_(&ct, &mmn, &n, &dm1, pA+n, &lda, z+n, &i1, &d1, z, &i1);
	dtrsv_(&cl, &ct, &cn, &n, pA, &lda, z, &i1);
#endif
	return;
	}



#else

#error : wrong LA choice

#endif
