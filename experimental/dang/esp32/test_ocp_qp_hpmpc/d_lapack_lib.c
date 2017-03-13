/**************************************************************************************************
* acados/external/blasfeo/blas/d_lapack_lib.c                                                                                                *
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
#include <math.h>

#if defined(LA_BLAS)
#if defined(LA_BLAS_OPENBLAS)
#include <f77blas.h>
#elif defined(LA_BLAS_BLIS)
#elif defined(LA_BLAS_NETLIB)
#include "d_blas.h"
#elif defined(LA_BLAS_MKL)
#include <mkl_blas.h>
#include <mkl_lapack.h>
#endif
#endif

#include "blasfeo_common.h"
#include "blasfeo_d_aux.h"



#if defined(LA_REFERENCE)



// dpotrf
void dpotrf_l_libstr(int m, int n, struct d_strmat *sC, int ci, int cj, struct d_strmat *sD, int di, int dj)
	{
	int ii, jj, kk;
	double
		f_00_inv,
		f_10, f_11_inv,
		c_00, c_01,
		c_10, c_11;
	int ldc = sC->m;
	int ldd = sD->m;
	double *pC = sC->pA + ci + cj*ldc;
	double *pD = sD->pA + di + dj*ldd;
	double *dD = sD->dA;
	if(di==0 & dj==0)
		sD->use_dA = 1;
	else
		sD->use_dA = 0;
	jj = 0;
	for(; jj<n-1; jj+=2)
		{
		// factorize diagonal
		c_00 = pC[jj+0+ldc*(jj+0)];;
		c_10 = pC[jj+1+ldc*(jj+0)];;
		c_11 = pC[jj+1+ldc*(jj+1)];;
		for(kk=0; kk<jj; kk++)
			{
			c_00 -= pD[jj+0+ldd*kk] * pD[jj+0+ldd*kk];
			c_10 -= pD[jj+1+ldd*kk] * pD[jj+0+ldd*kk];
			c_11 -= pD[jj+1+ldd*kk] * pD[jj+1+ldd*kk];
			}
		if(c_00>0)
			{
			f_00_inv = 1.0/sqrt(c_00);
			}
		else
			{
			f_00_inv = 0.0;
			}
		dD[jj+0] = f_00_inv;
		pD[jj+0+ldd*(jj+0)] = c_00 * f_00_inv;
		f_10 = c_10 * f_00_inv;
		pD[jj+1+ldd*(jj+0)] = f_10;
		c_11 -= f_10 * f_10;
		if(c_11>0)
			{
			f_11_inv = 1.0/sqrt(c_11);
			}
		else
			{
			f_11_inv = 0.0;
			}
		dD[jj+1] = f_11_inv;
		pD[jj+1+ldd*(jj+1)] = c_11 * f_11_inv;
		// solve lower
		ii = jj+2;
		for(; ii<m-1; ii+=2)
			{
			c_00 = pC[ii+0+ldc*(jj+0)];
			c_10 = pC[ii+1+ldc*(jj+0)];
			c_01 = pC[ii+0+ldc*(jj+1)];
			c_11 = pC[ii+1+ldc*(jj+1)];
			for(kk=0; kk<jj; kk++)
				{
				c_00 -= pD[ii+0+ldd*kk] * pD[jj+0+ldd*kk];
				c_10 -= pD[ii+1+ldd*kk] * pD[jj+0+ldd*kk];
				c_01 -= pD[ii+0+ldd*kk] * pD[jj+1+ldd*kk];
				c_11 -= pD[ii+1+ldd*kk] * pD[jj+1+ldd*kk];
				}
			c_00 *= f_00_inv;
			c_10 *= f_00_inv;
			pD[ii+0+ldd*(jj+0)] = c_00;
			pD[ii+1+ldd*(jj+0)] = c_10;
			c_01 -= c_00 * f_10;
			c_11 -= c_10 * f_10;
			pD[ii+0+ldd*(jj+1)] = c_01 * f_11_inv;
			pD[ii+1+ldd*(jj+1)] = c_11 * f_11_inv;
			}
		for(; ii<m; ii++)
			{
			c_00 = pC[ii+0+ldc*(jj+0)];
			c_01 = pC[ii+0+ldc*(jj+1)];
			for(kk=0; kk<jj; kk++)
				{
				c_00 -= pD[ii+0+ldd*kk] * pD[jj+0+ldd*kk];
				c_01 -= pD[ii+0+ldd*kk] * pD[jj+1+ldd*kk];
				}
			c_00 *= f_00_inv;
			pD[ii+0+ldd*(jj+0)] = c_00;
			c_01 -= c_00 * f_10;
			pD[ii+0+ldd*(jj+1)] = c_01 * f_11_inv;
			}
		}
	for(; jj<n; jj++)
		{
		// factorize diagonal
		c_00 = pC[jj+ldc*jj];;
		for(kk=0; kk<jj; kk++)
			{
			c_00 -= pD[jj+ldd*kk] * pD[jj+ldd*kk];
			}
		if(c_00>0)
			{
			f_00_inv = 1.0/sqrt(c_00);
			}
		else
			{
			f_00_inv = 0.0;
			}
		dD[jj] = f_00_inv;
		pD[jj+ldd*jj] = c_00 * f_00_inv;
		// solve lower
		for(ii=jj+1; ii<m; ii++)
			{
			c_00 = pC[ii+ldc*jj];
			for(kk=0; kk<jj; kk++)
				{
				c_00 -= pD[ii+ldd*kk] * pD[jj+ldd*kk];
				}
			pD[ii+ldd*jj] = c_00 * f_00_inv;
			}
		}
	return;
	}



// dsyrk dpotrf
void dsyrk_dpotrf_ln_libstr(int m, int n, int k, struct d_strmat *sA, int ai, int aj, struct d_strmat *sB, int bi, int bj, struct d_strmat *sC, int ci, int cj, struct d_strmat *sD, int di, int dj)
	{
	int ii, jj, kk;
	double
		f_00_inv,
		f_10, f_11_inv,
		c_00, c_01,
		c_10, c_11;
	int lda = sA->m;
	int ldb = sB->m;
	int ldc = sC->m;
	int ldd = sD->m;
	double *pA = sA->pA + ai + aj*lda;
	double *pB = sB->pA + bi + bj*ldb;
	double *pC = sC->pA + ci + cj*ldc;
	double *pD = sD->pA + di + dj*ldd;
	double *dD = sD->dA;
	if(di==0 & dj==0)
		sD->use_dA = 1;
	else
		sD->use_dA = 0;
	jj = 0;
	for(; jj<n-1; jj+=2)
		{
		// factorize diagonal
		c_00 = pC[jj+0+ldc*(jj+0)];;
		c_10 = pC[jj+1+ldc*(jj+0)];;
		c_11 = pC[jj+1+ldc*(jj+1)];;
		for(kk=0; kk<k; kk++)
			{
			c_00 += pA[jj+0+lda*kk] * pB[jj+0+ldb*kk];
			c_10 += pA[jj+1+lda*kk] * pB[jj+0+ldb*kk];
			c_11 += pA[jj+1+lda*kk] * pB[jj+1+ldb*kk];
			}
		for(kk=0; kk<jj; kk++)
			{
			c_00 -= pD[jj+0+ldd*kk] * pD[jj+0+ldd*kk];
			c_10 -= pD[jj+1+ldd*kk] * pD[jj+0+ldd*kk];
			c_11 -= pD[jj+1+ldd*kk] * pD[jj+1+ldd*kk];
			}
		if(c_00>0)
			{
			f_00_inv = 1.0/sqrt(c_00);
			}
		else
			{
			f_00_inv = 0.0;
			}
		dD[jj+0] = f_00_inv;
		pD[jj+0+ldd*(jj+0)] = c_00 * f_00_inv;
		f_10 = c_10 * f_00_inv;
		pD[jj+1+ldd*(jj+0)] = f_10;
		c_11 -= f_10 * f_10;
		if(c_11>0)
			{
			f_11_inv = 1.0/sqrt(c_11);
			}
		else
			{
			f_11_inv = 0.0;
			}
		dD[jj+1] = f_11_inv;
		pD[jj+1+ldd*(jj+1)] = c_11 * f_11_inv;
		// solve lower
		ii = jj+2;
		for(; ii<m-1; ii+=2)
			{
			c_00 = pC[ii+0+ldc*(jj+0)];
			c_10 = pC[ii+1+ldc*(jj+0)];
			c_01 = pC[ii+0+ldc*(jj+1)];
			c_11 = pC[ii+1+ldc*(jj+1)];
			for(kk=0; kk<k; kk++)
				{
				c_00 += pA[ii+0+lda*kk] * pB[jj+0+ldb*kk];
				c_10 += pA[ii+1+lda*kk] * pB[jj+0+ldb*kk];
				c_01 += pA[ii+0+lda*kk] * pB[jj+1+ldb*kk];
				c_11 += pA[ii+1+lda*kk] * pB[jj+1+ldb*kk];
				}
			for(kk=0; kk<jj; kk++)
				{
				c_00 -= pD[ii+0+ldd*kk] * pD[jj+0+ldd*kk];
				c_10 -= pD[ii+1+ldd*kk] * pD[jj+0+ldd*kk];
				c_01 -= pD[ii+0+ldd*kk] * pD[jj+1+ldd*kk];
				c_11 -= pD[ii+1+ldd*kk] * pD[jj+1+ldd*kk];
				}
			c_00 *= f_00_inv;
			c_10 *= f_00_inv;
			pD[ii+0+ldd*(jj+0)] = c_00;
			pD[ii+1+ldd*(jj+0)] = c_10;
			c_01 -= c_00 * f_10;
			c_11 -= c_10 * f_10;
			pD[ii+0+ldd*(jj+1)] = c_01 * f_11_inv;
			pD[ii+1+ldd*(jj+1)] = c_11 * f_11_inv;
			}
		for(; ii<m; ii++)
			{
			c_00 = pC[ii+0+ldc*(jj+0)];
			c_01 = pC[ii+0+ldc*(jj+1)];
			for(kk=0; kk<k; kk++)
				{
				c_00 += pA[ii+0+lda*kk] * pB[jj+0+ldb*kk];
				c_01 += pA[ii+0+lda*kk] * pB[jj+1+ldb*kk];
				}
			for(kk=0; kk<jj; kk++)
				{
				c_00 -= pD[ii+0+ldd*kk] * pD[jj+0+ldd*kk];
				c_01 -= pD[ii+0+ldd*kk] * pD[jj+1+ldd*kk];
				}
			c_00 *= f_00_inv;
			pD[ii+0+ldd*(jj+0)] = c_00;
			c_01 -= c_00 * f_10;
			pD[ii+0+ldd*(jj+1)] = c_01 * f_11_inv;
			}
		}
	for(; jj<n; jj++)
		{
		// factorize diagonal
		c_00 = pC[jj+ldc*jj];;
		for(kk=0; kk<k; kk++)
			{
			c_00 += pA[jj+lda*kk] * pB[jj+ldb*kk];
			}
		for(kk=0; kk<jj; kk++)
			{
			c_00 -= pD[jj+ldd*kk] * pD[jj+ldd*kk];
			}
		if(c_00>0)
			{
			f_00_inv = 1.0/sqrt(c_00);
			}
		else
			{
			f_00_inv = 0.0;
			}
		dD[jj] = f_00_inv;
		pD[jj+ldd*jj] = c_00 * f_00_inv;
		// solve lower
		for(ii=jj+1; ii<m; ii++)
			{
			c_00 = pC[ii+ldc*jj];
			for(kk=0; kk<k; kk++)
				{
				c_00 += pA[ii+lda*kk] * pB[jj+ldb*kk];
				}
			for(kk=0; kk<jj; kk++)
				{
				c_00 -= pD[ii+ldd*kk] * pD[jj+ldd*kk];
				}
			pD[ii+ldd*jj] = c_00 * f_00_inv;
			}
		}
	return;
	}



// dgetrf without pivoting
void dgetf2_nopivot(int m, int n, double *A, int lda, double *dA)
	{
	int ii, jj, kk, itmp0, itmp1;
	int iimax = m<n ? m : n;
	int i1 = 1;
	double dtmp;
	double dm1 = -1.0;

	for(ii=0; ii<iimax; ii++)
		{
		itmp0 = m-ii-1;
		dtmp = 1.0/A[ii+lda*ii];
		dA[ii] = dtmp;
		for(jj=0; jj<itmp0; jj++)
			{
			A[ii+1+jj+lda*ii] *= dtmp;
			}
		itmp1 = n-ii-1;
		for(jj=0; jj<itmp1; jj++)
			{
			for(kk=0; kk<itmp0; kk++)
				{
				A[(ii+1+kk)+lda*(ii+1+jj)] -= A[(ii+1+kk)+lda*ii] * A[ii+lda*(ii+1+jj)];
				}
			}
		}
	return;
	}



// dgetrf without pivoting
void dgetrf_nopivot_libstr(int m, int n, struct d_strmat *sC, int ci, int cj, struct d_strmat *sD, int di, int dj)
	{
	int ii, jj, kk;
//	int i1 = 1;
//	double d1 = 1.0;
	double
		d_00_inv, d_11_inv,
		d_00, d_01,
		d_10, d_11;
	int ldc = sC->m;
	int ldd = sD->m;
	double *pC = sC->pA + ci + cj*ldc;
	double *pD = sD->pA + di + dj*ldd;
	double *dD = sD->dA;
	if(di==0 & dj==0)
		sD->use_dA = 1;
	else
		sD->use_dA = 0;
#if 1
	jj = 0;
	for(; jj<n-1; jj+=2)
		{
		// upper
		ii = 0;
		for(; ii<jj-1; ii+=2)
			{
			// correct upper
			d_00 = pC[(ii+0)+ldc*(jj+0)];
			d_10 = pC[(ii+1)+ldc*(jj+0)];
			d_01 = pC[(ii+0)+ldc*(jj+1)];
			d_11 = pC[(ii+1)+ldc*(jj+1)];
			for(kk=0; kk<ii; kk++)
				{
				d_00 -= pD[(ii+0)+ldd*kk] * pD[kk+ldd*(jj+0)];
				d_10 -= pD[(ii+1)+ldd*kk] * pD[kk+ldd*(jj+0)];
				d_01 -= pD[(ii+0)+ldd*kk] * pD[kk+ldd*(jj+1)];
				d_11 -= pD[(ii+1)+ldd*kk] * pD[kk+ldd*(jj+1)];
				}
			// solve upper
			d_10 -= pD[(ii+1)+ldd*kk] * d_00;
			d_11 -= pD[(ii+1)+ldd*kk] * d_01;
			pD[(ii+0)+ldd*(jj+0)] = d_00;
			pD[(ii+1)+ldd*(jj+0)] = d_10;
			pD[(ii+0)+ldd*(jj+1)] = d_01;
			pD[(ii+1)+ldd*(jj+1)] = d_11;
			}
		for(; ii<jj; ii++)
			{
			// correct upper
			d_00 = pC[(ii+0)+ldc*(jj+0)];
			d_01 = pC[(ii+0)+ldc*(jj+1)];
			for(kk=0; kk<ii; kk++)
				{
				d_00 -= pD[(ii+0)+ldd*kk] * pD[kk+ldd*(jj+0)];
				d_01 -= pD[(ii+0)+ldd*kk] * pD[kk+ldd*(jj+1)];
				}
			// solve upper
			pD[(ii+0)+ldd*(jj+0)] = d_00;
			pD[(ii+0)+ldd*(jj+1)] = d_01;
			}
		// diagonal
		ii = jj;
		if(ii<m-1)
			{
			// correct diagonal
			d_00 = pC[(ii+0)+ldc*(jj+0)];
			d_10 = pC[(ii+1)+ldc*(jj+0)];
			d_01 = pC[(ii+0)+ldc*(jj+1)];
			d_11 = pC[(ii+1)+ldc*(jj+1)];
			for(kk=0; kk<ii; kk++)
				{
				d_00 -= pD[(ii+0)+ldd*kk] * pD[kk+ldd*(jj+0)];
				d_10 -= pD[(ii+1)+ldd*kk] * pD[kk+ldd*(jj+0)];
				d_01 -= pD[(ii+0)+ldd*kk] * pD[kk+ldd*(jj+1)];
				d_11 -= pD[(ii+1)+ldd*kk] * pD[kk+ldd*(jj+1)];
				}
			// factorize diagonal
			d_00_inv = 1.0/d_00;
			d_10 *= d_00_inv;
			d_11 -= d_10 * d_01;
			d_11_inv = 1.0/d_11;
			pD[(ii+0)+ldd*(jj+0)] = d_00;
			pD[(ii+1)+ldd*(jj+0)] = d_10;
			pD[(ii+0)+ldd*(jj+1)] = d_01;
			pD[(ii+1)+ldd*(jj+1)] = d_11;
			dD[ii+0] = d_00_inv;
			dD[ii+1] = d_11_inv;
			ii += 2;
			}
		else if(ii<m)
			{
			// correct diagonal
			d_00 = pC[(ii+0)+ldc*(jj+0)];
			d_01 = pC[(ii+0)+ldc*(jj+1)];
			for(kk=0; kk<ii; kk++)
				{
				d_00 -= pD[(ii+0)+ldd*kk] * pD[kk+ldd*(jj+0)];
				d_01 -= pD[(ii+0)+ldd*kk] * pD[kk+ldd*(jj+1)];
				}
			// factorize diagonal
			d_00_inv = 1.0/d_00;
			pD[(ii+0)+ldd*(jj+0)] = d_00;
			pD[(ii+0)+ldd*(jj+1)] = d_01;
			dD[ii+0] = d_00_inv;
			ii += 1;
			}
		// lower
		for(; ii<m-1; ii+=2)
			{
			// correct lower
			d_00 = pC[(ii+0)+ldc*(jj+0)];
			d_10 = pC[(ii+1)+ldc*(jj+0)];
			d_01 = pC[(ii+0)+ldc*(jj+1)];
			d_11 = pC[(ii+1)+ldc*(jj+1)];
			for(kk=0; kk<jj; kk++)
				{
				d_00 -= pD[(ii+0)+ldd*kk] * pD[kk+ldd*(jj+0)];
				d_10 -= pD[(ii+1)+ldd*kk] * pD[kk+ldd*(jj+0)];
				d_01 -= pD[(ii+0)+ldd*kk] * pD[kk+ldd*(jj+1)];
				d_11 -= pD[(ii+1)+ldd*kk] * pD[kk+ldd*(jj+1)];
				}
			// solve lower
			d_00 *= d_00_inv;
			d_10 *= d_00_inv;
			d_01 -= d_00 * pD[kk+ldd*(jj+1)];
			d_11 -= d_10 * pD[kk+ldd*(jj+1)];
			d_01 *= d_11_inv;
			d_11 *= d_11_inv;
			pD[(ii+0)+ldd*(jj+0)] = d_00;
			pD[(ii+1)+ldd*(jj+0)] = d_10;
			pD[(ii+0)+ldd*(jj+1)] = d_01;
			pD[(ii+1)+ldd*(jj+1)] = d_11;
			}
		for(; ii<m; ii++)
			{
			// correct lower
			d_00 = pC[(ii+0)+ldc*(jj+0)];
			d_01 = pC[(ii+0)+ldc*(jj+1)];
			for(kk=0; kk<jj; kk++)
				{
				d_00 -= pD[(ii+0)+ldd*kk] * pD[kk+ldd*(jj+0)];
				d_01 -= pD[(ii+0)+ldd*kk] * pD[kk+ldd*(jj+1)];
				}
			// solve lower
			d_00 *= d_00_inv;
			d_01 -= d_00 * pD[kk+ldd*(jj+1)];
			d_01 *= d_11_inv;
			pD[(ii+0)+ldd*(jj+0)] = d_00;
			pD[(ii+0)+ldd*(jj+1)] = d_01;
			}
		}
	for(; jj<n; jj++)
		{
		// upper
		ii = 0;
		for(; ii<jj-1; ii+=2)
			{
			// correct upper
			d_00 = pC[(ii+0)+ldc*jj];
			d_10 = pC[(ii+1)+ldc*jj];
			for(kk=0; kk<ii; kk++)
				{
				d_00 -= pD[(ii+0)+ldd*kk] * pD[kk+ldd*jj];
				d_10 -= pD[(ii+1)+ldd*kk] * pD[kk+ldd*jj];
				}
			// solve upper
			d_10 -= pD[(ii+1)+ldd*kk] * d_00;
			pD[(ii+0)+ldd*jj] = d_00;
			pD[(ii+1)+ldd*jj] = d_10;
			}
		for(; ii<jj; ii++)
			{
			// correct upper
			d_00 = pC[(ii+0)+ldc*jj];
			for(kk=0; kk<ii; kk++)
				{
				d_00 -= pD[(ii+0)+ldd*kk] * pD[kk+ldd*jj];
				}
			// solve upper
			pD[(ii+0)+ldd*jj] = d_00;
			}
		// diagonal
		ii = jj;
		if(ii<m-1)
			{
			// correct diagonal
			d_00 = pC[(ii+0)+ldc*jj];
			d_10 = pC[(ii+1)+ldc*jj];
			for(kk=0; kk<ii; kk++)
				{
				d_00 -= pD[(ii+0)+ldd*kk] * pD[kk+ldd*jj];
				d_10 -= pD[(ii+1)+ldd*kk] * pD[kk+ldd*jj];
				}
			// factorize diagonal
			d_00_inv = 1.0/d_00;
			d_10 *= d_00_inv;
			pD[(ii+0)+ldd*jj] = d_00;
			pD[(ii+1)+ldd*jj] = d_10;
			dD[ii+0] = d_00_inv;
			ii += 2;
			}
		else if(ii<m)
			{
			// correct diagonal
			d_00 = pC[(ii+0)+ldc*jj];
			for(kk=0; kk<ii; kk++)
				{
				d_00 -= pD[(ii+0)+ldd*kk] * pD[kk+ldd*jj];
				}
			// factorize diagonal
			d_00_inv = 1.0/d_00;
			pD[(ii+0)+ldd*jj] = d_00;
			dD[ii+0] = d_00_inv;
			ii += 1;
			}
		// lower
		for(; ii<m-1; ii+=2)
			{
			// correct lower
			d_00 = pC[(ii+0)+ldc*jj];
			d_10 = pC[(ii+1)+ldc*jj];
			for(kk=0; kk<jj; kk++)
				{
				d_00 -= pD[(ii+0)+ldd*kk] * pD[kk+ldd*jj];
				d_10 -= pD[(ii+1)+ldd*kk] * pD[kk+ldd*jj];
				}
			// solve lower
			d_00 *= d_00_inv;
			d_10 *= d_00_inv;
			pD[(ii+0)+ldd*jj] = d_00;
			pD[(ii+1)+ldd*jj] = d_10;
			}
		for(; ii<m; ii++)
			{
			// correct lower
			d_00 = pC[(ii+0)+ldc*jj];
			for(kk=0; kk<jj; kk++)
				{
				d_00 -= pD[(ii+0)+ldd*kk] * pD[kk+ldd*jj];
				}
			// solve lower
			d_00 *= d_00_inv;
			pD[(ii+0)+ldd*jj] = d_00;
			}
		}
#else
	if(pC!=pD)
		{
		for(jj=0; jj<n; jj++)
			{
			for(ii=0; ii<m; ii++)
				{
				pD[ii+ldd*jj] = pC[ii+ldc*jj];
				}
			}
		}
	dgetf2_nopivot(m, n, pD, ldd, dD);
#endif
	return;
	}



// dgetrf pivoting
void dgetrf_libstr(int m, int n, struct d_strmat *sC, int ci, int cj, struct d_strmat *sD, int di, int dj, int *ipiv)
	{
	int ii, i0, jj, kk, ip, itmp0, itmp1;
	double dtmp, dmax;
	double
		d_00_inv, d_11_inv,
		d_00, d_01,
		d_10, d_11;
	int i1 = 1;
	double d1 = 1.0;
	int ldc = sC->m;
	int ldd = sD->m;
	double *pC = sC->pA+ci+cj*ldc;
	double *pD = sD->pA+di+dj*ldd;
	double *dD = sD->dA;
	if(di==0 & dj==0)
		sD->use_dA = 1;
	else
		sD->use_dA = 0;
	// copy if needed
	if(pC!=pD)
		{
		for(jj=0; jj<n; jj++)
			{
			for(ii=0; ii<m; ii++)
				{
				pD[ii+ldd*jj] = pC[ii+ldc*jj];
				}
			}
		}
	// factorize
#if 1
	jj = 0;
	for(; jj<n-1; jj+=2)
		{
		ii = 0;
		for(; ii<jj-1; ii+=2)
			{
			// correct upper
			d_00 = pD[(ii+0)+ldd*(jj+0)];
			d_10 = pD[(ii+1)+ldd*(jj+0)];
			d_01 = pD[(ii+0)+ldd*(jj+1)];
			d_11 = pD[(ii+1)+ldd*(jj+1)];
			for(kk=0; kk<ii; kk++)
				{
				d_00 -= pD[(ii+0)+ldd*kk] * pD[kk+ldd*(jj+0)];
				d_10 -= pD[(ii+1)+ldd*kk] * pD[kk+ldd*(jj+0)];
				d_01 -= pD[(ii+0)+ldd*kk] * pD[kk+ldd*(jj+1)];
				d_11 -= pD[(ii+1)+ldd*kk] * pD[kk+ldd*(jj+1)];
				}
			// solve upper
			d_10 -= pD[(ii+1)+ldd*kk] * d_00;
			d_11 -= pD[(ii+1)+ldd*kk] * d_01;
			pD[(ii+0)+ldd*(jj+0)] = d_00;
			pD[(ii+1)+ldd*(jj+0)] = d_10;
			pD[(ii+0)+ldd*(jj+1)] = d_01;
			pD[(ii+1)+ldd*(jj+1)] = d_11;
			}
		for(; ii<jj; ii++)
			{
			// correct upper
			d_00 = pD[(ii+0)+ldd*(jj+0)];
			d_01 = pD[(ii+0)+ldd*(jj+1)];
			for(kk=0; kk<ii; kk++)
				{
				d_00 -= pD[(ii+0)+ldd*kk] * pD[kk+ldd*(jj+0)];
				d_01 -= pD[(ii+0)+ldd*kk] * pD[kk+ldd*(jj+1)];
				}
			// solve upper
			pD[(ii+0)+ldd*(jj+0)] = d_00;
			pD[(ii+0)+ldd*(jj+1)] = d_01;
			}
		// correct diagonal and lower and look for pivot
		// correct
		ii = jj;
		i0 = ii;
		for(; ii<m-1; ii+=2)
			{
			d_00 = pD[(ii+0)+ldd*(jj+0)];
			d_10 = pD[(ii+1)+ldd*(jj+0)];
			d_01 = pD[(ii+0)+ldd*(jj+1)];
			d_11 = pD[(ii+1)+ldd*(jj+1)];
			for(kk=0; kk<jj; kk++)
				{
				d_00 -= pD[(ii+0)+ldd*kk] * pD[kk+ldd*(jj+0)];
				d_10 -= pD[(ii+1)+ldd*kk] * pD[kk+ldd*(jj+0)];
				d_01 -= pD[(ii+0)+ldd*kk] * pD[kk+ldd*(jj+1)];
				d_11 -= pD[(ii+1)+ldd*kk] * pD[kk+ldd*(jj+1)];
				}
			pD[(ii+0)+ldd*(jj+0)] = d_00;
			pD[(ii+1)+ldd*(jj+0)] = d_10;
			pD[(ii+0)+ldd*(jj+1)] = d_01;
			pD[(ii+1)+ldd*(jj+1)] = d_11;
			}
		for(; ii<m; ii++)
			{
			d_00 = pD[(ii+0)+ldd*(jj+0)];
			d_01 = pD[(ii+0)+ldd*(jj+1)];
			for(kk=0; kk<jj; kk++)
				{
				d_00 -= pD[(ii+0)+ldd*kk] * pD[kk+ldd*(jj+0)];
				d_01 -= pD[(ii+0)+ldd*kk] * pD[kk+ldd*(jj+1)];
				}
			pD[(ii+0)+ldd*(jj+0)] = d_00;
			pD[(ii+0)+ldd*(jj+1)] = d_01;
			}
		// look for pivot & solve
		// left column
		ii = i0;
		dmax = 0;
		ip = ii;
		for(; ii<m-1; ii+=2)
			{
			d_00 = pD[(ii+0)+ldd*jj];
			d_10 = pD[(ii+1)+ldd*jj];
			dtmp = d_00>0 ? d_00 : -d_00;
			if(dtmp>dmax)
				{
				dmax = dtmp;
				ip = ii+0;
				}
			dtmp = d_10>0 ? d_10 : -d_10;
			if(dtmp>dmax)
				{
				dmax = dtmp;
				ip = ii+1;
				}
			}
		for(; ii<m; ii++)
			{
			d_00 = pD[(ii+0)+ldd*jj];
			dtmp = d_00>0 ? d_00 : -d_00;
			if(dtmp>dmax)
				{
				dmax = dtmp;
				ip = ii+0;
				}
			}
		// row swap
		ii = i0;
		ipiv[ii] = ip;
		if(ip!=ii)
			{
			for(kk=0; kk<n; kk++)
				{
				dtmp = pD[ii+ldd*kk];
				pD[ii+ldd*kk] = pD[ip+ldd*kk];
				pD[ip+ldd*kk] = dtmp;
				}
			}
		// factorize diagonal
		d_00 = pD[(ii+0)+ldd*(jj+0)];
		d_00_inv = 1.0/d_00;
		pD[(ii+0)+ldd*(jj+0)] = d_00;
		dD[ii] = d_00_inv;
		ii += 1;
		// solve & compute next pivot
		dmax = 0;
		ip = ii;
		for(; ii<m-1; ii+=2)
			{
			d_00 = pD[(ii+0)+ldd*(jj+0)];
			d_10 = pD[(ii+1)+ldd*(jj+0)];
			d_00 *= d_00_inv;
			d_10 *= d_00_inv;
			d_01 = pD[(ii+0)+ldd*(jj+1)];
			d_11 = pD[(ii+1)+ldd*(jj+1)];
			d_01 -= d_00 * pD[jj+ldd*(jj+1)];
			d_11 -= d_10 * pD[jj+ldd*(jj+1)];
			dtmp = d_01>0 ? d_01 : -d_01;
			if(dtmp>dmax)
				{
				dmax = dtmp;
				ip = ii+0;
				}
			dtmp = d_11>0 ? d_11 : -d_11;
			if(dtmp>dmax)
				{
				dmax = dtmp;
				ip = ii+1;
				}
			pD[(ii+0)+ldd*(jj+0)] = d_00;
			pD[(ii+1)+ldd*(jj+0)] = d_10;
			pD[(ii+0)+ldd*(jj+1)] = d_01;
			pD[(ii+1)+ldd*(jj+1)] = d_11;
			}
		for(; ii<m; ii++)
			{
			d_00 = pD[(ii+0)+ldd*(jj+0)];
			d_00 *= d_00_inv;
			d_01 = pD[(ii+0)+ldd*(jj+1)];
			d_01 -= d_00 * pD[jj+ldd*(jj+1)];
			dtmp = d_01>0 ? d_01 : -d_01;
			if(dtmp>dmax)
				{
				dmax = dtmp;
				ip = ii+0;
				}
			pD[(ii+0)+ldd*(jj+0)] = d_00;
			pD[(ii+0)+ldd*(jj+1)] = d_01;
			}
		// row swap
		ii = i0+1;
		ipiv[ii] = ip;
		if(ip!=ii)
			{
			for(kk=0; kk<n; kk++)
				{
				dtmp = pD[ii+ldd*kk];
				pD[ii+ldd*kk] = pD[ip+ldd*kk];
				pD[ip+ldd*kk] = dtmp;
				}
			}
		// factorize diagonal
		d_00 = pD[ii+ldd*(jj+1)];
		d_00_inv = 1.0/d_00;
		pD[ii+ldd*(jj+1)] = d_00;
		dD[ii] = d_00_inv;
		ii += 1;
		// solve lower
		for(; ii<m; ii++)
			{
			d_00 = pD[ii+ldd*(jj+1)];
			d_00 *= d_00_inv;
			pD[ii+ldd*(jj+1)] = d_00;
			}
		}
	for(; jj<n; jj++)
		{
		ii = 0;
		for(; ii<jj-1; ii+=2)
			{
			// correct upper
			d_00 = pD[(ii+0)+ldd*jj];
			d_10 = pD[(ii+1)+ldd*jj];
			for(kk=0; kk<ii; kk++)
				{
				d_00 -= pD[(ii+0)+ldd*kk] * pD[kk+ldd*jj];
				d_10 -= pD[(ii+1)+ldd*kk] * pD[kk+ldd*jj];
				}
			// solve upper
			d_10 -= pD[(ii+1)+ldd*kk] * d_00;
			pD[(ii+0)+ldd*jj] = d_00;
			pD[(ii+1)+ldd*jj] = d_10;
			}
		for(; ii<jj; ii++)
			{
			// correct upper
			d_00 = pD[ii+ldd*jj];
			for(kk=0; kk<ii; kk++)
				{
				d_00 -= pD[ii+ldd*kk] * pD[kk+ldd*jj];
				}
			// solve upper
			pD[ii+ldd*jj] = d_00;
			}
		i0 = ii;
		ii = jj;
		// correct diagonal and lower and look for pivot
		dmax = 0;
		ip = ii;
		for(; ii<m-1; ii+=2)
			{
			d_00 = pD[(ii+0)+ldd*jj];
			d_10 = pD[(ii+1)+ldd*jj];
			for(kk=0; kk<jj; kk++)
				{
				d_00 -= pD[(ii+0)+ldd*kk] * pD[kk+ldd*jj];
				d_10 -= pD[(ii+1)+ldd*kk] * pD[kk+ldd*jj];
				}
			dtmp = d_00>0 ? d_00 : -d_00;
			if(dtmp>dmax)
				{
				dmax = dtmp;
				ip = ii+0;
				}
			dtmp = d_10>0 ? d_10 : -d_10;
			if(dtmp>dmax)
				{
				dmax = dtmp;
				ip = ii+1;
				}
			pD[(ii+0)+ldd*jj] = d_00;
			pD[(ii+1)+ldd*jj] = d_10;
			}
		for(; ii<m; ii++)
			{
			d_00 = pD[(ii+0)+ldd*jj];
			for(kk=0; kk<jj; kk++)
				{
				d_00 -= pD[(ii+0)+ldd*kk] * pD[kk+ldd*jj];
				}
			dtmp = d_00>0 ? d_00 : -d_00;
			if(dtmp>dmax)
				{
				dmax = dtmp;
				ip = ii+0;
				}
			pD[(ii+0)+ldd*jj] = d_00;
			}
		// row swap
		ii = i0;
		ipiv[ii] = ip;
		if(ip!=ii)
			{
			for(kk=0; kk<n; kk++)
				{
				dtmp = pD[ii+ldd*kk];
				pD[ii+ldd*kk] = pD[ip+ldd*kk];
				pD[ip+ldd*kk] = dtmp;
				}
			}
		// factorize diagonal
		d_00 = pD[ii+ldd*jj];
		d_00_inv = 1.0/d_00;
		pD[ii+ldd*jj] = d_00;
		dD[ii] = d_00_inv;
		ii += 1;
		for(; ii<m; ii++)
			{
			// correct lower
			d_00 = pD[ii+ldd*jj];
			// solve lower
			d_00 *= d_00_inv;
			pD[ii+ldd*jj] = d_00;
			}
		}
#else
	int iimax = m<n ? m : n;
	for(ii=0; ii<iimax; ii++)
		{
		dmax = (pD[ii+ldd*ii]>0 ? pD[ii+ldd*ii] : -pD[ii+ldd*ii]);
		ip = ii;
		for(jj=1; jj<m-ii; jj++)
			{
			dtmp = pD[ii+jj+ldd*ii]>0 ? pD[ii+jj+ldd*ii] : -pD[ii+jj+ldd*ii];
			if(dtmp>dmax)
				{
				dmax = dtmp;
				ip = ii+jj;
				}
			}
		ipiv[ii] = ip;
		if(ip!=ii)
			{
			for(jj=0; jj<n; jj++)
				{
				dtmp = pD[ii+jj*ldd];
				pD[ii+jj*ldd] = pD[ip+jj*ldd];
				pD[ip+jj*ldd] = dtmp;
				}
			}
		itmp0 = m-ii-1;
		dtmp = 1.0/pD[ii+ldd*ii];
		dD[ii] = dtmp;
		for(jj=0; jj<itmp0; jj++)
			{
			pD[ii+1+jj+ldd*ii] *= dtmp;
			}
		itmp1 = n-ii-1;
		for(jj=0; jj<itmp1; jj++)
			{
			for(kk=0; kk<itmp0; kk++)
				{
				pD[(ii+1+kk)+ldd*(ii+1+jj)] -= pD[(ii+1+kk)+ldd*ii] * pD[ii+ldd*(ii+1+jj)];
				}
			}
		}
#endif
	return;
	}



#elif defined(LA_BLAS)



// dpotrf
void dpotrf_l_libstr(int m, int n, struct d_strmat *sC, int ci, int cj, struct d_strmat *sD, int di, int dj)
	{
	int jj;
	char cl = 'l';
	char cn = 'n';
	char cr = 'r';
	char ct = 't';
	int i1 = 1;
	int mmn = m-n;
	int info;
	double d1 = 1.0;
	double *pC = sC->pA+ci+cj*sC->m;
	double *pD = sD->pA+di+dj*sD->m;
	if(!(pC==pD))
		{
		for(jj=0; jj<n; jj++)
#if defined(LA_BLAS_MKL)
			dcopy(&m, pC+jj*sC->m, &i1, pD+jj*sD->m, &i1);
#else
			dcopy_(&m, pC+jj*sC->m, &i1, pD+jj*sD->m, &i1);
#endif
		}
#if defined(LA_BLAS_MKL)
	dpotrf(&cl, &n, pD, &(sD->m), &info);
	dtrsm(&cr, &cl, &ct, &cn, &mmn, &n, &d1, pD, &(sD->m), pD+n, &(sD->m));
#else
	dpotrf_(&cl, &n, pD, &(sD->m), &info);
	dtrsm_(&cr, &cl, &ct, &cn, &mmn, &n, &d1, pD, &(sD->m), pD+n, &(sD->m));
#endif
	return;
	}



// dgetrf without pivoting
void dgetf2_nopivot(int m, int n, double *A, int lda)
	{
	if(m<=0 | n<=0)
		return;
	int i, j, itmp0, itmp1;
	int jmax = m<n ? m : n;
	int i1 = 1;
	double dtmp;
	double dm1 = -1.0;

	for(j=0; j<jmax; j++)
		{
		itmp0 = m-j-1;
		dtmp = 1.0/A[j+lda*j];
#if defined(LA_BLAS_MKL)
		dscal(&itmp0, &dtmp, &A[(j+1)+lda*j], &i1);
#else
		dscal_(&itmp0, &dtmp, &A[(j+1)+lda*j], &i1);
#endif
		itmp1 = n-j-1;
#if defined(LA_BLAS_MKL)
		dger(&itmp0, &itmp1, &dm1, &A[(j+1)+lda*j], &i1, &A[j+lda*(j+1)], &lda, &A[(j+1)+lda*(j+1)], &lda);
#else
		dger_(&itmp0, &itmp1, &dm1, &A[(j+1)+lda*j], &i1, &A[j+lda*(j+1)], &lda, &A[(j+1)+lda*(j+1)], &lda);
#endif
		}

	return;

	}

// dsyrk dpotrf
void dsyrk_dpotrf_ln_libstr(int m, int n, int k, struct d_strmat *sA, int ai, int aj, struct d_strmat *sB, int bi, int bj, struct d_strmat *sC, int ci, int cj, struct d_strmat *sD, int di, int dj)
	{
	int jj;
	char cl = 'l';
	char cn = 'n';
	char cr = 'r';
	char ct = 't';
	char cu = 'u';
	int i1 = 1;
	double d1 = 1.0;
	int mmn = m-n;
	int info;
	double *pA = sA->pA+ai+aj*sA->m;
	double *pB = sB->pA+bi+bj*sB->m;
	double *pC = sC->pA+ci+cj*sC->m;
	double *pD = sD->pA+di+dj*sD->m;
	if(!(pC==pD))
		{
		for(jj=0; jj<n; jj++)
#if defined(LA_BLAS_MKL)
			dcopy(&m, pC+jj*sC->m, &i1, pD+jj*sD->m, &i1);
#else
			dcopy_(&m, pC+jj*sC->m, &i1, pD+jj*sD->m, &i1);
#endif
		}
#if defined(LA_BLAS_MKL)
	dgemm(&cn, &ct, &m, &n, &k, &d1, pA, &(sA->m), pB, &(sB->m), &d1, pD, &(sD->m));
	dpotrf(&cl, &n, pD, &(sD->m), &info);
	dtrsm(&cr, &cl, &ct, &cn, &mmn, &n, &d1, pD, &(sD->m), pD+n, &(sD->m));
#else
	dgemm_(&cn, &ct, &m, &n, &k, &d1, pA, &(sA->m), pB, &(sB->m), &d1, pD, &(sD->m));
	dpotrf_(&cl, &n, pD, &(sD->m), &info);
	dtrsm_(&cr, &cl, &ct, &cn, &mmn, &n, &d1, pD, &(sD->m), pD+n, &(sD->m));
#endif
	return;
	}



// dgetrf without pivoting
void dgetrf_nopivot_libstr(int m, int n, struct d_strmat *sC, int ci, int cj, struct d_strmat *sD, int di, int dj)
	{
	// TODO with custom level 2 LAPACK + level 3 BLAS
//	printf("\nfeature not implemented yet\n\n");
//	exit(1);
	int jj;
	int i1 = 1;
	double d1 = 1.0;
	double *pC = sC->pA+ci+cj*sC->m;
	double *pD = sD->pA+di+dj*sD->m;
	if(!(pC==pD))
		{
		for(jj=0; jj<n; jj++)
#if defined(LA_BLAS_MKL)
			dcopy(&m, pC+jj*sC->m, &i1, pD+jj*sD->m, &i1);
#else
			dcopy_(&m, pC+jj*sC->m, &i1, pD+jj*sD->m, &i1);
#endif
		}
	dgetf2_nopivot(m, n, pD, sD->m);
	return;
	}



// dgetrf pivoting
void dgetrf_libstr(int m, int n, struct d_strmat *sC, int ci, int cj, struct d_strmat *sD, int di, int dj, int *ipiv)
	{
	// TODO with custom level 2 LAPACK + level 3 BLAS
//	printf("\nfeature not implemented yet\n\n");
//	exit(1);
	int jj;
	int i1 = 1;
	double d1 = 1.0;
	double *pC = sC->pA+ci+cj*sC->m;
	double *pD = sD->pA+di+dj*sD->m;
	if(!(pC==pD))
		{
		for(jj=0; jj<n; jj++)
#if defined(LA_BLAS_MKL)
			dcopy(&m, pC+jj*sC->m, &i1, pD+jj*sD->m, &i1);
#else
			dcopy_(&m, pC+jj*sC->m, &i1, pD+jj*sD->m, &i1);
#endif
		}
	int info;
#if defined(LA_BLAS_MKL)
	dgetrf(&m, &n, pD, &(sD->m), ipiv, &info);
#else
	dgetrf_(&m, &n, pD, &(sD->m), ipiv, &info);
#endif
	int tmp = m<n ? m : n;
	for(jj=0; jj<tmp; jj++)
		ipiv[jj] -= 1;
	return;
	}



#else

#error : wrong LA choice

#endif
