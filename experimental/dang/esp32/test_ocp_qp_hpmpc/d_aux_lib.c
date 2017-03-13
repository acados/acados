/**************************************************************************************************
* acados/external/blasfeo/auxiliary/d_aux_lib.c                                                                                                *
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

#include "blasfeo_common.h"

#if defined(LA_REFERENCE) | defined(LA_BLAS)



// return memory size (in bytes) needed for a strmat
int d_size_strmat(int m, int n)
	{
#if defined(LA_REFERENCE)
	int tmp = m<n ? m : n; // al(min(m,n)) // XXX max ???
	int size = (m*n+tmp)*sizeof(double);
#else
	int size = (m*n)*sizeof(double);
#endif
	return size;
	}



// return memory size (in bytes) needed for the diagonal of a strmat
int d_size_diag_strmat(int m, int n)
	{
	int size = 0;
#if defined(LA_REFERENCE)
	int tmp = m<n ? m : n; // al(min(m,n)) // XXX max ???
	size = tmp*sizeof(double);
#endif
	return size;
	}



// create a matrix structure for a matrix of size m*n by using memory passed by a pointer
void d_create_strmat(int m, int n, struct d_strmat *sA, void *memory)
	{
	sA->m = m;
	sA->n = n;
	double *ptr = (double *) memory;
	sA->pA = ptr;
	ptr += m*n;
#if defined(LA_REFERENCE)
	int tmp = m<n ? m : n; // al(min(m,n)) // XXX max ???
	sA->dA = ptr;
	ptr += tmp;
	sA->use_dA = 0;
	sA->memory_size = (m*n+tmp)*sizeof(double);
#else
	sA->memory_size = (m*n)*sizeof(double);
#endif
	return;
	}



// return memory size (in bytes) needed for a strvec
int d_size_strvec(int m)
	{
	int size = m*sizeof(double);
	return size;
	}



// create a matrix structure for a matrix of size m*n by using memory passed by a pointer
void d_create_strvec(int m, struct d_strvec *sa, void *memory)
	{
	sa->m = m;
	double *ptr = (double *) memory;
	sa->pa = ptr;
//	ptr += m * n;
	sa->memory_size = m*sizeof(double);
	return;
	}



// convert a matrix into a matrix structure
void d_cvt_mat2strmat(int m, int n, double *A, int lda, struct d_strmat *sA, int ai, int aj)
	{
	int ii, jj;
	int lda2 = sA->m;
	double *pA = sA->pA + ai + aj*lda2;
	for(jj=0; jj<n; jj++)
		{
		ii = 0;
		for(; ii<m-3; ii+=4)
			{
			pA[ii+0+jj*lda2] = A[ii+0+jj*lda];
			pA[ii+1+jj*lda2] = A[ii+1+jj*lda];
			pA[ii+2+jj*lda2] = A[ii+2+jj*lda];
			pA[ii+3+jj*lda2] = A[ii+3+jj*lda];
			}
		for(; ii<m; ii++)
			{
			pA[ii+0+jj*lda2] = A[ii+0+jj*lda];
			}
		}
	return;
	}



// convert and transpose a matrix into a matrix structure
void d_cvt_tran_mat2strmat(int m, int n, double *A, int lda, struct d_strmat *sA, int ai, int aj)
	{
	int ii, jj;
	int lda2 = sA->m;
	double *pA = sA->pA + ai + aj*lda2;
	for(jj=0; jj<n; jj++)
		{
		ii = 0;
		for(; ii<m-3; ii+=4)
			{
			pA[jj+(ii+0)*lda2] = A[ii+0+jj*lda];
			pA[jj+(ii+1)*lda2] = A[ii+1+jj*lda];
			pA[jj+(ii+2)*lda2] = A[ii+2+jj*lda];
			pA[jj+(ii+3)*lda2] = A[ii+3+jj*lda];
			}
		for(; ii<m; ii++)
			{
			pA[jj+(ii+0)*lda2] = A[ii+0+jj*lda];
			}
		}
	return;
	}



// convert a vector into a vector structure
void d_cvt_vec2strvec(int m, double *a, struct d_strvec *sa, int ai)
	{
	double *pa = sa->pa + ai;
	int ii;
	for(ii=0; ii<m; ii++)
		pa[ii] = a[ii];
	return;
	}



// convert a matrix structure into a matrix
void d_cvt_strmat2mat(int m, int n, struct d_strmat *sA, int ai, int aj, double *A, int lda)
	{
	int ii, jj;
	int lda2 = sA->m;
	double *pA = sA->pA + ai + aj*lda2;
	for(jj=0; jj<n; jj++)
		{
		ii = 0;
		for(; ii<m-3; ii+=4)
			{
			A[ii+0+jj*lda] = pA[ii+0+jj*lda2];
			A[ii+1+jj*lda] = pA[ii+1+jj*lda2];
			A[ii+2+jj*lda] = pA[ii+2+jj*lda2];
			A[ii+3+jj*lda] = pA[ii+3+jj*lda2];
			}
		for(; ii<m; ii++)
			{
			A[ii+0+jj*lda] = pA[ii+0+jj*lda2];
			}
		}
	return;
	}



// convert and transpose a matrix structure into a matrix
void d_cvt_tran_strmat2mat(int m, int n, struct d_strmat *sA, int ai, int aj, double *A, int lda)
	{
	int ii, jj;
	int lda2 = sA->m;
	double *pA = sA->pA + ai + aj*lda2;
	for(jj=0; jj<n; jj++)
		{
		ii = 0;
		for(; ii<m-3; ii+=4)
			{
			A[jj+(ii+0)*lda] = pA[ii+0+jj*lda2];
			A[jj+(ii+1)*lda] = pA[ii+1+jj*lda2];
			A[jj+(ii+2)*lda] = pA[ii+2+jj*lda2];
			A[jj+(ii+3)*lda] = pA[ii+3+jj*lda2];
			}
		for(; ii<m; ii++)
			{
			A[jj+(ii+0)*lda] = pA[ii+0+jj*lda2];
			}
		}
	return;
	}



// convert a vector structure into a vector
void d_cvt_strvec2vec(int m, struct d_strvec *sa, int ai, double *a)
	{
	double *pa = sa->pa + ai;
	int ii;
	for(ii=0; ii<m; ii++)
		a[ii] = pa[ii];
	return;
	}



// cast a matrix into a matrix structure
void d_cast_mat2strmat(double *A, struct d_strmat *sA)
	{
	sA->pA = A;
	return;
	}



// cast a matrix into the diagonal of a matrix structure
void d_cast_diag_mat2strmat(double *dA, struct d_strmat *sA)
	{
#if defined(LA_REFERENCE)
	sA->dA = dA;
#endif
	return;
	}



// cast a vector into a vector structure
void d_cast_vec2vecmat(double *a, struct d_strvec *sa)
	{
	sa->pa = a;
	return;
	}



// insert element into strmat
void dmatin1_libstr(double a, struct d_strmat *sA, int ai, int aj)
	{
	int lda = sA->m;
	double *pA = sA->pA + ai + aj*lda;
	pA[0] = a;
	return;
	}



// extract element from strmat
double dmatex1_libstr(struct d_strmat *sA, int ai, int aj)
	{
	int lda = sA->m;
	double *pA = sA->pA + ai + aj*lda;
	return pA[0];
	}



// insert element into strvec
void dvecin1_libstr(double a, struct d_strvec *sx, int xi)
	{
	double *x = sx->pa + xi;
	x[0] = a;
	return;
	}



// extract element from strvec
double dvecex1_libstr(struct d_strvec *sx, int xi)
	{
	double *x = sx->pa + xi;
	return x[0];
	}



// set all elements of a strmat to a value
void dmatse_libstr(int m, int n, double alpha, struct d_strmat *sA, int ai, int aj)
	{
	int lda = sA->m;
	double *pA = sA->pA + ai + aj*lda;
	int ii, jj;
	for(jj=0; jj<n; jj++)
		{
		for(ii=0; ii<m; ii++)
			{
			pA[ii+lda*jj] = alpha;
			}
		}
	return;
	}



// set all elements of a strvec to a value
void dvecse_libstr(int m, double alpha, struct d_strvec *sx, int xi)
	{
	double *x = sx->pa + xi;
	int ii;
	for(ii=0; ii<m; ii++)
		x[ii] = alpha;
	return;
	}



// insert a vector into diagonal
void ddiain_libstr(int kmax, double alpha, struct d_strvec *sx, int xi, struct d_strmat *sA, int ai, int aj)
	{
	int lda = sA->m;
	double *pA = sA->pA + ai + aj*lda;
	double *x = sx->pa + xi;
	int ii;
	for(ii=0; ii<kmax; ii++)
		pA[ii*(lda+1)] = alpha*x[ii];
	return;
	}



// extract a row into a vector
void drowex_libstr(int kmax, double alpha, struct d_strmat *sA, int ai, int aj, struct d_strvec *sx, int xi)
	{
	int lda = sA->m;
	double *pA = sA->pA + ai + aj*lda;
	double *x = sx->pa + xi;
	int ii;
	for(ii=0; ii<kmax; ii++)
		x[ii] = alpha*pA[ii*lda];
	return;
	}



// insert a vector  into a row
void drowin_libstr(int kmax, double alpha, struct d_strvec *sx, int xi, struct d_strmat *sA, int ai, int aj)
	{
	int lda = sA->m;
	double *pA = sA->pA + ai + aj*lda;
	double *x = sx->pa + xi;
	int ii;
	for(ii=0; ii<kmax; ii++)
		pA[ii*lda] = alpha*x[ii];
	return;
	}



// add a vector to a row
void drowad_libstr(int kmax, double alpha, struct d_strvec *sx, int xi, struct d_strmat *sA, int ai, int aj)
	{
	int lda = sA->m;
	double *pA = sA->pA + ai + aj*lda;
	double *x = sx->pa + xi;
	int ii;
	for(ii=0; ii<kmax; ii++)
		pA[ii*lda] += alpha*x[ii];
	return;
	}



// swap two rows of a matrix struct
void drowsw_libstr(int kmax, struct d_strmat *sA, int ai, int aj, struct d_strmat *sC, int ci, int cj)
	{
	int lda = sA->m;
	double *pA = sA->pA + ai + aj*lda;
	int ldc = sC->m;
	double *pC = sC->pA + ci + cj*lda;
	int ii;
	double tmp;
	for(ii=0; ii<kmax; ii++)
		{
		tmp = pA[ii*lda];
		pA[ii*lda] = pC[ii*ldc];
		pC[ii*ldc] = tmp;
		}
	return;
	}



// permute the rows of a matrix struct
void drowpe_libstr(int kmax, int *ipiv, struct d_strmat *sA)
	{
	int ii;
	for(ii=0; ii<kmax; ii++)
		{
		if(ipiv[ii]!=ii)
			drowsw_libstr(sA->n, sA, ii, 0, sA, ipiv[ii], 0);
		}
	return;
	}



// swap two cols of a matrix struct
void dcolsw_libstr(int kmax, struct d_strmat *sA, int ai, int aj, struct d_strmat *sC, int ci, int cj)
	{
	int lda = sA->m;
	double *pA = sA->pA + ai + aj*lda;
	int ldc = sC->m;
	double *pC = sC->pA + ci + cj*lda;
	int ii;
	double tmp;
	for(ii=0; ii<kmax; ii++)
		{
		tmp = pA[ii];
		pA[ii] = pC[ii];
		pC[ii] = tmp;
		}
	return;
	}



// permute the cols of a matrix struct
void dcolpe_libstr(int kmax, int *ipiv, struct d_strmat *sA)
	{
	int ii;
	for(ii=0; ii<kmax; ii++)
		{
		if(ipiv[ii]!=ii)
			dcolsw_libstr(sA->m, sA, 0, ii, sA, 0, ipiv[ii]);
		}
	return;
	}



// copy a generic strmat into a generic strmat
void dgecp_libstr(int m, int n, double alpha, struct d_strmat *sA, int ai, int aj, struct d_strmat *sC, int ci, int cj)
	{
	int lda = sA->m;
	double *pA = sA->pA + ai + aj*lda;
	int ldc = sC->m;
	double *pC = sC->pA + ci + cj*ldc;
	int ii, jj;
	for(jj=0; jj<n; jj++)
		{
		ii = 0;
		for(; ii<m-3; ii+=4)
			{
			pC[ii+0+jj*ldc] = alpha*pA[ii+0+jj*lda];
			pC[ii+1+jj*ldc] = alpha*pA[ii+1+jj*lda];
			pC[ii+2+jj*ldc] = alpha*pA[ii+2+jj*lda];
			pC[ii+3+jj*ldc] = alpha*pA[ii+3+jj*lda];
			}
		for(; ii<m; ii++)
			{
			pC[ii+0+jj*ldc] = alpha*pA[ii+0+jj*lda];
			}
		}
	return;
	}



// copy a strvec into a strvec
void dveccp_libstr(int m, double alpha, struct d_strvec *sa, int ai, struct d_strvec *sc, int ci)
	{
	double *pa = sa->pa + ai;
	double *pc = sc->pa + ci;
	int ii;
	ii = 0;
	for(; ii<m-3; ii+=4)
		{
		pc[ii+0] = alpha*pa[ii+0];
		pc[ii+1] = alpha*pa[ii+1];
		pc[ii+2] = alpha*pa[ii+2];
		pc[ii+3] = alpha*pa[ii+3];
		}
	for(; ii<m; ii++)
		{
		pc[ii+0] = alpha*pa[ii+0];
		}
	return;
	}



// copy a lower triangular strmat into a lower triangular strmat
void dtrcp_l_libstr(int m, double alpha, struct d_strmat *sA, int ai, int aj, struct d_strmat *sC, int ci, int cj)
	{
	int lda = sA->m;
	double *pA = sA->pA + ai + aj*lda;
	int ldc = sC->m;
	double *pC = sC->pA + ci + cj*ldc;
	int ii, jj;
	for(jj=0; jj<m; jj++)
		{
		ii = jj;
		for(; ii<m; ii++)
			{
			pC[ii+0+jj*ldc] = alpha*pA[ii+0+jj*lda];
			}
		}
	return;
	}



// scale and add a generic strmat into a generic strmat
void dgead_libstr(int m, int n, double alpha, struct d_strmat *sA, int ai, int aj, struct d_strmat *sC, int ci, int cj)
	{
	int lda = sA->m;
	double *pA = sA->pA + ai + aj*lda;
	int ldc = sC->m;
	double *pC = sC->pA + ci + cj*ldc;
	int ii, jj;
	for(jj=0; jj<n; jj++)
		{
		ii = 0;
		for(; ii<m-3; ii+=4)
			{
			pC[ii+0+jj*ldc] += alpha*pA[ii+0+jj*lda];
			pC[ii+1+jj*ldc] += alpha*pA[ii+1+jj*lda];
			pC[ii+2+jj*ldc] += alpha*pA[ii+2+jj*lda];
			pC[ii+3+jj*ldc] += alpha*pA[ii+3+jj*lda];
			}
		for(; ii<m; ii++)
			{
			pC[ii+0+jj*ldc] += alpha*pA[ii+0+jj*lda];
			}
		}
	return;
	}



// copy and transpose a generic strmat into a generic strmat
void dgetr_libstr(int m, int n, double alpha, struct d_strmat *sA, int ai, int aj, struct d_strmat *sC, int ci, int cj)
	{
	int lda = sA->m;
	double *pA = sA->pA + ai + aj*lda;
	int ldc = sC->m;
	double *pC = sC->pA + ci + cj*ldc;
	int ii, jj;
	for(jj=0; jj<n; jj++)
		{
		ii = 0;
		for(; ii<m-3; ii+=4)
			{
			pC[jj+(ii+0)*ldc] = alpha * pA[ii+0+jj*lda];
			pC[jj+(ii+1)*ldc] = alpha * pA[ii+1+jj*lda];
			pC[jj+(ii+2)*ldc] = alpha * pA[ii+2+jj*lda];
			pC[jj+(ii+3)*ldc] = alpha * pA[ii+3+jj*lda];
			}
		for(; ii<m; ii++)
			{
			pC[jj+(ii+0)*ldc] = alpha * pA[ii+0+jj*lda];
			}
		}
	return;
	}



// copy and transpose a lower triangular strmat into an upper triangular strmat
void dtrtr_l_libstr(int m, double alpha, struct d_strmat *sA, int ai, int aj, struct d_strmat *sC, int ci, int cj)
	{
	int lda = sA->m;
	double *pA = sA->pA + ai + aj*lda;
	int ldc = sC->m;
	double *pC = sC->pA + ci + cj*ldc;
	int ii, jj;
	for(jj=0; jj<m; jj++)
		{
		ii = jj;
		for(; ii<m; ii++)
			{
			pC[jj+(ii+0)*ldc] = alpha * pA[ii+0+jj*lda];
			}
		}
	return;
	}



// copy and transpose an upper triangular strmat into a lower triangular strmat
void dtrtr_u_libstr(int m, double alpha, struct d_strmat *sA, int ai, int aj, struct d_strmat *sC, int ci, int cj)
	{
	int lda = sA->m;
	double *pA = sA->pA + ai + aj*lda;
	int ldc = sC->m;
	double *pC = sC->pA + ci + cj*ldc;
	int ii, jj;
	for(jj=0; jj<m; jj++)
		{
		ii = 0;
		for(; ii<=jj; ii++)
			{
			pC[jj+(ii+0)*ldc] = alpha * pA[ii+0+jj*lda];
			}
		}
	return;
	}



// insert a strvec to the diagonal of a strmat, sparse formulation
void ddiain_libspstr(int kmax, int *idx, double alpha, struct d_strvec *sx, int xi, struct d_strmat *sD, int di, int dj)
	{
	double *x = sx->pa + xi;
	int ldd = sD->m;
	double *pD = sD->pA + di + dj*ldd;
	int ii, jj;
	for(jj=0; jj<kmax; jj++)
		{
		ii = idx[jj];
		pD[ii*(ldd+1)] = alpha * x[jj];
		}
	return;
	}



// extract the diagonal of a strmat from a strvec , sparse formulation
void ddiaex_libspstr(int kmax, int *idx, double alpha, struct d_strmat *sD, int di, int dj, struct d_strvec *sx, int xi)
	{
	double *x = sx->pa + xi;
	int ldd = sD->m;
	double *pD = sD->pA + di + dj*ldd;
	int ii, jj;
	for(jj=0; jj<kmax; jj++)
		{
		ii = idx[jj];
		x[jj] = alpha * pD[ii*(ldd+1)];
		}
	return;
	}



// add scaled strvec to another strvec and insert to diagonal of strmat, sparse formulation
void ddiaad_libspstr(int kmax, int *idx, double alpha, struct d_strvec *sx, int xi, struct d_strmat *sD, int di, int dj)
	{
	double *x = sx->pa + xi;
	int ldd = sD->m;
	double *pD = sD->pA + di + dj*ldd;
	int ii, jj;
	for(jj=0; jj<kmax; jj++)
		{
		ii = idx[jj];
		pD[ii*(ldd+1)] += alpha * x[jj];
		}
	return;
	}



// add scaled strvec to another strvec and insert to diagonal of strmat, sparse formulation
void ddiaadin_libspstr(int kmax, int *idx, double alpha, struct d_strvec *sx, int xi, struct d_strvec *sy, int yi, struct d_strmat *sD, int di, int dj)
	{
	double *x = sx->pa + xi;
	double *y = sy->pa + yi;
	int ldd = sD->m;
	double *pD = sD->pA + di + dj*ldd;
	int ii, jj;
	for(jj=0; jj<kmax; jj++)
		{
		ii = idx[jj];
		pD[ii*(ldd+1)] = y[jj] + alpha * x[jj];
		}
	return;
	}



// add scaled strvec to row of strmat, sparse formulation
void drowad_libspstr(int kmax, int *idx, double alpha, struct d_strvec *sx, int xi, struct d_strmat *sD, int di, int dj)
	{
	double *x = sx->pa + xi;
	int ldd = sD->m;
	double *pD = sD->pA + di + dj*ldd;
	int ii, jj;
	for(jj=0; jj<kmax; jj++)
		{
		ii = idx[jj];
		pD[ii*ldd] += alpha * x[jj];
		}
	return;
	}




// adds strvec to strvec, sparse formulation
void dvecad_libspstr(int kmax, int *idx, double alpha, struct d_strvec *sx, int xi, struct d_strvec *sy, int yi)
	{
	int jj;
	double *x = sx->pa + xi;
	double *y = sy->pa + yi;
//	dvecad_libsp(kmax, idx, alpha, x, y);
	for(jj=0; jj<kmax; jj++)
		{
		y[idx[jj]] += alpha * x[jj];
		}
	return;
	}



#else

#error : wrong LA choice

#endif
