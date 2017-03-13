/**************************************************************************************************
* /home/dang/acados/external/hpmpc/auxiliary/d_aux_extern_depend_lib4.c                                                                                                *
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
#include <string.h>

#include "block_size.h"
#include "kernel_d_lib4.h"

int posix_memalign(void **memptr, size_t alignment, size_t size);

// 2017.03.13 Dang add target.h
#include "target.h"


#if ! defined(BLASFEO)

/* creates a zero matrix aligned */
void d_zeros(double **pA, int row, int col)
	{
	void *temp = malloc((row*col)*sizeof(double));
	*pA = temp;
	double *A = *pA;
	int i;
	for(i=0; i<row*col; i++) A[i] = 0.0;
	}



/* creates a zero matrix aligned to a cache line */
void d_zeros_align(double **pA, int row, int col)
	{
#if defined(OS_WINDOWS)
	*pA = (double *) _aligned_malloc( (row*col)*sizeof(double), 64 );
#else
	void *temp;
	int err = posix_memalign(&temp, 64, (row*col)*sizeof(double));
	if(err!=0)
		{
		printf("\nMemory allocation error\n");
		exit(1);
		}
	*pA = temp;
#endif
//	double *A = *pA;
//	int i;
//	for(i=0; i<row*col; i++) A[i] = 0.0;
	memset(temp, 0, row*col*sizeof(double));
	}



/* allocates memory aligned to a cache line */
void v_zeros_align(void **pA, int size_in_bytes)
	{
#if defined(OS_WINDOWS)
	*pA = _aligned_malloc( size_in_bytes, 64 );
#else
	void *temp;
	int err = posix_memalign(pA, 64, size_in_bytes);
	if(err!=0)
		{
		printf("Memory allocation error");
		exit(1);
		}
#endif
	memset(*pA, 0, size_in_bytes);
	}



/* frees memory */
void d_free(double *pA)
	{
	free( pA );
	}



/* frees aligned memory */
void d_free_align(double *pA)
	{
#if defined(OS_WINDOWS)
	_aligned_free( pA );
#else
	free( pA );
#endif
	}



/* frees aligned memory */
void v_free_align(void *pA)
	{
#if defined(OS_WINDOWS)
	_aligned_free( pA );
#else
	free( pA );
#endif
	}



/* creates a zero matrix aligned */
void d_ones(double **pA, int row, int col)
	{
	void *temp = malloc((row*col)*sizeof(double));
	*pA = temp;
	double *A = *pA;
	int i;
	for(i=0; i<row*col; i++) A[i] = 1.0;
	}



/* creates a zero matrix aligned to a cache line */
void d_ones_align(double **pA, int row, int col)
	{
#if defined(OS_WINDOWS)
	*pA = (double *) _aligned_malloc( (row*col)*sizeof(double), 64 );
#else
	void *temp;
	int err = posix_memalign(&temp, 64, (row*col)*sizeof(double));
	if(err!=0)
		{
		printf("Memory allocation error");
		exit(1);
		}
	*pA = temp;
#endif
	double *A = *pA;
	int i;
	for(i=0; i<row*col; i++) A[i] = 1.0;
	}

#endif



/* creates a zero matrix aligned */
void d_eye(double **pA, int row)
	{
	void *temp = malloc((row*row)*sizeof(double));
	*pA = temp;
	double *A = *pA;
	int i;
	for(i=0; i<row*row; i++) A[i] = 0.0;
	for(i=0; i<row; i++) A[i*(row+1)] = 1.0;
	}



#if ! defined(BLASFEO)

/* prints a matrix */
void d_print_mat(int row, int col, double *A, int lda)
	{
	int i, j;
	for(i=0; i<row; i++)
		{
		for(j=0; j<col; j++)
			{
//			printf("%5.2f ", *(A+i+j*lda));
//			printf("%7.3f ", *(A+i+j*lda));
			printf("%9.5f ", *(A+i+j*lda));
//			printf("%11.7f ", *(A+i+j*lda));
//			printf("%13.9f ", *(A+i+j*lda));
//			printf("%19.15f ", *(A+i+j*lda));
//			printf("%e\t", *(A+i+j*lda));
			}
		printf("\n");
		}
	printf("\n");
	}

void d_print_mat_e(int row, int col, double *A, int lda)
	{
	int i, j;
	for(i=0; i<row; i++)
		{
		for(j=0; j<col; j++)
			{
			printf("%e\t", *(A+i+j*lda));
			}
		printf("\n");
		}
	printf("\n");
	}



/* prints a packed matrix */
void d_print_pmat(int row, int col, int bs, double *A, int sda)
	{

	int ii, i, j, row2;

	for(ii=0; ii<row; ii+=bs)
		{
		row2 = row-ii; if(bs<row2) row2=bs;
		for(i=0; i<row2; i++)
			{
			for(j=0; j<col; j++)
				{
//				printf("%5.2f ", *(A+i+j*lda));
//				printf("%7.3f ", *(A+i+j*bs+ii*sda));
				printf("%9.5f ", *(A+i+j*bs+ii*sda));
//				printf("%11.7f ", *(A+i+j*bs+ii*sda));
//				printf("%13.9f ", *(A+i+j*bs+ii*sda));
//				printf("%19.15f ", *(A+i+j*bs+ii*sda));
//				printf("%e\t", *(A+i+j*lda));
				}
			printf("\n");
			}
//		printf("\n");
		}
	printf("\n");

	}

void d_print_pmat_e(int row, int col, int bs, double *A, int sda)
	{

	int ii, i, j, row2;

	for(ii=0; ii<row; ii+=bs)
		{
		row2 = row-ii; if(bs<row2) row2=bs;
		for(i=0; i<row2; i++)
			{
			for(j=0; j<col; j++)
				{
				printf("%12.4e ", *(A+i+j*bs+ii*sda));
				}
			printf("\n");
			}
//		printf("\n");
		}
	printf("\n");

	}



// prints a banded matrix of band size nb, stored panel-wise in squared blocks of size bs
void d_print_bmat(int nt, int nd, double *pB)
	{

	const int bs = 4;

	int ii, jj, j0, j1;

	double *ptrB = pB;
	double *ptrT = ptrB;

	for(ii=0; ii<nt; ii+=bs)
		{
		j0 = ii-nd+bs; if(j0<0) j0=0;
		j1 = ii+nd; if(j1>nt) j1=nt;

		ptrB = ptrT;
		for(jj=0; jj<j0; jj+=bs)
			{
			// print a row of zeros block
			printf("%9.5f %9.5f %9.5f %9.5f ", 0.0, 0.0, 0.0, 0.0);
			}
		for(jj=j0; jj<j1; jj+=bs)
			{
			// print a row of a block
			printf("%9.5f %9.5f %9.5f %9.5f ", ptrB[0+bs*0], ptrB[0+bs*1], ptrB[0+bs*2], ptrB[0+bs*3]);
			ptrB += bs*bs;
			}
		for(jj=j1; jj<nt; jj+=bs)
			{
			// print a row of zeros block
			printf("%9.5f %9.5f %9.5f %9.5f ", 0.0, 0.0, 0.0, 0.0);
			}
		printf("\n");

		ptrB = ptrT;
		for(jj=0; jj<j0; jj+=bs)
			{
			// print a row of zeros block
			printf("%9.5f %9.5f %9.5f %9.5f ", 0.0, 0.0, 0.0, 0.0);
			}
		for(jj=j0; jj<j1; jj+=bs)
			{
			// print a row of a block
			printf("%9.5f %9.5f %9.5f %9.5f ", ptrB[1+bs*0], ptrB[1+bs*1], ptrB[1+bs*2], ptrB[1+bs*3]);
			ptrB += bs*bs;
			}
		for(jj=j1; jj<nt; jj+=bs)
			{
			// print a row of zeros block
			printf("%9.5f %9.5f %9.5f %9.5f ", 0.0, 0.0, 0.0, 0.0);
			}
		printf("\n");

		ptrB = ptrT;
		for(jj=0; jj<j0; jj+=bs)
			{
			// print a row of zeros block
			printf("%9.5f %9.5f %9.5f %9.5f ", 0.0, 0.0, 0.0, 0.0);
			}
		for(jj=j0; jj<j1; jj+=bs)
			{
			// print a row of a block
			printf("%9.5f %9.5f %9.5f %9.5f ", ptrB[2+bs*0], ptrB[2+bs*1], ptrB[2+bs*2], ptrB[2+bs*3]);
			ptrB += bs*bs;
			}
		for(jj=j1; jj<nt; jj+=bs)
			{
			// print a row of zeros block
			printf("%9.5f %9.5f %9.5f %9.5f ", 0.0, 0.0, 0.0, 0.0);
			}
		printf("\n");

		ptrB = ptrT;
		for(jj=0; jj<j0; jj+=bs)
			{
			// print a row of zeros block
			printf("%9.5f %9.5f %9.5f %9.5f ", 0.0, 0.0, 0.0, 0.0);
			}
		for(jj=j0; jj<j1; jj+=bs)
			{
			// print a row of a block
			printf("%9.5f %9.5f %9.5f %9.5f ", ptrB[3+bs*0], ptrB[3+bs*1], ptrB[3+bs*2], ptrB[3+bs*3]);
			ptrB += bs*bs;
			}
		for(jj=j1; jj<nt; jj+=bs)
			{
			// print a row of zeros block
			printf("%9.5f %9.5f %9.5f %9.5f ", 0.0, 0.0, 0.0, 0.0);
			}
		printf("\n");

		ptrT = ptrB;

		}
	printf("\n");

	}

#endif
