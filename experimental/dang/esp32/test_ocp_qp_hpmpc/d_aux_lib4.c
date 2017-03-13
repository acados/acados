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

#include <math.h>

#include "block_size.h"
#include "kernel_d_lib4.h"
// 2017.03.13 Dang add target.h
#include "target.h"

#if defined(TARGET_X64_AVX2) || defined(TARGET_X64_AVX)
#include <mmintrin.h>
#include <xmmintrin.h>  // SSE
#include <emmintrin.h>  // SSE2
#include <pmmintrin.h>  // SSE3
#include <smmintrin.h>  // SSE4
#include <immintrin.h>  // AVX
#endif


/* copies a matrix */
void d_rep_mat(int reps, int row, int col, double *A, int lda, double *B, int ldb)
	{

	int i, j, l;

	for(l=0; l<reps; l++)
		{
		for(j=0; j<col; j++)
			{
			for(i=0; i<row; i++)
				{
				B[i+j*ldb+l*row*col] = A[i+j*lda];
				}
			}
		}

	}



/* copies and scales a matrix */
void dadd_mat(int row, int col, double alpha, double *A, int lda, double *B, int ldb)
	{

	int i, j;

	for(j=0; j<col; j++)
		{
		i = 0;
		for(; i<row-3; i+=4)
			{
			B[i+0+j*ldb] += alpha*A[i+0+j*lda];
			B[i+1+j*ldb] += alpha*A[i+1+j*lda];
			B[i+2+j*ldb] += alpha*A[i+2+j*lda];
			B[i+3+j*ldb] += alpha*A[i+3+j*lda];
			}
		for(; i<row; i++)
			{
			B[i+j*ldb] += alpha*A[i+j*lda];
			}
		}

	}



void dax_mat(int row, int col, double alpha, double *A, int lda, double *B, int ldb)
	{

	int i, j;

	for(j=0; j<col; j++)
		{
		i = 0;
		for(; i<row-3; i+=4)
			{
			B[i+0+j*ldb] = alpha*A[i+0+j*lda];
			B[i+1+j*ldb] = alpha*A[i+1+j*lda];
			B[i+2+j*ldb] = alpha*A[i+2+j*lda];
			B[i+3+j*ldb] = alpha*A[i+3+j*lda];
			}
		for(; i<row; i++)
			{
			B[i+j*ldb] = alpha*A[i+j*lda];
			}
		}

	}



float d_max_mat(int row, int col, double *A, int lda)
	{

	if(row<=0 || col<=0 )
		return 0.0;

	int i, j;

	static float max_vec[4];

	max_vec[0] = A[0];
	max_vec[1] = A[0];
	max_vec[2] = A[0];
	max_vec[3] = A[0];

	for(j=0; j<col; j++)
		{
		i = 0;
		for(; i<row-3; i+=4)
			{
			max_vec[0] = fmax( max_vec[0], A[i+0+j*lda] );
			max_vec[1] = fmax( max_vec[1], A[i+1+j*lda] );
			max_vec[2] = fmax( max_vec[2], A[i+2+j*lda] );
			max_vec[3] = fmax( max_vec[3], A[i+3+j*lda] );
			}
		for(; i<row; i++)
			{
			max_vec[0] = fmax( max_vec[0], A[i+j*lda] );
			}
		}

	max_vec[0] = fmax( max_vec[0], max_vec[2] );
	max_vec[1] = fmax( max_vec[1], max_vec[3] );
	max_vec[0] = fmax( max_vec[0], max_vec[1] );

	return max_vec[0];

	}



float d_min_mat(int row, int col, double *A, int lda)
	{

	if(row<=0 || col<=0 )
		return 0.0;

	int i, j;

	static float min_vec[4];

	min_vec[0] = A[0];
	min_vec[1] = A[0];
	min_vec[2] = A[0];
	min_vec[3] = A[0];

	for(j=0; j<col; j++)
		{
		i = 0;
		for(; i<row-3; i+=4)
			{
			min_vec[0] = fmin( min_vec[0], A[i+0+j*lda] );
			min_vec[1] = fmin( min_vec[1], A[i+1+j*lda] );
			min_vec[2] = fmin( min_vec[2], A[i+2+j*lda] );
			min_vec[3] = fmin( min_vec[3], A[i+3+j*lda] );
			}
		for(; i<row; i++)
			{
			min_vec[0] = fmin( min_vec[0], A[i+j*lda] );
			}
		}

	min_vec[0] = fmin( min_vec[0], min_vec[2] );
	min_vec[1] = fmin( min_vec[1], min_vec[3] );
	min_vec[0] = fmin( min_vec[0], min_vec[1] );

	return min_vec[0];

	}



void d_set_mat(int row, int col, double alpha, double *A, int lda)
	{

	int i, j;

	for(j=0; j<col; j++)
		{
		i = 0;
		for(; i<row-3; i+=4)
			{
			A[i+0+j*lda] = alpha;
			A[i+1+j*lda] = alpha;
			A[i+2+j*lda] = alpha;
			A[i+3+j*lda] = alpha;
			}
		for(; i<row; i++)
			{
			A[i+j*lda] = alpha;
			}
		}

	}



void dgeset_lib(int row, int col, double alpha, int offset, double *pA, int sda)
	{

	const int bs = 4;

	int ii, jj;

	int na = (bs-offset%bs)%bs;
	na = row<na ? row : na;

	ii = 0;
	if(na>0)
		{
		for(; ii<na; ii++)
			{
			for(jj=0; jj<col; jj++)
				{
				pA[jj*bs] = alpha;
				}
			pA += 1;
			}
		pA += bs*(sda-1);
		}
	for( ; ii<row-3; ii+=4)
		{
#if defined(TARGET_X64_AVX2) || defined(TARGET_X64_AVX)
		kernel_dgeset_4_lib4(col, alpha, pA);
#else
		for(jj=0; jj<col; jj++)
			{
			pA[0+jj*bs] = alpha;
			pA[1+jj*bs] = alpha;
			pA[2+jj*bs] = alpha;
			pA[3+jj*bs] = alpha;
			}
#endif
		pA += bs*sda;
		}
	for(; ii<row; ii++)
		{
		for(jj=0; jj<col; jj++)
			{
			pA[jj*bs] = alpha;
			}
		pA += 1;
		}

	}



void dtrset_lib(int row, double alpha, int offset, double *pA, int sda)
	{

	const int bs = 4;

	int ii, jj;

	int col = row;

	int na = (bs-offset%bs)%bs;
	na = row<na ? row : na;

	ii = 0;
	if(na>0)
		{
		for(; ii<na; ii++)
			{
			for(jj=0; jj<=ii; jj++)
				{
				pA[jj*bs] = alpha;
				}
			pA += 1;
			}
		pA += bs*(sda-1);
		}
	for( ; ii<row-3; ii+=4)
		{
#if defined(TARGET_X64_AVX2) || defined(TARGET_X64_AVX)
		kernel_dtrset_4_lib4(ii, alpha, pA);
#else
		for(jj=0; jj<ii; jj++)
			{
			pA[0+jj*bs] = alpha;
			pA[1+jj*bs] = alpha;
			pA[2+jj*bs] = alpha;
			pA[3+jj*bs] = alpha;
			}
		pA[0+(jj+0)*bs] = alpha;
		pA[1+(jj+0)*bs] = alpha;
		pA[2+(jj+0)*bs] = alpha;
		pA[3+(jj+0)*bs] = alpha;
		pA[1+(jj+1)*bs] = alpha;
		pA[2+(jj+1)*bs] = alpha;
		pA[3+(jj+1)*bs] = alpha;
		pA[2+(jj+2)*bs] = alpha;
		pA[3+(jj+2)*bs] = alpha;
		pA[3+(jj+3)*bs] = alpha;
#endif
		pA += bs*sda;
		}
	for(; ii<row; ii++)
		{
		for(jj=0; jj<ii+1; jj++)
			{
			pA[jj*bs] = alpha;
			}
		pA += 1;
		}

	}



void d_scale_mat(int row, int col, double alpha, double *A, int lda)
	{

	int i, j;

	for(j=0; j<col; j++)
		{
		i = 0;
		for(; i<row-3; i+=4)
			{
			A[i+0+j*lda] *= alpha;
			A[i+1+j*lda] *= alpha;
			A[i+2+j*lda] *= alpha;
			A[i+3+j*lda] *= alpha;
			}
		for(; i<row; i++)
			{
			A[i+j*lda] *= alpha;
			}
		}

	}



void d_scale_pmat(int row, int col, double alpha, int offset, double *pA, int sda)
	{

	const int bs = 4;

	int ii, jj;

	int na = (bs-offset%bs)%bs;
	na = row<na ? row : na;

	ii = 0;
	if(na>0)
		{
		for(; ii<na; ii++)
			{
			for(jj=0; jj<col; jj++)
				{
				pA[jj*bs] *= alpha;
				}
			pA += 1;
			}
		pA += bs*(sda-1);
		}
	for( ; ii<row-3; ii+=4)
		{
		for(jj=0; jj<col; jj++)
			{
			pA[0+jj*bs] *= alpha;
			pA[1+jj*bs] *= alpha;
			pA[2+jj*bs] *= alpha;
			pA[3+jj*bs] *= alpha;
			}
		pA += bs*sda;
		}
	for(; ii<row; ii++)
		{
		for(jj=0; jj<col; jj++)
			{
			pA[jj*bs] *= alpha;
			}
		pA += 1;
		}

	}



/* copies a matrix */
void d_copy_mat(int row, int col, double *A, int lda, double *B, int ldb)
	{

	int i, j;

	for(j=0; j<col; j++)
		{
		i = 0;
		for(; i<row-3; i+=4)
			{
			B[i+0+j*ldb] = A[i+0+j*lda];
			B[i+1+j*ldb] = A[i+1+j*lda];
			B[i+2+j*ldb] = A[i+2+j*lda];
			B[i+3+j*ldb] = A[i+3+j*lda];
			}
		for(; i<row; i++)
			{
			B[i+j*ldb] = A[i+j*lda];
			}
		}

	}



/* copies the transpose of a matrix */
void d_tran_mat(int row, int col, double *A, int lda, double *B, int ldb)
	{

	int i, j;

	for(j=0; j<col; j++)
		{
		i = 0;
		for(; i<row-3; i+=4)
			{
			B[j+(i+0)*ldb] = A[i+0+j*lda];
			B[j+(i+1)*ldb] = A[i+1+j*lda];
			B[j+(i+2)*ldb] = A[i+2+j*lda];
			B[j+(i+3)*ldb] = A[i+3+j*lda];
			}
		for(; i<row; i++)
			{
			B[j+i*ldb] = A[i+j*lda];
			}
		}

	}



/* copies a packed matrix */
void d_copy_pmat(int row, int col, int bs_dummy, double *A, int sda, double *B, int sdb)
	{

	const int bs = 4;

	int i, ii, j, row2;

	ii = 0;
	for(; ii<row-3; ii+=bs)
		{
		for(j=0; j<col; j++)
			{
			B[0+j*bs+ii*sdb] = A[0+j*bs+ii*sda];
			B[1+j*bs+ii*sdb] = A[1+j*bs+ii*sda];
			B[2+j*bs+ii*sdb] = A[2+j*bs+ii*sda];
			B[3+j*bs+ii*sdb] = A[3+j*bs+ii*sda];
			}
		}
	if(ii<row)
		{
		row2 = row-ii;
		for(j=0; j<col; j++)
			{
			for(i=0; i<row2; i++)
				{
				B[i+j*bs+ii*sdb] = A[i+j*bs+ii*sda];
				}
			}
		}

	}



#if 0
// copies a packed matrix into a packed matrix
void d_copy_pmat_general(int m, int n, int offset_A, double *A, int sda, int offset_B, double *B, int sdb)
	{

	if(m<=0 || n<=0)
		return;

	const int bs = D_MR;

	int mna, ii;

	int offA = offset_A%bs;
	int offB = offset_B%bs;

	// A at the beginning of the block
	A -= offA;

	// A at the beginning of the block
	B -= offB;

	// same alignment
	if(offA==offB)
		{
		ii = 0;
		// clean up at the beginning
		mna = (4-offB)%bs;
		if(mna>0)
			{
			if(m<mna) // mna<=3  ==>  m = { 1, 2 }
				{
				if(m==1)
					{
					kernel_align_panel_1_0_lib4(n, A+offA, B+offB);
					return;
					}
				else //if(m==2 && mna==3)
					{
					kernel_align_panel_2_0_lib4(n, A+offA, B+offB);
					return;
					}
				}
			if(mna==1)
				{
				kernel_align_panel_1_0_lib4(n, A+offA, B+offB);
				A += 4*sda;
				B += 4*sdb;
				ii += 1;
				}
			else if(mna==2)
				{
				kernel_align_panel_2_0_lib4(n, A+offA, B+offB);
				A += 4*sda;
				B += 4*sdb;
				ii += 2;
				}
			else // if(mna==3)
				{
				kernel_align_panel_3_0_lib4(n, A+offA, B+offB);
				A += 4*sda;
				B += 4*sdb;
				ii += 3;
				}
			}
		// main loop
#if defined(TARGET_X64_AVX) || defined(TARGET_X64_AVX2)
		for(; ii<m-7; ii+=8)
			{
			kernel_align_panel_8_0_lib4(n, A, sda, B, sdb);
			A += 8*sda;
			B += 8*sdb;
			}
#endif
		for(; ii<m-3; ii+=4)
			{
			kernel_align_panel_4_0_lib4(n, A, B);
			A += 4*sda;
			B += 4*sdb;
			}
		// clean up at the end
		if(ii<m)
			{
			if(m-ii==1)
				kernel_align_panel_1_0_lib4(n, A, B);
			else if(m-ii==2)
				kernel_align_panel_2_0_lib4(n, A, B);
			else // if(m-ii==3)
				kernel_align_panel_3_0_lib4(n, A, B);
			}
		}
	// skip one element of A
	else if(offA==(offB+1)%bs)
		{
		ii = 0;
		// clean up at the beginning
		mna = (4-offB)%bs;
		if(mna>0)
			{
			if(m<mna) // mna<=3  ==>  m = { 1, 2 }
				{
				if(m==1)
					{
					kernel_align_panel_1_0_lib4(n, A+offA, B+offB);
					return;
					}
				else //if(m==2 && mna==3)
					{
					kernel_align_panel_2_0_lib4(n, A+offA, B+offB);
					return;
					}
				}
			if(mna==1)
				{
				kernel_align_panel_1_0_lib4(n, A+offA, B+offB);
				//A += 4*sda;
				B += 4*sdb;
				ii += 1;
				}
			else if(mna==2)
				{
				kernel_align_panel_2_3_lib4(n, A, sda, B+2);
				A += 4*sda;
				B += 4*sdb;
				ii += 2;
				}
			else // if(mna==3)
				{
				kernel_align_panel_3_2_lib4(n, A, sda, B+1);
				A += 4*sda;
				B += 4*sdb;
				ii += 3;
				}
			}
		// main loop
#if defined(TARGET_X64_AVX) || defined(TARGET_X64_AVX2)
		for( ; ii<m-7; ii+=8)
			{
			kernel_align_panel_8_1_lib4(n, A, sda, B, sdb);
			A += 8*sda;
			B += 8*sdb;
			}
#endif
		for( ; ii<m-3; ii+=4)
			{
			kernel_align_panel_4_1_lib4(n, A, sda, B);
			A += 4*sda;
			B += 4*sdb;
			}
		// clean up at the end
		if(ii<m)
			{
			if(m-ii==1)
				kernel_align_panel_1_0_lib4(n, A+1, B);
			else if(m-ii==2)
				kernel_align_panel_2_0_lib4(n, A+1, B);
			else // if(m-ii==3)
				kernel_align_panel_3_0_lib4(n, A+1, B);
			}
		}
	// skip 2 elements of A
	else if(offA==(offB+2)%bs)
		{
		ii = 0;
		// clean up at the beginning
		mna = (4-offB)%bs;
		if(mna>0)
			{
			if(m<mna)
				{
				if(m==1)
					{
					kernel_align_panel_1_0_lib4(n, A+offA, B+offB);
					return;
					}
				else // if(m==2 && mna==3)
					{
					kernel_align_panel_2_3_lib4(n, A, sda, B+1);
					return;
					}
				}
			if(mna==1)
				{
				kernel_align_panel_1_0_lib4(n, A+1, B+3);
				// A += 4*sda;
				B += 4*sdb;
				ii += 1;
				}
			else if(mna==2)
				{
				kernel_align_panel_2_0_lib4(n, A, B+2);
				// A += 4*sda;
				B += 4*sdb;
				ii += 2;
				}
			else // if(mna==3)
				{
				kernel_align_panel_3_3_lib4(n, A, sda, B+1);
				A += 4*sda;
				B += 4*sdb;
				ii += 3;
				}
			}
		// main loop
#if defined(TARGET_X64_AVX) || defined(TARGET_X64_AVX2)
		for(; ii<m-7; ii+=8)
			{
			kernel_align_panel_8_2_lib4(n, A, sda, B, sdb);
			A += 8*sda;
			B += 8*sdb;
			}
#endif
		for(; ii<m-3; ii+=4)
			{
			kernel_align_panel_4_2_lib4(n, A, sda, B);
			A += 4*sda;
			B += 4*sdb;
			}
		// clean up at the end
		if(ii<m)
			{
			if(m-ii==1)
				kernel_align_panel_1_0_lib4(n, A+2, B);
			else if(m-ii==2)
				kernel_align_panel_2_0_lib4(n, A+2, B);
			else // if(m-ii==3)
				kernel_align_panel_3_2_lib4(n, A, sda, B);
			}
		}
	// skip 3 elements of A
	else // if(offA==(offB+3)%bs)
		{
		ii = 0;
		// clean up at the beginning
		mna = (4-offB)%bs;
		if(mna>0)
			{
			if(m<mna)
				{
				if(m==1)
					{
					kernel_align_panel_1_0_lib4(n, A+offA, B+offB);
					return;
					}
				else // if(m==2 && mna==3)
					{
					kernel_align_panel_2_0_lib4(n, A+offA, B+offB);
					return;
					}
				}
			if(mna==1)
				{
				kernel_align_panel_1_0_lib4(n, A+offA, B+offB);
				// A += 4*sda;
				B += 4*sdb;
				ii += 1;
				}
			else if(mna==2)
				{
				kernel_align_panel_2_0_lib4(n, A+offA, B+offB);
				// A += 4*sda;
				B += 4*sdb;
				ii += 2;
				}
			else // if(mna==3)
				{
				kernel_align_panel_3_0_lib4(n, A+offA, B+offB);
				// A += 4*sda;
				B += 4*sdb;
				ii += 3;
				}
			}
		// main loop
#if defined(TARGET_X64_AVX) || defined(TARGET_X64_AVX2)
		for(; ii<m-7; ii+=8)
			{
			kernel_align_panel_8_3_lib4(n, A, sda, B, sdb);
			A += 8*sda;
			B += 8*sdb;
			}
#endif
		for(; ii<m-3; ii+=4)
			{
			kernel_align_panel_4_3_lib4(n, A, sda, B);
			A += 4*sda;
			B += 4*sdb;
			}
		// clean up at the end
		if(ii<m)
			{
			if(m-ii==1)
				kernel_align_panel_1_0_lib4(n, A+3, B);
			else if(m-ii==2)
				kernel_align_panel_2_3_lib4(n, A, sda, B);
			else // if(m-ii==3)
				kernel_align_panel_3_3_lib4(n, A, sda, B);
			}
		}

	}
#endif



/* copies 1 to 4 rows from a block of a packed matrix into a (misalinged) packed matrix */
void d_copy_pmat_panel(int row, int col, int offset, double *A, double *B, int sdb)
	{

	if(row<=0 || col<=0)
		return;

	const int bs = 4;

	int i, ii;

	int row0 = bs - offset%bs;
	row0 = row0>row ? row : row0;
	int row1 = row - row0;

	double *(ptrB[4]);
	for(i=0; i<row0; i++)
		ptrB[i] = B + i;
	for(i=0; i<row1; i++)
		ptrB[row0+i] = B + row0 + (sdb-1)*bs + i;


	if(row==1)
		{
		i = 0;
		for(; i<col-3; i+=4)
			{
			ptrB[0][(i+0)*bs] = A[0+(i+0)*bs];
			ptrB[0][(i+1)*bs] = A[0+(i+1)*bs];
			ptrB[0][(i+2)*bs] = A[0+(i+2)*bs];
			ptrB[0][(i+3)*bs] = A[0+(i+3)*bs];
			}
		for(; i<col; i++)
			{
			ptrB[0][(i+0)*bs] = A[0+(i+0)*bs];
			}
		}
	else if(row==2)
		{
		i = 0;
		for(; i<col-3; i+=4)
			{
			ptrB[0][(i+0)*bs] = A[0+(i+0)*bs];
			ptrB[1][(i+0)*bs] = A[1+(i+0)*bs];

			ptrB[0][(i+1)*bs] = A[0+(i+1)*bs];
			ptrB[1][(i+1)*bs] = A[1+(i+1)*bs];

			ptrB[0][(i+2)*bs] = A[0+(i+2)*bs];
			ptrB[1][(i+2)*bs] = A[1+(i+2)*bs];

			ptrB[0][(i+3)*bs] = A[0+(i+3)*bs];
			ptrB[1][(i+3)*bs] = A[1+(i+3)*bs];
			}
		for(; i<col; i++)
			{
			ptrB[0][(i+0)*bs] = A[0+(i+0)*bs];
			ptrB[1][(i+0)*bs] = A[1+(i+0)*bs];
			}
		}
	else if(row==3)
		{
		i = 0;
		for(; i<col-3; i+=4)
			{
			ptrB[0][(i+0)*bs] = A[0+(i+0)*bs];
			ptrB[1][(i+0)*bs] = A[1+(i+0)*bs];
			ptrB[2][(i+0)*bs] = A[2+(i+0)*bs];

			ptrB[0][(i+1)*bs] = A[0+(i+1)*bs];
			ptrB[1][(i+1)*bs] = A[1+(i+1)*bs];
			ptrB[2][(i+1)*bs] = A[2+(i+1)*bs];

			ptrB[0][(i+2)*bs] = A[0+(i+2)*bs];
			ptrB[1][(i+2)*bs] = A[1+(i+2)*bs];
			ptrB[2][(i+2)*bs] = A[2+(i+2)*bs];

			ptrB[0][(i+3)*bs] = A[0+(i+3)*bs];
			ptrB[1][(i+3)*bs] = A[1+(i+3)*bs];
			ptrB[2][(i+3)*bs] = A[2+(i+3)*bs];
			}
		for(; i<col; i++)
			{
			ptrB[0][(i+0)*bs] = A[0+(i+0)*bs];
			ptrB[1][(i+0)*bs] = A[1+(i+0)*bs];
			ptrB[2][(i+0)*bs] = A[2+(i+0)*bs];
			}
		}
	else //if(row==4)
		{
		i = 0;
		for(; i<col-3; i+=4)
			{
			ptrB[0][(i+0)*bs] = A[0+(i+0)*bs];
			ptrB[1][(i+0)*bs] = A[1+(i+0)*bs];
			ptrB[2][(i+0)*bs] = A[2+(i+0)*bs];
			ptrB[3][(i+0)*bs] = A[3+(i+0)*bs];

			ptrB[0][(i+1)*bs] = A[0+(i+1)*bs];
			ptrB[1][(i+1)*bs] = A[1+(i+1)*bs];
			ptrB[2][(i+1)*bs] = A[2+(i+1)*bs];
			ptrB[3][(i+0)*bs] = A[3+(i+0)*bs];

			ptrB[0][(i+2)*bs] = A[0+(i+2)*bs];
			ptrB[1][(i+2)*bs] = A[1+(i+2)*bs];
			ptrB[2][(i+2)*bs] = A[2+(i+2)*bs];
			ptrB[3][(i+0)*bs] = A[3+(i+0)*bs];

			ptrB[0][(i+3)*bs] = A[0+(i+3)*bs];
			ptrB[1][(i+3)*bs] = A[1+(i+3)*bs];
			ptrB[2][(i+3)*bs] = A[2+(i+3)*bs];
			ptrB[3][(i+0)*bs] = A[3+(i+0)*bs];
			}
		for(; i<col; i++)
			{
			ptrB[0][(i+0)*bs] = A[0+(i+0)*bs];
			ptrB[1][(i+0)*bs] = A[1+(i+0)*bs];
			ptrB[2][(i+0)*bs] = A[2+(i+0)*bs];
			ptrB[3][(i+0)*bs] = A[3+(i+0)*bs];
			}
		}


	}



/* copies 1 to 4 rows from a (misaligned) packed matrix into a block of a packed matrix */
void d_align_pmat_panel(int row, int col, int offset, double *A, int sda, double *B)
	{

	if(row<=0 || col<=0)
		return;

	const int bs = 4;

	int i, ii;

	int row0 = bs - offset%bs;
	row0 = row0>row ? row : row0;
	int row1 = row - row0;

	double *(ptrA[4]);
	for(i=0; i<row0; i++)
		ptrA[i] = A + i;
	for(i=0; i<row1; i++)
		ptrA[row0+i] = A + row0 + (sda-1)*bs + i;


	if(row==1)
		{
		i = 0;
		for(; i<col-3; i+=4)
			{
			B[0+(i+0)*bs] = ptrA[0][(i+0)*bs];
			B[0+(i+1)*bs] = ptrA[0][(i+1)*bs];
			B[0+(i+2)*bs] = ptrA[0][(i+2)*bs];
			B[0+(i+3)*bs] = ptrA[0][(i+3)*bs];
			}
		for(; i<col; i++)
			{
			B[0+(i+0)*bs] = ptrA[0][(i+0)*bs];
			}
		}
	else if(row==2)
		{
		i = 0;
		for(; i<col-3; i+=4)
			{
			B[0+(i+0)*bs] = ptrA[0][(i+0)*bs];
			B[1+(i+0)*bs] = ptrA[1][(i+0)*bs];

			B[0+(i+1)*bs] = ptrA[0][(i+1)*bs];
			B[1+(i+1)*bs] = ptrA[1][(i+1)*bs];

			B[0+(i+2)*bs] = ptrA[0][(i+2)*bs];
			B[1+(i+2)*bs] = ptrA[1][(i+2)*bs];

			B[0+(i+3)*bs] = ptrA[0][(i+3)*bs];
			B[1+(i+3)*bs] = ptrA[1][(i+3)*bs];
			}
		for(; i<col; i++)
			{
			B[0+(i+0)*bs] = ptrA[0][(i+0)*bs];
			B[1+(i+0)*bs] = ptrA[1][(i+0)*bs];
			}
		}
	else if(row==3)
		{
		i = 0;
		for(; i<col-3; i+=4)
			{
			B[0+(i+0)*bs] = ptrA[0][(i+0)*bs];
			B[1+(i+0)*bs] = ptrA[1][(i+0)*bs];
			B[2+(i+0)*bs] = ptrA[2][(i+0)*bs];

			B[0+(i+1)*bs] = ptrA[0][(i+1)*bs];
			B[1+(i+1)*bs] = ptrA[1][(i+1)*bs];
			B[2+(i+1)*bs] = ptrA[2][(i+1)*bs];

			B[0+(i+2)*bs] = ptrA[0][(i+2)*bs];
			B[1+(i+2)*bs] = ptrA[1][(i+2)*bs];
			B[2+(i+2)*bs] = ptrA[2][(i+2)*bs];

			B[0+(i+3)*bs] = ptrA[0][(i+3)*bs];
			B[1+(i+3)*bs] = ptrA[1][(i+3)*bs];
			B[2+(i+3)*bs] = ptrA[2][(i+3)*bs];
			}
		for(; i<col; i++)
			{
			B[0+(i+0)*bs] = ptrA[0][(i+0)*bs];
			B[1+(i+0)*bs] = ptrA[1][(i+0)*bs];
			B[2+(i+0)*bs] = ptrA[2][(i+0)*bs];
			}
		}
	else //if(row==4)
		{
		i = 0;
		for(; i<col-3; i+=4)
			{
			B[0+(i+0)*bs] = ptrA[0][(i+0)*bs];
			B[1+(i+0)*bs] = ptrA[1][(i+0)*bs];
			B[2+(i+0)*bs] = ptrA[2][(i+0)*bs];
			B[3+(i+0)*bs] = ptrA[3][(i+0)*bs];

			B[0+(i+1)*bs] = ptrA[0][(i+1)*bs];
			B[1+(i+1)*bs] = ptrA[1][(i+1)*bs];
			B[2+(i+1)*bs] = ptrA[2][(i+1)*bs];
			B[3+(i+0)*bs] = ptrA[3][(i+0)*bs];

			B[0+(i+2)*bs] = ptrA[0][(i+2)*bs];
			B[1+(i+2)*bs] = ptrA[1][(i+2)*bs];
			B[2+(i+2)*bs] = ptrA[2][(i+2)*bs];
			B[3+(i+0)*bs] = ptrA[3][(i+0)*bs];

			B[0+(i+3)*bs] = ptrA[0][(i+3)*bs];
			B[1+(i+3)*bs] = ptrA[1][(i+3)*bs];
			B[2+(i+3)*bs] = ptrA[2][(i+3)*bs];
			B[3+(i+0)*bs] = ptrA[3][(i+0)*bs];
			}
		for(; i<col; i++)
			{
			B[0+(i+0)*bs] = ptrA[0][(i+0)*bs];
			B[1+(i+0)*bs] = ptrA[1][(i+0)*bs];
			B[2+(i+0)*bs] = ptrA[2][(i+0)*bs];
			B[3+(i+0)*bs] = ptrA[3][(i+0)*bs];
			}
		}


	}



// copies a lower triangular packed matrix
void d_copy_pmat_l(int row, int bs_dummy, double *A, int sda, double *B, int sdb)
	{

	const int bs = 4;

	int i, ii, j, row2, row0;

	ii = 0;
	for(; ii<row-3; ii+=bs)
		{
		j = 0;
		for(; j<ii; j++)
			{
			B[0+j*bs+ii*sdb] = A[0+j*bs+ii*sda];
			B[1+j*bs+ii*sdb] = A[1+j*bs+ii*sda];
			B[2+j*bs+ii*sdb] = A[2+j*bs+ii*sda];
			B[3+j*bs+ii*sdb] = A[3+j*bs+ii*sda];
			}
		for(; j<ii+bs; j++)
			{
			row0 = j-ii;
			if(row0<0) row0=0;
			for(i=row0; i<bs; i++)
				{
				B[i+j*bs+ii*sdb] = A[i+j*bs+ii*sda];
				}
			}
		}
	if(ii<row)
		{
		row2 = row-ii;
		if(bs<row2) row2 = bs;
		for(j=0; j<ii+row2; j++)
			{
			row0 = j-ii;
			if(row0<0) row0=0;
			for(i=row0; i<row2; i++)
				{
				B[i+j*bs+ii*sdb] = A[i+j*bs+ii*sda];
				}
			}
		}

	}



/* transposes a lower triangular packed matrix */
/*void d_transpose_pmat_lo(int row, int offset, double *A, int sda, double *B, int sdb)*/
/*	{*/
/*	*/
/*	const int bs = 4;*/

/*	int i, j, jj;*/
/*	*/
/*	int row0, row1, row2, row3;*/
/*	row0 = (bs-offset%bs)%bs; // row2 < bs !!!*/
/*	*/
/*	double *pA, *pB;*/

/*	jj = 0;*/
/*	for(; jj<row-3; jj+=4)*/
/*		{*/
/*		row1 = row - jj;*/
/*		pA = A + jj*bs + jj*sda;*/
/*		pB = B + jj*bs + jj*sdb;*/
/*		row2 = row0; // row2 < bs !!!*/
/*		if(row1<row2)*/
/*			row2 = row1;*/
/*		i = 0;*/
/*		if(row2>0)*/
/*			{*/
/*			for(; i<row2; i++)*/
/*				{*/
/*				for(j=0; j<=i; j++)*/
/*					{*/
/*					pB[j] = pA[j*bs];*/
/*					}*/
/*				pA += 1;*/
/*				pB += bs;*/
/*				}*/
/*			pA += (sda-1)*bs;*/
/*			}*/
/*		row3 = row2 + 4;*/
/*		if(row1<row3)*/
/*			row3 = row1;*/
/*		row2 = 4;*/
/*		if(row1<row2)*/
/*			row2 = row1;*/
/*		for(; i<row2; i++)*/
/*			{*/
/*			for(j=0; j<=i; j++)*/
/*				{*/
/*				pB[j] = pA[j*bs];*/
/*				}*/
/*			pA += 1;*/
/*			pB += bs;*/
/*			}*/
/*		for(; i<row3; i++)*/
/*			{*/
/*			pB[0] = pA[0*bs];*/
/*			pB[1] = pA[1*bs];*/
/*			pB[2] = pA[2*bs];*/
/*			pB[3] = pA[3*bs];*/
/*			pA += 1;*/
/*			pB += bs;*/
/*			}*/
/*		pA += (sda-1)*bs;*/
/*		for(; i<row1-3; i+=4)*/
/*			{*/
/*			// buildin_prefetch*/
/*			// unroll 0*/
/*			pB[0+0*bs] = pA[0+0*bs];*/
/*			pB[1+0*bs] = pA[0+1*bs];*/
/*			pB[2+0*bs] = pA[0+2*bs];*/
/*			pB[3+0*bs] = pA[0+3*bs];*/
/*			// unroll 1*/
/*			pB[0+1*bs] = pA[1+0*bs];*/
/*			pB[1+1*bs] = pA[1+1*bs];*/
/*			pB[2+1*bs] = pA[1+2*bs];*/
/*			pB[3+1*bs] = pA[1+3*bs];*/
/*			// unroll 2*/
/*			pB[0+2*bs] = pA[2+0*bs];*/
/*			pB[1+2*bs] = pA[2+1*bs];*/
/*			pB[2+2*bs] = pA[2+2*bs];*/
/*			pB[3+2*bs] = pA[2+3*bs];*/
/*			// unroll 3*/
/*			pB[0+3*bs] = pA[3+0*bs];*/
/*			pB[1+3*bs] = pA[3+1*bs];*/
/*			pB[2+3*bs] = pA[3+2*bs];*/
/*			pB[3+3*bs] = pA[3+3*bs];*/
/*			pA += sda*bs;*/
/*			pB += 4*bs;*/
/*			}*/
/*		for(; i<row1; i++)*/
/*			{*/
/*			pB[0] = pA[0*bs];*/
/*			pB[1] = pA[1*bs];*/
/*			pB[2] = pA[2*bs];*/
/*			pB[3] = pA[3*bs];*/
/*			pA += 1;*/
/*			pB += bs;*/
/*			}*/
/*		}*/
/*	if(jj<row)*/
/*		{*/
/*		row1 = row - jj;*/
/*		pA = A + jj*bs + jj*sda;*/
/*		pB = B + jj*bs + jj*sdb;*/
/*		row2 = row0; // row2 < bs !!!*/
/*		if(row1<row2)*/
/*			row2 = row1;*/
/*		i = 0;*/
/*		if(row2>0)*/
/*			{*/
/*			for(; i<row2; i++)*/
/*				{*/
/*				for(j=0; j<=i; j++)*/
/*					{*/
/*					pB[j] = pA[j*bs];*/
/*					}*/
/*				pA += 1;*/
/*				pB += bs;*/
/*				}*/
/*			pA += (sda-1)*bs;*/
/*			}*/
/*		row2 = 4;*/
/*		if(row1<row2)*/
/*			row2 = row1;*/
/*		for(; i<row2; i++)*/
/*			{*/
/*			for(j=0; j<=i; j++)*/
/*				{*/
/*				pB[j] = pA[j*bs];*/
/*				}*/
/*			pA += 1;*/
/*			pB += bs;*/
/*			}*/
/*		}*/
/*	*/
/*	}*/



// copies a packed matrix into an aligned packed matrix ; A has to be aligned at the beginning of the current block : the offset takes care of the row to be copied TODO IMPROVE IMPLEMENTATION !!!
void d_align_pmat(int row, int col, int offset, int bs_dummy, double *A, int sda, double *B, int sdb)
	{

	const int bs = 4;

	int i, j;

	double *ptrA, *ptrB;

	int miss_align = offset%bs;

	i = 0;
	if(miss_align==0) // B is aligned
		{
		for(; i<row-3; i+=4)
			{
			ptrA = A + (i+offset)*sda;
			ptrB = B + i*sdb;
			j = 0;
			for(; j<col-3; j+=4)
				{
				ptrB[0+(j+0)*bs] = ptrA[0+(j+0)*bs];
				ptrB[1+(j+0)*bs] = ptrA[1+(j+0)*bs];
				ptrB[2+(j+0)*bs] = ptrA[2+(j+0)*bs];
				ptrB[3+(j+0)*bs] = ptrA[3+(j+0)*bs];

				ptrB[0+(j+1)*bs] = ptrA[0+(j+1)*bs];
				ptrB[1+(j+1)*bs] = ptrA[1+(j+1)*bs];
				ptrB[2+(j+1)*bs] = ptrA[2+(j+1)*bs];
				ptrB[3+(j+1)*bs] = ptrA[3+(j+1)*bs];

				ptrB[0+(j+2)*bs] = ptrA[0+(j+2)*bs];
				ptrB[1+(j+2)*bs] = ptrA[1+(j+2)*bs];
				ptrB[2+(j+2)*bs] = ptrA[2+(j+2)*bs];
				ptrB[3+(j+2)*bs] = ptrA[3+(j+2)*bs];

				ptrB[0+(j+3)*bs] = ptrA[0+(j+3)*bs];
				ptrB[1+(j+3)*bs] = ptrA[1+(j+3)*bs];
				ptrB[2+(j+3)*bs] = ptrA[2+(j+3)*bs];
				ptrB[3+(j+3)*bs] = ptrA[3+(j+3)*bs];
				}
			for(; j<col; j++)
				{
				ptrB[0+j*bs] = ptrA[0+j*bs];
				ptrB[1+j*bs] = ptrA[1+j*bs];
				ptrB[2+j*bs] = ptrA[2+j*bs];
				ptrB[3+j*bs] = ptrA[3+j*bs];
				}
			}
		}
	for(; i<row; i++)
		{
		ptrA = A + ((offset+i)/bs)*bs*sda + ((offset+i)%bs);
		ptrB = B + (i/bs)*bs*sdb + (i%bs);
		j = 0;
		for(; j<col-3; j+=4)
			{
			ptrB[(j+0)*bs] = ptrA[(j+0)*bs];
			ptrB[(j+1)*bs] = ptrA[(j+1)*bs];
			ptrB[(j+2)*bs] = ptrA[(j+2)*bs];
			ptrB[(j+3)*bs] = ptrA[(j+3)*bs];
			}
		for(; j<col; j++)
			{
			ptrB[j*bs] = ptrA[j*bs];
			}
		}

	}



#if ! defined(BLASFEO)
/* converts a matrix into a packed matrix */
void d_cvt_mat2pmat(int row, int col, double *A, int lda, int offset, double *pA, int sda)
	{

	const int bs = 4;

	int
		i, ii, j, jj, row0, row1, row2;

	double
		*B, *pB;

#if defined(TARGET_X64_AVX2) || defined(TARGET_X64_AVX)
	__m256d
		tmp;
#endif

	row0 = (bs-offset%bs)%bs;
	if(row0>row)
		row0 = row;
	row1 = row - row0;

#if 1
	jj = 0;
	for( ; jj<col-3; jj+=4)
		{

		B  =  A + jj*lda;
		pB = pA + jj*bs;

		ii = 0;
		if(row0>0)
			{
			for( ; ii<row0; ii++)
				{
				pB[ii+bs*0] = B[ii+lda*0];
				pB[ii+bs*1] = B[ii+lda*1];
				pB[ii+bs*2] = B[ii+lda*2];
				pB[ii+bs*3] = B[ii+lda*3];
				}
			B  += row0;
			pB += row0 + bs*(sda-1);
			}
#if defined(TARGET_X64_AVX2) || defined(TARGET_X64_AVX)
		for( ; ii<row-3; ii+=4)
			{
			tmp = _mm256_loadu_pd( &B[0+lda*0] );
			_mm256_store_pd( &pB[0+bs*0], tmp );
			tmp = _mm256_loadu_pd( &B[0+lda*1] );
			_mm256_store_pd( &pB[0+bs*1], tmp );
			tmp = _mm256_loadu_pd( &B[0+lda*2] );
			_mm256_store_pd( &pB[0+bs*2], tmp );
			tmp = _mm256_loadu_pd( &B[0+lda*3] );
			_mm256_store_pd( &pB[0+bs*3], tmp );
			// update
			B  += 4;
			pB += bs*sda;
			}
#else
		for( ; ii<row-3; ii+=4)
			{
			// col 0
			pB[0+bs*0] = B[0+lda*0];
			pB[1+bs*0] = B[1+lda*0];
			pB[2+bs*0] = B[2+lda*0];
			pB[3+bs*0] = B[3+lda*0];
			// col 1
			pB[0+bs*1] = B[0+lda*1];
			pB[1+bs*1] = B[1+lda*1];
			pB[2+bs*1] = B[2+lda*1];
			pB[3+bs*1] = B[3+lda*1];
			// col 2
			pB[0+bs*2] = B[0+lda*2];
			pB[1+bs*2] = B[1+lda*2];
			pB[2+bs*2] = B[2+lda*2];
			pB[3+bs*2] = B[3+lda*2];
			// col 3
			pB[0+bs*3] = B[0+lda*3];
			pB[1+bs*3] = B[1+lda*3];
			pB[2+bs*3] = B[2+lda*3];
			pB[3+bs*3] = B[3+lda*3];
			// update
			B  += 4;
			pB += bs*sda;
			}
#endif
		for( ; ii<row; ii++)
			{
			// col 0
			pB[0+bs*0] = B[0+lda*0];
			// col 1
			pB[0+bs*1] = B[0+lda*1];
			// col 2
			pB[0+bs*2] = B[0+lda*2];
			// col 3
			pB[0+bs*3] = B[0+lda*3];
			// update
			B  += 1;
			pB += 1;
			}
		}
	for( ; jj<col; jj++)
		{

		B  =  A + jj*lda;
		pB = pA + jj*bs;

		ii = 0;
		if(row0>0)
			{
			for( ; ii<row0; ii++)
				{
				pB[ii+bs*0] = B[ii+lda*0];
				}
			B  += row0;
			pB += row0 + bs*(sda-1);
			}
		for( ; ii<row-3; ii+=4)
			{
			// col 0
			pB[0+bs*0] = B[0+lda*0];
			pB[1+bs*0] = B[1+lda*0];
			pB[2+bs*0] = B[2+lda*0];
			pB[3+bs*0] = B[3+lda*0];
			// update
			B  += 4;
			pB += bs*sda;
			}
		for( ; ii<row; ii++)
			{
			// col 0
			pB[0+bs*0] = B[0+lda*0];
			// update
			B  += 1;
			pB += 1;
			}
		}

#else
	ii = 0;
	if(row0>0)
		{
		for(j=0; j<col; j++)
			{
			for(i=0; i<row0; i++)
				{
				pA[i+j*bs+ii*sda] = A[i+ii+j*lda];
				}
			}

		A  += row0;
		pA += row0 + bs*(sda-1);
		}

	ii = 0;
	for(; ii<row1-3; ii+=bs)
		{
		j=0;
		for(; j<col-3; j+=4)
			{
#if defined(TARGET_X64_AVX)
			tmp = _mm256_loadu_pd( &A[ii+0+(j+0)*lda] );
			_mm256_store_pd( &pA[0+(j+0)*bs+ii*sda], tmp );
			tmp = _mm256_loadu_pd( &A[ii+0+(j+1)*lda] );
			_mm256_store_pd( &pA[0+(j+1)*bs+ii*sda], tmp );
			tmp = _mm256_loadu_pd( &A[ii+0+(j+2)*lda] );
			_mm256_store_pd( &pA[0+(j+2)*bs+ii*sda], tmp );
			tmp = _mm256_loadu_pd( &A[ii+0+(j+3)*lda] );
			_mm256_store_pd( &pA[0+(j+3)*bs+ii*sda], tmp );
#else
			// unroll 0
			pA[0+(j+0)*bs+ii*sda] = A[ii+0+(j+0)*lda];
			pA[1+(j+0)*bs+ii*sda] = A[ii+1+(j+0)*lda];
			pA[2+(j+0)*bs+ii*sda] = A[ii+2+(j+0)*lda];
			pA[3+(j+0)*bs+ii*sda] = A[ii+3+(j+0)*lda];
			// unroll 1
			pA[0+(j+1)*bs+ii*sda] = A[ii+0+(j+1)*lda];
			pA[1+(j+1)*bs+ii*sda] = A[ii+1+(j+1)*lda];
			pA[2+(j+1)*bs+ii*sda] = A[ii+2+(j+1)*lda];
			pA[3+(j+1)*bs+ii*sda] = A[ii+3+(j+1)*lda];
			// unroll 2
			pA[0+(j+2)*bs+ii*sda] = A[ii+0+(j+2)*lda];
			pA[1+(j+2)*bs+ii*sda] = A[ii+1+(j+2)*lda];
			pA[2+(j+2)*bs+ii*sda] = A[ii+2+(j+2)*lda];
			pA[3+(j+2)*bs+ii*sda] = A[ii+3+(j+2)*lda];
			// unroll 3
			pA[0+(j+3)*bs+ii*sda] = A[ii+0+(j+3)*lda];
			pA[1+(j+3)*bs+ii*sda] = A[ii+1+(j+3)*lda];
			pA[2+(j+3)*bs+ii*sda] = A[ii+2+(j+3)*lda];
			pA[3+(j+3)*bs+ii*sda] = A[ii+3+(j+3)*lda];
#endif
			}
		for(; j<col; j++)
			{
			pA[0+(j+0)*bs+ii*sda] = A[ii+0+(j+0)*lda];
			pA[1+(j+0)*bs+ii*sda] = A[ii+1+(j+0)*lda];
			pA[2+(j+0)*bs+ii*sda] = A[ii+2+(j+0)*lda];
			pA[3+(j+0)*bs+ii*sda] = A[ii+3+(j+0)*lda];
			}
		}
	if(ii<row1)
		{
		row2 = row1-ii;
		if(bs<row2) row2 = bs;
		for(j=0; j<col; j++)
			{
			for(i=0; i<row2; i++)
				{
				pA[i+j*bs+ii*sda] = A[i+ii+j*lda];
				}
			}
		}
#endif

	}
#endif



#if ! defined(BLASFEO)
/* converts a matrix into a packed matrix */
// row and col of the source matrix, offsett in the destination matrix
void d_cvt_tran_mat2pmat(int row, int col, double *A, int lda, int offset, double *pA, int sda)
	{

	const int bs = 4;

	int i, ii, j, row0, row1, row2;

	double
		*B, *pB;

#if defined(TARGET_X64_AVX2) || defined(TARGET_X64_AVX)
	__m256d
		v0, v1, v2, v3,
		v4, v5, v6, v7;
#endif

	row0 = (bs-offset%bs)%bs;
	if(row0>col)
		row0 = col;
	row1 = col - row0;

	ii = 0;
	if(row0>0)
		{
		for(j=0; j<row; j++)
			{
			for(i=0; i<row0; i++)
				{
				pA[i+j*bs+ii*sda] = A[j+(i+ii)*lda];
				}
			}

		A  += row0*lda;
		pA += row0 + bs*(sda-1);
		}

	ii = 0;
	for(; ii<row1-3; ii+=bs)
		{
		j=0;
		B  = A + ii*lda;
		pB = pA + ii*sda;
#if defined(TARGET_X64_AVX2) || defined(TARGET_X64_AVX)
		for(; j<row-3; j+=4)
			{
			v0 = _mm256_loadu_pd( &B[0+0*lda] ); // 00 10 20 30
			v1 = _mm256_loadu_pd( &B[0+1*lda] ); // 01 11 21 31
			v4 = _mm256_unpacklo_pd( v0, v1 ); // 00 01 20 21
			v5 = _mm256_unpackhi_pd( v0, v1 ); // 10 11 30 31
			v2 = _mm256_loadu_pd( &B[0+2*lda] ); // 02 12 22 32
			v3 = _mm256_loadu_pd( &B[0+3*lda] ); // 03 13 23 33
			v6 = _mm256_unpacklo_pd( v2, v3 ); // 02 03 22 23
			v7 = _mm256_unpackhi_pd( v2, v3 ); // 12 13 32 33

			B += 4;

			v0 = _mm256_permute2f128_pd( v4, v6, 0x20 ); // 00 01 02 03
			_mm256_store_pd( &pB[0+bs*0], v0 );
			v2 = _mm256_permute2f128_pd( v4, v6, 0x31 ); // 20 21 22 23
			_mm256_store_pd( &pB[0+bs*2], v2 );
			v1 = _mm256_permute2f128_pd( v5, v7, 0x20 ); // 10 11 12 13
			_mm256_store_pd( &pB[0+bs*1], v1 );
			v3 = _mm256_permute2f128_pd( v5, v7, 0x31 ); // 30 31 32 33
			_mm256_store_pd( &pB[0+bs*3], v3 );

			pB += 4*bs;
			}
#else
		for(; j<row-3; j+=4)
			{
#if 1
			// unroll 0
			pB[0+0*bs] = B[0+0*lda];
			pB[1+0*bs] = B[0+1*lda];
			pB[2+0*bs] = B[0+2*lda];
			pB[3+0*bs] = B[0+3*lda];
			// unroll 1
			pB[0+1*bs] = B[1+0*lda];
			pB[1+1*bs] = B[1+1*lda];
			pB[2+1*bs] = B[1+2*lda];
			pB[3+1*bs] = B[1+3*lda];
			// unroll 2
			pB[0+2*bs] = B[2+0*lda];
			pB[1+2*bs] = B[2+1*lda];
			pB[2+2*bs] = B[2+2*lda];
			pB[3+2*bs] = B[2+3*lda];
			// unroll 3
			pB[0+3*bs] = B[3+0*lda];
			pB[1+3*bs] = B[3+1*lda];
			pB[2+3*bs] = B[3+2*lda];
			pB[3+3*bs] = B[3+3*lda];
			B  += 4;
			pB += 4*bs;
#else
			// unroll 0
			pA[0+(j+0)*bs+ii*sda] = A[j+0+(ii+0)*lda];
			pA[1+(j+0)*bs+ii*sda] = A[j+0+(ii+1)*lda];
			pA[2+(j+0)*bs+ii*sda] = A[j+0+(ii+2)*lda];
			pA[3+(j+0)*bs+ii*sda] = A[j+0+(ii+3)*lda];
			// unroll 1
			pA[0+(j+1)*bs+ii*sda] = A[j+1+(ii+0)*lda];
			pA[1+(j+1)*bs+ii*sda] = A[j+1+(ii+1)*lda];
			pA[2+(j+1)*bs+ii*sda] = A[j+1+(ii+2)*lda];
			pA[3+(j+1)*bs+ii*sda] = A[j+1+(ii+3)*lda];
			// unroll 2
			pA[0+(j+2)*bs+ii*sda] = A[j+2+(ii+0)*lda];
			pA[1+(j+2)*bs+ii*sda] = A[j+2+(ii+1)*lda];
			pA[2+(j+2)*bs+ii*sda] = A[j+2+(ii+2)*lda];
			pA[3+(j+2)*bs+ii*sda] = A[j+2+(ii+3)*lda];
			// unroll 3
			pA[0+(j+3)*bs+ii*sda] = A[j+3+(ii+0)*lda];
			pA[1+(j+3)*bs+ii*sda] = A[j+3+(ii+1)*lda];
			pA[2+(j+3)*bs+ii*sda] = A[j+3+(ii+2)*lda];
			pA[3+(j+3)*bs+ii*sda] = A[j+3+(ii+3)*lda];
#endif
			}
#endif
		for(; j<row; j++)
			{
#if 1
			// unroll 0
			pB[0+0*bs] = B[0+0*lda];
			pB[1+0*bs] = B[0+1*lda];
			pB[2+0*bs] = B[0+2*lda];
			pB[3+0*bs] = B[0+3*lda];
			B  += 1;
			pB += 1*bs;
#else
			pA[0+(j+0)*bs+ii*sda] = A[j+0+(ii+0)*lda];
			pA[1+(j+0)*bs+ii*sda] = A[j+0+(ii+1)*lda];
			pA[2+(j+0)*bs+ii*sda] = A[j+0+(ii+2)*lda];
			pA[3+(j+0)*bs+ii*sda] = A[j+0+(ii+3)*lda];
#endif
			}
		}
	if(ii<row1)
		{
		row2 = row1-ii;
		if(bs<row2) row2 = bs;
		for(j=0; j<row; j++)
			{
			for(i=0; i<row2; i++)
				{
				pA[i+j*bs+ii*sda] = A[j+(i+ii)*lda];
				}
			}
		}

	}
#endif



#if ! defined(BLASFEO)
/* converts a packed matrix into a matrix */
void d_cvt_pmat2mat(int row, int col, int offset, double *pA, int sda, double *A, int lda)
	{

	const int bs = 4;

	int i, ii, jj;

	int row0 = (bs-offset%bs)%bs;

	double *ptr_pA;


	jj=0;
	for(; jj<col; jj++)
		{
		ptr_pA = pA + jj*bs;
		ii = 0;
		if(row0>0)
			{
			for(; ii<row0; ii++)
				{
				A[ii+lda*jj] = ptr_pA[0];
				ptr_pA++;
				}
			ptr_pA += (sda-1)*bs;
			}
		for(; ii<row-bs+1; ii+=bs)
			{
			i=0;
			for(; i<bs; i++)
				{
				A[i+ii+lda*jj] = ptr_pA[0];
				ptr_pA++;
				}
			ptr_pA += (sda-1)*bs;
			}
		for(; ii<row; ii++)
			{
			A[ii+lda*jj] = ptr_pA[0];
			ptr_pA++;
			}
		}

	}
#endif



#if ! defined(BLASFEO)
/* converts a packed matrix into a matrix */
void d_cvt_tran_pmat2mat(int row, int col, int offset, double *pA, int sda, double *A, int lda)
	{

	const int bs = 4;

	int i, ii, jj;

	int row0 = (bs-offset%bs)%bs;

	double *ptr_pA;


	jj=0;
	for(; jj<col; jj++)
		{
		ptr_pA = pA + jj*bs;
		ii = 0;
		if(row0>0)
			{
			for(; ii<row0; ii++)
				{
				A[jj+lda*ii] = ptr_pA[0];
				ptr_pA++;
				}
			ptr_pA += (sda-1)*bs;
			}
		for(; ii<row-bs+1; ii+=bs)
			{
			i=0;
			for(; i<bs; i++)
				{
				A[jj+lda*(i+ii)] = ptr_pA[0];
				ptr_pA++;
				}
			ptr_pA += (sda-1)*bs;
			}
		for(; ii<row; ii++)
			{
			A[jj+lda*ii] = ptr_pA[0];
			ptr_pA++;
			}
		}

	}
#endif
