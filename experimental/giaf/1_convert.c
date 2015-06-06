/******************************************************************************
*
* Example of the use of some LA routines in HPMPC
*
* Author: Gianluca Frison
*
******************************************************************************/

#include <stdlib.h>
#include <stdio.h>

// include HPMPC headers
#include <hpmpc/target.h>
#include <hpmpc/block_size.h>
#include <hpmpc/aux_d.h>
#include <hpmpc/blas_d.h>



int main()
	{

	// architecture-dependent quantities
	const int bs = D_MR; // number of rows in a panel
	const int ncl = D_NCL; // number of columns to make the following panel cache-aligned



	// loop index
	int ii, jj;



	// matrix size
	int n = 10;



	// create data matrices in column-major format (not-aligned memory allocation)
	double *A; d_zeros(&A, n, n);
	double *B; d_zeros(&B, n, n);
	double *C; d_zeros(&C, n, n);
	double *D; d_zeros(&D, n, n);

	for(ii=0; ii<n*n; ii++)
		A[ii] = ii;
	
	for(ii=0; ii<n; ii++)
		B[ii*(n+1)] = 2.0;
	
	// print matrices in column-major format
	// void d_print_mat(int rows, int cols, double *A, int lda)
	// lda is Leading Dimension of matrix A, that is the size in ROWS of the matrix containing A
	printf("\nprint matrix in column-major order\n");
	d_print_mat(n, n, A, n);
	d_print_mat(n, n, B, n);



	// number of rows in the packed matrix format
	int pn = (n+bs-1)/bs*bs;
	// number of columns in the packed matrix format
	int cn = (n+ncl-1)/ncl*ncl;

	// create matrices in panel-wise format (aligned memory allocation)
	double *pA; d_zeros_align(&pA, pn, cn);
	double *pAt; d_zeros_align(&pAt, pn, cn);
	double *pB; d_zeros_align(&pB, pn, cn);
	double *pC; d_zeros_align(&pC, pn, cn);

	// convert from column-major to panel wise
	// conver A
	d_cvt_mat2pmat(n, n, A, n, 0, pA, cn);
	// conver & transpose A
	d_cvt_tran_mat2pmat(n, n, A, n, 0, pAt, cn);
	// conver B
	d_cvt_mat2pmat(n, n, B, n, 0, pB, cn);

	// print matrices in panel-wise format
	// void d_print_pmat(int rows, int cols, int bs, double *pA, int sda)
	// sda is Second Dimension of matrix pA, that is the size in COLUMNS for the pmatrix (panel-wise matrix) containint A
	printf("\nprint matrix in panel-wise order\n");
	d_print_pmat(n, n, bs, pA, cn);
	d_print_pmat(n, n, bs, pAt, cn);
	d_print_pmat(n, n, bs, pB, cn);

	printf("\nprint matrix in panel-wise order, one panel at a time\n");
	for(ii=0; ii<n; ii+=bs)
		d_print_mat(bs, cn, pA+ii*cn, bs);
	

	// convert from panel-wise to column-major
	d_cvt_pmat2mat(n, n, 0, pA, cn, C, n);
	d_cvt_tran_pmat2mat(n, n, 0, pA, cn, D, n);
	printf("\nprint matrix in column-major order\n");
	d_print_mat(n, n, C, n);
	d_print_mat(n, n, D, n);



	// free memory
	free(A);
	free(B);
	free(C);
	free(D);
	free(pA);
	free(pAt);
	free(pB);
	free(pC);

	// return
	return 0;

	}
