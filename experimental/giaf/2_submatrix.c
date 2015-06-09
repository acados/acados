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
	int p2n = (2*n+bs-1)/bs*bs;
	// number of columns in the packed matrix format
	int cn = (n+ncl-1)/ncl*ncl;
	int c2n = (2*n+ncl-1)/ncl*ncl;

	// create matrices in panel-wise format (aligned memory allocation)
	double *pA; d_zeros_align(&pA, pn, cn);
	double *pAt; d_zeros_align(&pAt, pn, cn);
	double *pB; d_zeros_align(&pB, pn, cn);
	double *pC; d_zeros_align(&pC, pn, cn);
	double *pD; d_zeros_align(&pD, pn, cn);
	double *pAB; d_zeros_align(&pAB, p2n, c2n);

	// convert from column-major to panel wise (sub-matrices)
	// d_cvt_mat2pmat( int rows, int cols, double *A, int lda, int offset, double *pA, int sda);
	// offset gives the position in the panel of the first element of the matrix
	// conver A in the top-left corner of pAB: offset=0
	d_cvt_mat2pmat(n, n, A, n, 0, pAB, c2n);
	// conver B in the bottom-right corner of pAB: offset=n (i.e. n rows from the top of the matrix)
	// B has to be copied at row n: this means panel number n/bs _times_ bs*c2n size of a panel _plus_ n%bs position in the panel
	// B has to be copied at col n: this means _plus_ n*bs (n cols * bs leading dimension of each panel)
	d_cvt_mat2pmat(n, n, B, n, n, pAB+n/bs*bs*c2n+n%bs+n*bs, c2n);

	// print matrices in panel-wise format
	// void d_print_pmat(int rows, int cols, int bs, double *pA, int sda)
	// sda is Second Dimension of matrix pA, that is the size in COLUMNS for the pmatrix (panel-wise matrix) containint A
	printf("\nprint matrix in panel-wise order\n");
	d_print_pmat(2*n, 2*n, bs, pAB, c2n);

	printf("\nprint matrix in panel-wise order, one panel at a time\n");
	for(ii=0; ii<2*n; ii+=bs)
		d_print_mat(bs, c2n, pAB+ii*c2n, bs);
	
	

	// copy a submatrix into a submatrix
	// pC(1:5,3:9) = pAB(2:6,1:7)
	// dgecp_lib( int rows, int cols, int offsetA, double *pA, int sda, int offsetC, double *C, int sdc);
	dgecp_lib(4, 6, 2, pAB+2/bs*bs*c2n+2%bs+1*bs, c2n, 1, pC+1/bs*bs*cn+1%bs+3*bs, cn);
	printf("\nprint matrix in panel-wise orderr\n");
	d_print_pmat(n, n, bs, pC, cn);

	printf("\nprint matrix in panel-wise order, one panel at a time\n");
	for(ii=0; ii<n; ii+=bs)
		d_print_mat(bs, cn, pC+ii*cn, bs);


	// transpose a submatrix into a submatrix
	// pC(1:7,3:7) = (pAB(2:6,1:7))'
	// dgetr_lib( int rows, int cols, int offsetA, double *pA, int sda, int offsetC, double *C, int sdc);
	// rows and cols are referred to pA
	dgetr_lib(4, 6, 2, pAB+2/bs*bs*c2n+2%bs+1*bs, c2n, 1, pD+1/bs*bs*cn+1%bs+3*bs, cn);
	printf("\nprint matrix in panel-wise orderr\n");
	d_print_pmat(n, n, bs, pD, cn);

	printf("\nprint matrix in panel-wise order, one panel at a time\n");
	for(ii=0; ii<n; ii+=bs)
		d_print_mat(bs, cn, pD+ii*cn, bs);


	// free memory
	free(A);
	free(B);
	free(C);
	free(D);
	free(pA);
	free(pAt);
	free(pB);
	free(pAB);
	free(pC);
	free(pD);

	// return
	return 0;

	}
