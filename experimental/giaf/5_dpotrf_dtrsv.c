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

	for(ii=0; ii<n; ii++)
		A[ii+n*ii] = 2.0;
	

	// number of rows in the packed matrix format
	int pn = (n+bs-1)/bs*bs;
	// number of columns in the packed matrix format
	int cn = (n+ncl-1)/ncl*ncl;

	// create matrices in panel-wise format (aligned memory allocation)
	double *pA; d_zeros_align(&pA, pn, cn);
	double *pB; d_zeros_align(&pB, pn, cn);

	// convert from column-major to panel wise (sub-matrices)
	d_cvt_mat2pmat(n, n, A, n, 0, pA, cn);

	// print pmat
	d_print_pmat(n, n, bs, pA, cn);

	// create vector holding the inverse of diagonal of B
	double *diag; d_zeros(&diag, n, 1);

	// create vectors 
	double *x; d_zeros_align(&x, pn, 1); // x has to be aligned to SIMD boundaries
	double *y; d_zeros_align(&y, pn, 1); // y has to be aligned to SIMD boundaries

	for(ii=0; ii<n; ii++) x[ii] = 1;

	// print vector
	d_print_mat(n, 1, x, n);

	// void dpotrf_lib(int m, int n, double *pA, ind sda, double *pB, int sdb, double *diag);
	//
	// If m==n, it computes in pB the lower triangular Cholesky factor of the n times n matrix pA
	// pB = chol(pA, 'lower')
	//
	// If m>n, it computes in the top n timex n sub-matrix of pB the lower Cholesky factor of the top sub-matrix of pA;
	// and it computes in the bottom (m-n) times n sub-matrix of pB the substitution with matrix RHS using the newly computed lower Cholesky factor of the top sub-matrix of pA
	// pB(1:n,1:n) = chol(pA(1:n,1:n),'lower')
	// pB(n+1:m,1:n) = pA(n+1:m,1:n) * inv(pB(1:n,1:n)')
	//

	// void dtrsv_n_lib(int m, int n, double *pA, int sda, int use_inv_diag_A, double *inv_diag_A, double *x, double *y)
	// void dtrsv_t_lib(int m, int n, double *pA, int sda, int use_inv_diag_A, double *inv_diag_A, double *x, double *y)
	//
	// If m==n, it computes the substitution with vector RHS
	// y = inv(A) * x     ('n' version)
	// y = inv(A') * x    ('t' version)
	//
	// If m>n, if A, x and y are decomposed as
	// A = [A0; A1];
	// x = [x0; x1];
	// y = [y0; y1];
	// it computes ('n' version)
	// [y0; y1] = [inv(A0)*x0; x1-A1*y0];
	// and ('t' version)
	// [y0; y1] = [inv(A0')*(x0-A1'*x1; x1];
	//
	// The vector inv_diag_A can be used to provide the inverse of the diagonal of A (as e.g. returned by the dpotrf_lib routine), avoiding the divisions
	//
	d_set_pmat(n, n, 0.0, 0, pB, cn); // set matrix to 0
	dpotrf_lib(n, n, pA, cn, pB, cn, diag);
	d_print_pmat(n, n, bs, pA, cn);
	d_print_pmat(n, n, bs, pB, cn);

	double *dummy;
	dtrsv_n_lib(n, n, pB, cn, 0, dummy, x, y); // the routine uses divisions
	d_print_mat(n, 1, y, n);

	dtrsv_n_lib(n, n, pB, cn, 1, diag, x, y); // the routine does not use divisions
	d_print_mat(n, 1, y, n);



	// free memory
	free(A);
	free(pA);
	free(pB);
	free(diag);
	free(x);
	free(y);

	// return
	return 0;

	}
