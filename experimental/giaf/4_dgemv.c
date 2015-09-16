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

	for(ii=0; ii<n*n; ii++)
		A[ii] = ii;
	

	// number of rows in the packed matrix format
	int pn = (n+bs-1)/bs*bs;
	// number of columns in the packed matrix format
	int cn = (n+ncl-1)/ncl*ncl;

	// create matrices in panel-wise format (aligned memory allocation)
	double *pA; d_zeros_align(&pA, pn, cn);

	// convert from column-major to panel wise (sub-matrices)
	d_cvt_mat2pmat(n, n, A, n, 0, pA, cn);

	// print pmat
	d_print_pmat(n, n, bs, pA, cn);

	// create vectors 
	double *x; d_zeros_align(&x, pn, 1); // x has to be aligned to SIMD boundaries
	double *y; d_zeros(&y, pn, 1);
	double *z; d_zeros(&z, pn, 1);

	for(ii=0; ii<n; ii++) x[ii] = 1;
	for(ii=0; ii<n; ii++) y[ii] = 1e3;

	// print vector
	d_print_mat(n, 1, x, n);
	d_print_mat(n, 1, y, n);

	// void dgemv_n_lib(int m, int n, double *pA, int sda, double *x, int alg, double *y, double *z)
	// void dgemv_t_lib(int m, int n, double *pA, int sda, double *x, int alg, double *y, double *z)
	//
	// Vector elements are supposed to be stored contiguously (i.e. inc==1 in standard BLAS). Therefore it is not possible to use there routines to operate on columns or rows of matrices.
	// There are 1 matrix and 3 vector arguments (that may be useful to avoid a matrix copy), and there are only 3 combinations of alpha and beta (controlled by 'alg'):
	// if alg==0  : z = pA * x
	// if alg==1  : z = y + pA * x
	// if alg==-1 : z = y - pA * x
	// In my code I never had to use other alpha and beta combinations, but in case of need they can be easily added.
	//
	// The sizes m and n refer to the size of the matrix pA (so in the transposed version, pA' is n times m)
	//
	// The matrix pA have to be at the top of a panel (it is not possible to have them at arbitrary vertical position in a matrix).
	// Changes to this would require a quite some programming work and performance loss
	// The vector x has to be aligned to SIMD boundaries (this is assumed for performance purposes); y and z can be arbitrarly aligned
	//
	dgemv_n_lib(n, n, pA, cn, x, 1, y, z);
	d_print_mat(n, 1, z, n);

	dgemv_t_lib(n, n, pA, cn, x, 1, y, z);
	d_print_mat(n, 1, z, n);



	// free memory
	free(A);
	free(pA);
	free(x);
	free(y);
	free(z);

	// return
	return 0;

	}
