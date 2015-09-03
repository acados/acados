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

	for(ii=0; ii<n*n; ii++)
		A[ii] = ii;
	
	for(ii=0; ii<n; ii++)
		B[ii*(n+1)] = 2.0;
	

	// number of rows in the packed matrix format
	int pn = (n+bs-1)/bs*bs;
	// number of columns in the packed matrix format
	int cn = (n+ncl-1)/ncl*ncl;

	// create matrices in panel-wise format (aligned memory allocation)
	double *pA; d_zeros_align(&pA, pn, cn);
	double *pB; d_zeros_align(&pB, pn, cn);
	double *pC; d_zeros_align(&pC, pn, cn);
	double *pD; d_zeros_align(&pD, pn, cn);

	// convert from column-major to panel wise (sub-matrices)
	d_cvt_mat2pmat(n, n, A, n, 0, pA, cn);
	d_cvt_mat2pmat(n, n, B, n, 0, pB, cn);

	// print pmat
	d_print_pmat(n, n, bs, pA, cn);
	d_print_pmat(n, n, bs, pB, cn);

	// void dgemm_nt_lib(int m, int n, int k, double *pA, int sda, double *pB, int sdb, int alg, double *pC, int sdc, double *pD, int sdd, int tc, int td)
	//
	// Compared to the dgemm in BLAS, here the pA matrix is always normal and the pB matrix is always transposed: this gives optimal performance.
	// There are 4 matrix arguments (that may be useful to avoid a matrix copy), and there are only 3 combinations of alpha and beta (controlled by 'alg'):
	// if alg==0  : pD = pA * pB'
	// if alg==1  : pD = pC + pA * pB'
	// if alg==-1 : pD = pC - pA * pB'
	// In my code I never had to use other alpha and beta combinations, but in case of need they can be easily added.
	//
	// Furthermore tc and td tells if the pC and pD matrices are supposed to be transposed, as follows (example in the case of alg==1):
	// if tc==0 & td==0 : pD =  pC  + pA * pB'
	// if tc==0 & td==1 : pD = (pC  + pA * pB')'
	// if tc==1 & td==0 : pD =  pC' + pA * pB'
	// if tc==1 & td==1 : pD = (pC' + pA * pB')'
	// This can avoid the need to explicitly transpose matrices, and to permute pA and pB to have the best performance (if m, n and k are not of the same size, the best performane is obtained for k > m > n).
	//
	// The sizes m and n refer to the size of the product (pA * pB'), so if e.g. tc==0 & td==1, pC is of size m times n, while pD is of size n times m
	//
	// All matrices pA, pB, pC, pD have to be at the top of a panel (it is not possible to have them at arbitrary vertical position in a matrix).
	// Changes to this would require a lot of programming work and performance loss: therefore the best option is probably the combination of dgemm + dgecp.
	//
	for(ii=0; ii<cn*pn; ii++) pD[ii] = 0.0;
	dgemm_nt_lib(5, 7, n, pA, cn, pB, cn, 0, pC, cn, pD+1*bs, cn, 0, 0);
	d_print_pmat(n, n, bs, pD, cn);

	for(ii=0; ii<cn*pn; ii++) pD[ii] = 0.0;
	dgemm_nt_lib(5, 7, n, pA, cn, pB, cn, 0, pC, cn, pD+1*bs, cn, 0, 1);
	d_print_pmat(n, n, bs, pD, cn);

	for(ii=0; ii<cn*pn; ii++) pD[ii] = 0.0;
	dgemm_nt_lib(5, 7, n, pB, cn, pA, cn, -1, pC, cn, pD+1*bs, cn, 0, 0);
	d_print_pmat(n, n, bs, pD, cn);


	// void dgemm_nn_lib(int m, int n, int k, double *pA, int sda, double *pB, int sdb, double *pC, int sdc, double *pD, int sdd, int alg, int tc, int td)
	//
	// This routine is identical to the dgemm_nt, with the difference that both pA and pB are normal (i.e. not transposed); however, by its nature this routine has worse performance, so dgemm_nt should be used if possible.
	// Furthermore, a.t.m. this routine is not well optimized.
	//
	// XXX Pull the latest HPMPC version !!!
	//
	for(ii=0; ii<cn*pn; ii++) pD[ii] = 0.0;
	dgemm_nn_lib(5, 7, n, pB, cn, pA, cn, -1, pC, cn, pD+1*bs, cn, 0, 0);
	d_print_pmat(n, n, bs, pD, cn);



	// free memory
	free(A);
	free(B);
	free(pA);
	free(pB);
	free(pC);
	free(pD);

	// return
	return 0;

	}
