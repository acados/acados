/*
 * Copyright (c) The acados authors.
 *
 * This file is part of acados.
 *
 * Licensed under the 2-Clause BSD License.
 */

#ifndef ACADOS_UTILS_MATH_H_
#define ACADOS_UTILS_MATH_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "acados/utils/types.h"
#include "blasfeo_common.h"

#if defined(__MABX2__)
double fmax(double a, double b);
int isnan(double x);
#endif
#if defined(_DS1104)
double fmax(double a, double b);
int isnan(double x);
#endif

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

void dgemm_nn_3l(int m, int n, int k, double *A, int lda, double *B, int ldb, double *C, int ldc);
// void dgemv_n_3l(int m, int n, double *A, int lda, double *x, double *y);
// void dgemv_t_3l(int m, int n, double *A, int lda, double *x, double *y);
// void dcopy_3l(int n, double *x, int incx, double *y, int incy);
void daxpy_3l(int n, double da, double *dx, double *dy);
void dscal_3l(int n, double da, double *dx);
double twonormv(int n, double *ptrv);

/* copies a matrix into another matrix */
void dmcopy(int row, int col, double *ptrA, int lda, double *ptrB, int ldb);

/* solution of a system of linear equations */
void dgesv_3l(int n, int nrhs, double *A, int lda, int *ipiv, double *B, int ldb, int *info);

/* matrix exponential */
void expm(int row, double *A);

int idamax_3l(int n, double *x);

void dswap_3l(int n, double *x, int incx, double *y, int incy);

void dger_3l(int m, int n, double alpha, double *x, int incx, double *y, int incy, double *A,
             int lda);

void dgetf2_3l(int m, int n, double *A, int lda, int *ipiv, int *info);

void dlaswp_3l(int n, double *A, int lda, int k1, int k2, int *ipiv);

void dtrsm_l_l_n_u_3l(int m, int n, double *A, int lda, double *B, int ldb);

void dgetrs_3l(int n, int nrhs, double *A, int lda, int *ipiv, double *B, int ldb);

void dgesv_3l(int n, int nrhs, double *A, int lda, int *ipiv, double *B, int ldb, int *info);

double onenorm(int row, int col, double *ptrA);

// double twonormv(int n, double *ptrv);

void padeapprox(int m, int row, double *A);

void expm(int row, double *A);

// void d_compute_qp_size_ocp2dense_rev(int N, int *nx, int *nu, int *nb, int **hidxb, int *ng,
//                                      int *nvd, int *ned, int *nbd, int *ngd);

// Eigendecomposition of matrix A: reads only the lower triangular of A
void acados_eigen_decomposition(int dim, double *A, double *V, double *d, double *e);

double minimum_of_doubles(double *x, int n);

void neville_algorithm(double xx, int n, double *x, double *Q, double *out);

// assumes A symmetric, stored as lower triangular
void compute_gershgorin_max_abs_eig_estimate(int n, struct blasfeo_dmat *A, double *out);

// assumes A symmetric, stored as lower triangular
void compute_gershgorin_min_eig_estimate(int n, struct blasfeo_dmat *A, double *out);


#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_UTILS_MATH_H_
