/**************************************************************************************************
* acados/acados/utils/tools.h                                                                                                *
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
#ifndef ACADOS_UTILS_TOOLS_H_
#define ACADOS_UTILS_TOOLS_H_

#ifdef __cplusplus
extern "C" {
#endif

void dgemm_nn_3l(int m, int n, int k, double *A, int lda , double *B, int ldb, double *C, int ldc);
void dgemv_n_3l(int m, int n, double *A, int lda, double *x, double *y);
void dgemv_t_3l(int m, int n, double *A, int lda, double *x, double *y);
void dcopy_3l(int n, double *x, int incx, double *y, int incy);
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

void dger_3l(int m, int n, double alpha, double *x, int incx, double *y,
             int incy, double *A, int lda);

void dgetf2_3l(int m, int n, double *A, int lda, int *ipiv, int *info);

void dlaswp_3l(int n, double *A, int lda, int k1, int k2, int *ipiv);

void dtrsm_l_l_n_u_3l(int m, int n, double *A, int lda, double *B, int ldb);

void dgetrs_3l(int n, int nrhs, double *A, int lda, int *ipiv, double *B,
               int ldb);

void dgesv_3l(int n, int nrhs, double *A, int lda, int *ipiv, double *B,
                             int ldb, int *info);

double onenorm(int row, int col, double *ptrA);

double twonormv(int n, double *ptrv);

void padeapprox(int m, int row, double *A);

void expm(int row, double *A);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_UTILS_TOOLS_H_
