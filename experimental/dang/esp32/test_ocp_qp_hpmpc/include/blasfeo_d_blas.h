/**************************************************************************************************
* acados/external/blasfeo/include/blasfeo_d_blas.h                                                                                                *
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



#ifdef __cplusplus
extern "C" {
#endif



//
// level 1 BLAS
//

// y = y + alpha*x
void daxpy_lib(int kmax, double alpha, double *x, double *y);
void daxpy_libstr(int kmax, double alpha, struct d_strvec *sx, int xi, struct d_strvec *sy, int yi);
// z = y, y = y + alpha*x
void daxpy_bkp_lib(int kmax, double alpha, double *x, double *y, double *z);
void daxpy_bkp_libstr(int m, double alpha, struct d_strvec *sx, int xi, struct d_strvec *sy, int yi, struct d_strvec *sz, int zi);



//
// level 2 BLAS
//

// z <= beta * y + alpha * A * x
void dgemv_n_lib(int m, int n, double alpha, double *pA, int sda, double *x, double beta, double *y, double *z);
void dgemv_n_libstr(int m, int n, double alpha, struct d_strmat *sA, int ai, int aj, struct d_strvec *sx, int xi, double beta, struct d_strvec *sy, int yi, struct d_strvec *sz, int zi);
// z <= beta * y + alpha * A' * x
void dgemv_t_lib(int m, int n, double alpha, double *pA, int sda, double *x, double beta, double *y, double *z);
void dgemv_t_libstr(int m, int n, double alpha, struct d_strmat *sA, int ai, int aj, struct d_strvec *sx, int xi, double beta, struct d_strvec *sy, int yi, struct d_strvec *sz, int zi);
// y <= inv( A ) * x
void dtrsv_ln_inv_lib(int m, int n, double *pA, int sda, double *inv_diag_A, double *x, double *y);
void dtrsv_lnn_libstr(int m, int n, struct d_strmat *sA, int ai, int aj, struct d_strvec *sx, int xi, struct d_strvec *sz, int zi);
// y <= inv( A' ) * x
void dtrsv_lt_inv_lib(int m, int n, double *pA, int sda, double *inv_diag_A, double *x, double *y);
void dtrsv_ltn_libstr(int m, int n, struct d_strmat *sA, int ai, int aj, struct d_strvec *sx, int xi, struct d_strvec *sz, int zi);
// z <= beta * y + alpha * A * x ; A upper triangular
void dtrmv_un_lib(int m, double *pA, int sda, double *x, int alg, double *y, double *z);
void dtrmv_unn_libstr(int m, struct d_strmat *sA, int ai, int aj, struct d_strvec *sx, int xi, struct d_strvec *sz, int zi);
// z <= beta * y + alpha * A' * x ; A upper triangular
void dtrmv_ut_lib(int m, double *pA, int sda, double *x, int alg, double *y, double *z);
void dtrmv_utn_libstr(int m, struct d_strmat *sA, int ai, int aj, struct d_strvec *sx, int xi, struct d_strvec *sz, int zi);
// z <= beta * y + alpha * A * x ; A lower triangular
void dtrmv_lnn_libstr(int m, struct d_strmat *sA, int ai, int aj, struct d_strvec *sx, int xi, struct d_strvec *sz, int zi);
// z <= beta * y + alpha * A' * x ; A lower triangular
void dtrmv_ltn_libstr(int m, struct d_strmat *sA, int ai, int aj, struct d_strvec *sx, int xi, struct d_strvec *sz, int zi);
// z_n <= beta_n * y_n + alpha_n * A  * x_n
// z_t <= beta_t * y_t + alpha_t * A' * x_t
void dgemv_nt_lib(int m, int n, double alpha_n, double alpha_t, double *pA, int sda, double *x_n, double *x_t, double beta_n, double beta_t, double *y_n, double *y_t, double *z_n, double *z_t);
void dgemv_nt_libstr(int m, int n, double alpha_n, double alpha_t, struct d_strmat *sA, int ai, int aj, struct d_strvec *sx_n, int xi_n, struct d_strvec *sx_t, int xi_t, double beta_n, double beta_t, struct d_strvec *sy_n, int yi_n, struct d_strvec *sy_t, int yi_t, struct d_strvec *sz_n, int zi_n, struct d_strvec *sz_t, int zi_t);
// z <= beta * y + alpha * A * x, where A is symmetric and only the lower triangular patr of A is accessed
void dsymv_l_lib(int m, int n, double alpha, double *pA, int sda, double *x, double beta, double *y, double *z);
void dsymv_l_libstr(int m, int n, double alpha, struct d_strmat *sA, int ai, int aj, struct d_strvec *sx, int xi, double beta, struct d_strvec *sy, int yi, struct d_strvec *sz, int zi);



//
// level 3 BLAS
//

// dense

// D <= beta * C + alpha * A * B'
void dgemm_nt_lib(int m, int n, int k, double alpha, double *pA, int sda, double *pB, int sdb, double beta, double *pC, int sdc, double *pD, int sdd);
void dgemm_nt_libstr(int m, int n, int k, double alpha, struct d_strmat *sA, int ai, int aj, struct d_strmat *sB, int bi, int bj, double beta, struct d_strmat *sC, int ci, int cj, struct d_strmat *sD, int di, int dj);
// D <= beta * C + alpha * A * B
void dgemm_nn_lib(int m, int n, int k, double alpha, double *pA, int sda, double *pB, int sdb, double beta, double *pC, int sdc, double *pD, int sdd);
void dgemm_nn_libstr(int m, int n, int k, double alpha, struct d_strmat *sA, int ai, int aj, struct d_strmat *sB, int bi, int bj, double beta, struct d_strmat *sC, int ci, int cj, struct d_strmat *sD, int di, int dj);
// D <= beta * C + alpha * A * B' ; C, D lower triangular
void dsyrk_nt_l_lib(int m, int n, int k, double alpha, double *pA, int sda, double *pB, int sdb, double beta, double *pC, int sdc, double *pD, int sdd);
void dsyrk_ln_libstr(int m, int n, int k, double alpha, struct d_strmat *sA, int ai, int aj, struct d_strmat *sB, int bi, int bj, double beta, struct d_strmat *sC, int ci, int cj, struct d_strmat *sD, int di, int dj);
// D <= beta * C + alpha * A * B' ; B upper triangular
void dtrmm_nt_ru_lib(int m, int n, double alpha, double *pA, int sda, double *pB, int sdb, double beta, double *pC, int sdc, double *pD, int sdd);
void dtrmm_rutn_libstr(int m, int n, double alpha, struct d_strmat *sA, int ai, int aj, struct d_strmat *sB, int bi, int bj, double beta, struct d_strmat *sC, int ci, int cj, struct d_strmat *sD, int di, int dj);
// D <= beta * C + alpha * A * B ; B lower triangular
//void dtrmm_nt_ru_lib(int m, int n, double alpha, double *pA, int sda, double *pB, int sdb, double beta, double *pC, int sdc, double *pD, int sdd);
void dtrmm_rlnn_libstr(int m, int n, double alpha, struct d_strmat *sA, int ai, int aj, struct d_strmat *sB, int bi, int bj, double beta, struct d_strmat *sC, int ci, int cj, struct d_strmat *sD, int di, int dj);
// D <= B * A^{-T} , with A lower triangular employeing explicit inverse of diagonal
void dtrsm_nt_rl_inv_lib(int m, int n, double *pA, int sda, double *inv_diag_A, double *pB, int sdb, double *pD, int sdd);
void dtrsm_rltn_libstr(int m, int n, double alpha, struct d_strmat *sA, int ai, int aj, struct d_strmat *sB, int bi, int bj, struct d_strmat *sD, int di, int dj);
// D <= B * A^{-T} , with A lower triangular with unit diagonal
void dtrsm_nt_rl_one_lib(int m, int n, double *pA, int sda, double *pB, int sdb, double *pD, int sdd);
void dtrsm_rltu_libstr(int m, int n, double alpha, struct d_strmat *sA, int ai, int aj, struct d_strmat *sB, int bi, int bj, struct d_strmat *sD, int di, int dj);
// D <= B * A^{-T} , with A upper triangular employing explicit inverse of diagonal
void dtrsm_nt_ru_inv_lib(int m, int n, double *pA, int sda, double *inv_diag_A, double *pB, int sdb, double *pD, int sdd);
void dtrsm_rutn_libstr(int m, int n, double alpha, struct d_strmat *sA, int ai, int aj, struct d_strmat *sB, int bi, int bj, struct d_strmat *sD, int di, int dj);
// D <= A^{-1} * B , with A lower triangular with unit diagonal
void dtrsm_nn_ll_one_lib(int m, int n, double *pA, int sda, double *pB, int sdb, double *pD, int sdd);
void dtrsm_llnu_libstr(int m, int n, double alpha, struct d_strmat *sA, int ai, int aj, struct d_strmat *sB, int bi, int bj, struct d_strmat *sD, int di, int dj);
// D <= A^{-1} * B , with A upper triangular employing explicit inverse of diagonal
void dtrsm_nn_lu_inv_lib(int m, int n, double *pA, int sda, double *inv_diag_A, double *pB, int sdb, double *pD, int sdd);
void dtrsm_lunn_libstr(int m, int n, double alpha, struct d_strmat *sA, int ai, int aj, struct d_strmat *sB, int bi, int bj, struct d_strmat *sD, int di, int dj);

// diagonal

// D <= alpha * A * B + beta * C, with A diagonal (stored as strvec)
void dgemm_diag_left_ib(int m, int n, double alpha, double *dA, double *pB, int sdb, double beta, double *pC, int sdc, double *pD, int sdd);
void dgemm_l_diag_libstr(int m, int n, double alpha, struct d_strvec *sA, int ai, struct d_strmat *sB, int bi, int bj, double beta, struct d_strmat *sC, int ci, int cj, struct d_strmat *sD, int di, int dj);
// D <= alpha * A * B + beta * C, with B diagonal (stored as strvec)
void dgemm_diag_right_lib(int m, int n, double alpha, double *pA, int sda, double *dB, double beta, double *pC, int sdc, double *pD, int sdd);
void dgemm_r_diag_libstr(int m, int n, double alpha, struct d_strmat *sA, int ai, int aj, struct d_strvec *sB, int bi, double beta, struct d_strmat *sC, int ci, int cj, struct d_strmat *sD, int di, int dj);



//
// LAPACK
//

// D <= chol( C ) ; C, D lower triangular
void dpotrf_nt_l_lib(int m, int n, double *pC, int sdc, double *pD, int sdd, double *inv_diag_D);
void dpotrf_l_libstr(int m, int n, struct d_strmat *sC, int ci, int cj, struct d_strmat *sD, int di, int dj);
// D <= chol( C + A * B' ) ; C, D lower triangular
void dsyrk_dpotrf_nt_l_lib(int m, int n, int k, double *pA, int sda, double *pB, int sdb, double *pC, int sdc, double *pD, int sdd, double *inv_diag_D);
void dsyrk_dpotrf_ln_libstr(int m, int n, int k, struct d_strmat *sA, int ai, int aj, struct d_strmat *sB, int bi, int bj, struct d_strmat *sC, int ci, int cj, struct d_strmat *sD, int di, int dj);
// D <= lu( C ) ; no pivoting
void dgetrf_nn_nopivot_lib(int m, int n, double *pC, int sdc, double *pD, int sdd, double *inv_diag_D);
void dgetrf_nopivot_libstr(int m, int n, struct d_strmat *sC, int ci, int cj, struct d_strmat *sD, int di, int dj);
// D <= lu( C ) ; pivoting
void dgetrf_nn_lib(int m, int n, double *pC, int sdc, double *pD, int sdd, double *inv_diag_D, int *ipiv);
void dgetrf_libstr(int m, int n, struct d_strmat *sC, int ci, int cj, struct d_strmat *sD, int di, int dj, int *ipiv);



#ifdef __cplusplus
}
#endif
