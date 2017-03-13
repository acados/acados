/**************************************************************************************************
* acados/external/blasfeo/include/blasfeo_d_kernel.h                                                                                                *
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



// level 2 BLAS
// 12
void kernel_dgemv_n_12_lib4(int k, double *alpha, double *A, int sda, double *x, double *beta, double *y, double *z);
void kernel_dgemv_t_12_lib4(int k, double *alpha, double *A, int sda, double *x, double *beta, double *y, double *z);
// 8
void kernel_dgemv_n_8_lib4(int k, double *alpha, double *A, int sda, double *x, double *beta, double *y, double *z);
void kernel_dgemv_t_8_lib4(int k, double *alpha, double *A, int sda, double *x, double *beta, double *y, double *z);
void kernel_dtrmv_un_8_lib4(int k, double *A, int sda, double *x, int alg, double *y, double *z);
// 4
void kernel_dgemv_n_4_lib4(int k, double *alpha, double *A, double *x, double *beta, double *y, double *z);
void kernel_dgemv_n_4_vs_lib4(int k, double *alpha, double *A, double *x, double *beta, double *y, double *z, int km);
void kernel_dgemv_t_4_lib4(int k, double *alpha, double *A, int sda, double *x, double *beta, double *y, double *z);
void kernel_dgemv_t_4_vs_lib4(int k, double *alpha, double *A, int sda, double *x, double *beta, double *C, double *D, int km);
void kernel_dtrsv_ln_inv_4_lib4(int k, double *A, double *inv_diag_A, double *x, double *y, double *z);
void kernel_dtrsv_ln_inv_4_vs_lib4(int k, double *A, double *inv_diag_A, double *x, double *y, double *z, int km, int kn);
void kernel_dtrsv_lt_inv_4_lib4(int k, double *A, int sda, double *inv_diag_A, double *x, double *y, double *z);
void kernel_dtrsv_lt_inv_3_lib4(int k, double *A, int sda, double *inv_diag_A, double *x, double *y, double *z);
void kernel_dtrsv_lt_inv_2_lib4(int k, double *A, int sda, double *inv_diag_A, double *x, double *y, double *z);
void kernel_dtrsv_lt_inv_1_lib4(int k, double *A, int sda, double *inv_diag_A, double *x, double *y, double *z);
void kernel_dtrmv_un_4_lib4(int k, double *A, double *x, int alg, double *y, double *z);
void kernel_dtrmv_ut_4_lib4(int k, double *A, int sda, double *x, int alg, double *y, double *z);
void kernel_dtrmv_ut_4_vs_lib4(int k, double *A, int sda, double *x, int alg, double *C, double *D, int km);
void kernel_dgemv_nt_6_lib4(int kmax, double *alpha_n, double *alpha_t, double *A, int sda, double *x_n, double *x_t, double *beta_t, double *y_t, double *z_n, double *z_t);
void kernel_dgemv_nt_4_lib4(int kmax, double *alpha_n, double *alpha_t, double *A, int sda, double *x_n, double *x_t, double *beta_t, double *y_t, double *z_n, double *z_t);
void kernel_dgemv_nt_4_vs_lib4(int kmax, double *alpha_n, double *alpha_t, double *A, int sda, double *x_n, double *x_t, double *beta_t, double *y_t, double *z_n, double *z_t, int km);
void kernel_dsymv_l_4_lib4(int kmax, double *alpha, double *A, int sda, double *x_n, double *x_t, double *z_n, double *z_t);
void kernel_dsymv_l_4_vs_lib4(int kmax, double *alpha, double *A, int sda, double *x_n, double *x_t, double *z_n, double *z_t, int km);



// level 3 BLAS
// 12x4
void kernel_dgemm_nt_12x4_a0_lib4(int k, double *alpha, double *A, int sda, double *B, double *D, int sdd); //
void kernel_dgemm_nt_12x4_lib4(int k, double *alpha, double *A, int sda, double *B, double *beta, double *C, int sdc, double *D, int sdd); //
void kernel_dgemm_nt_12x4_vs_lib4(int k, double *alpha, double *A, int sda, double *B, double *beta, double *C, int sdc, double *D, int sdd, int km, int kn); //
void kernel_dgemm_nn_12x4_lib4(int k, double *alpha, double *A, int sda, double *B, int sdb, double *beta, double *C, int sdc, double *D, int sdd); //
void kernel_dgemm_nn_12x4_vs_lib4(int k, double *alpha, double *A, int sda, double *B, int sdb, double *beta, double *C, int sdc, double *D, int sdd, int km, int kn); //
void kernel_dsyrk_nt_l_12x4_lib4(int k, double *alpha, double *A, int sda, double *B, double *beta, double *C, int sdc, double *D, int sdd); //
void kernel_dsyrk_nt_l_12x4_vs_lib4(int k, double *alpha, double *A, int sda, double *B, double *beta, double *C, int sdc, double *D, int sdd, int km, int kn); //
void kernel_dtrmm_nt_ru_12x4_lib4(int k, double *alpha, double *A, int sda, double *B, double *beta, double *C, int sdc, double *D, int sdd); //
void kernel_dtrmm_nt_ru_12x4_vs_lib4(int k, double *alpha, double *A, int sda, double *B, double *beta, double *C, int sdc, double *D, int sdd, int km, int kn); //
void kernel_dtrsm_nt_rl_inv_12x4_vs_lib4(int k, double *A, int sda, double *B, double *C, int sdc, double *D, int sdd, double *E, double *inv_diag_E, int km, int kn);
void kernel_dtrsm_nt_rl_inv_12x4_lib4(int k, double *A, int sda, double *B, double *C, int sdc, double *D, int sdd, double *E, double *inv_diag_E);
void kernel_dtrsm_nt_rl_one_12x4_lib4(int k, double *A, int sda, double *B, double *C, int sdc, double *D, int sdd, double *E);
void kernel_dtrsm_nt_rl_one_12x4_vs_lib4(int k, double *A, int sda, double *B, double *C, int sdc, double *D, int sdd, double *E, int km, int kn);
void kernel_dtrsm_nt_ru_inv_12x4_lib4(int k, double *A, int sda, double *B, double *C, int sdc, double *D, int sdd, double *E, double *inv_diag_E);
void kernel_dtrsm_nt_ru_inv_12x4_vs_lib4(int k, double *A, int sda, double *B, double *C, int sdc, double *D, int sdd, double *E, double *inv_diag_E, int km, int kn);
void kernel_dtrsm_nn_ru_inv_12x4_lib4(int k, double *A, int sda, double *B, int sdb, double *C, int sdc, double *D, int sdd, double *E, double *inv_diag_E);
void kernel_dtrsm_nn_ru_inv_12x4_vs_lib4(int k, double *A, int sda, double *B, int sdb, double *C, int sdc, double *D, int sdd, double *E, double *inv_diag_E, int km, int kn);
void kernel_dtrsm_nn_ll_one_12x4_lib4(int k, double *A, int sda, double *B, int sdb, double *C, int sdc, double *D, int sdd, double *E, int sde);
void kernel_dtrsm_nn_ll_one_12x4_vs_lib4(int k, double *A, int sda, double *B, int sdb, double *C, int sdc, double *D, int sdd, double *E, int sde, int km, int kn);
void kernel_dtrsm_nn_lu_inv_12x4_lib4(int kmax, double *A, int sda, double *B, int sdb, double *C, int sdc, double *D, int sdd, double *E, int sde, double *inv_diag_E);
void kernel_dtrsm_nn_lu_inv_12x4_vs_lib4(int kmax, double *A, int sda, double *B, int sdb, double *C, int sdc, double *D, int sdd, double *E, int sde, double *inv_diag_E, int km, int kn);
// 8x4
void kernel_dgemm_nt_8x4_a0_lib4(int k, double *alpha, double *A, int sda, double *B, double *D, int sdd); //
void kernel_dgemm_nt_8x4_lib4(int k, double *alpha, double *A, int sda, double *B, double *beta, double *C, int sdc, double *D, int sdd); //
void kernel_dgemm_nt_8x4_vs_lib4(int k, double *alpha, double *A, int sda, double *B, double *beta, double *C, int sdc, double *D, int sdd, int km, int kn); //
void kernel_dgemm_nt_8x4_gen_lib4(int k, double *alpha, double *A, int sda, double *B, double *beta, int offsetC, double *C, int sdc, int offsetD, double *D, int sdd, int m0, int m1, int k0, int k1);
void kernel_dgemm_nn_8x4_lib4(int k, double *alpha, double *A, int sda, double *B, int sdb, double *beta, double *C, int sdc, double *D, int sdd); //
void kernel_dgemm_nn_8x4_vs_lib4(int k, double *alpha, double *A, int sda, double *B, int sdb, double *beta, double *C, int sdc, double *D, int sdd, int km, int kn); //
void kernel_dsyrk_nt_l_8x4_lib4(int k, double *alpha, double *A, int sda, double *B, double *beta, double *C, int sdc, double *D, int sdd); //
void kernel_dsyrk_nt_l_8x4_vs_lib4(int k, double *alpha, double *A, int sda, double *B, double *beta, double *C, int sdc, double *D, int sdd, int km, int kn); //
void kernel_dtrmm_nt_ru_8x4_lib4(int k, double *alpha, double *A, int sda, double *B, double *beta, double *C, int sdc, double *D, int sdd); //
void kernel_dtrmm_nt_ru_8x4_vs_lib4(int k, double *alpha, double *A, int sda, double *B, double *beta, double *C, int sdc, double *D, int sdd, int km, int kn); //
void kernel_dtrsm_nt_rl_inv_8x4_vs_lib4(int k, double *A, int sda, double *B, double *C, int sdc, double *D, int sdd, double *E, double *inv_diag_E, int km, int kn);
void kernel_dtrsm_nt_rl_inv_8x4_lib4(int k, double *A, int sda, double *B, double *C, int sdc, double *D, int sdd, double *E, double *inv_diag_E);
void kernel_dtrsm_nt_rl_one_8x4_lib4(int k, double *A, int sda, double *B, double *C, int sdc, double *D, int sdd, double *E);
void kernel_dtrsm_nt_rl_one_8x4_vs_lib4(int k, double *A, int sda, double *B, double *C, int sdc, double *D, int sdd, double *E, int km, int kn);
void kernel_dtrsm_nt_ru_inv_8x4_lib4(int k, double *A, int sda, double *B, double *C, int sdc, double *D, int sdd, double *E, double *inv_diag_E);
void kernel_dtrsm_nt_ru_inv_8x4_vs_lib4(int k, double *A, int sda, double *B, double *C, int sdc, double *D, int sdd, double *E, double *inv_diag_E, int km, int kn);
void kernel_dtrsm_nn_ru_inv_8x4_lib4(int k, double *A, int sda, double *B, int sdb, double *C, int sdc, double *D, int sdd, double *E, double *inv_diag_E);
void kernel_dtrsm_nn_ru_inv_8x4_vs_lib4(int k, double *A, int sda, double *B, int sdb, double *C, int sdc, double *D, int sdd, double *E, double *inv_diag_E, int km, int kn);
void kernel_dtrsm_nn_ll_one_8x4_lib4(int k, double *A, int sda, double *B, int sdb, double *C, int sdc, double *D, int sdd, double *E, int sde);
void kernel_dtrsm_nn_ll_one_8x4_vs_lib4(int k, double *A, int sda, double *B, int sdb, double *C, int sdc, double *D, int sdd, double *E, int sde, int km, int kn);
void kernel_dtrsm_nn_lu_inv_8x4_lib4(int kmax, double *A, int sda, double *B, int sdb, double *C, int sdc, double *D, int sdd, double *E, int sde, double *inv_diag_E);
void kernel_dtrsm_nn_lu_inv_8x4_vs_lib4(int kmax, double *A, int sda, double *B, int sdb, double *C, int sdc, double *D, int sdd, double *E, int sde, double *inv_diag_E, int km, int kn);
// 4x4
void kernel_dgemm_nt_4x4_a0_lib4(int k, double *alpha, double *A, double *B, double *D); //
void kernel_dgemm_nt_4x4_lib4(int k, double *alpha, double *A, double *B, double *beta, double *C, double *D); //
void kernel_dgemm_nt_4x4_vs_lib4(int k, double *alpha, double *A, double *B, double *beta, double *C, double *D, int km, int kn); //
void kernel_dgemm_nt_4x4_gen_lib4(int k, double *alpha, double *A, double *B, double *beta, int offsetC, double *C, int sdc, int offsetD, double *D, int sdd, int m0, int m1, int k0, int k1);
void kernel_dgemm_nn_4x4_lib4(int k, double *alpha, double *A, double *B, int sdb, double *beta, double *C, double *D); //
void kernel_dgemm_nn_4x4_vs_lib4(int k, double *alpha, double *A, double *B, int sdb, double *beta, double *C, double *D, int km, int kn); //
void kernel_dsyrk_nt_l_4x4_lib4(int k, double *alpha, double *A, double *B, double *beta, double *C, double *D); //
void kernel_dsyrk_nt_l_4x4_vs_lib4(int k, double *alpha, double *A, double *B, double *beta, double *C, double *D, int km, int kn); //
void kernel_dtrmm_nt_ru_4x4_lib4(int k, double *alpha, double *A, double *B, double *beta, double *C, double *D); //
void kernel_dtrmm_nt_ru_4x4_vs_lib4(int k, double *alpha, double *A, double *B, double *beta, double *C, double *D, int km, int kn); //
void kernel_dtrsm_nt_rl_inv_4x4_lib4(int k, double *A, double *B, double *C, double *D, double *E, double *inv_diag_E);
void kernel_dtrsm_nt_rl_inv_4x4_vs_lib4(int k, double *A, double *B, double *C, double *D, double *E, double *inv_diag_E, int km, int kn);
void kernel_dtrsm_nt_rl_one_4x4_lib4(int k, double *A, double *B, double *C, double *D, double *E);
void kernel_dtrsm_nt_rl_one_4x4_vs_lib4(int k, double *A, double *B, double *C, double *D, double *E, int km, int kn);
void kernel_dtrsm_nt_ru_inv_4x4_lib4(int k, double *A, double *B, double *C, double *D, double *E, double *inv_diag_E);
void kernel_dtrsm_nt_ru_inv_4x4_vs_lib4(int k, double *A, double *B, double *C, double *D, double *E, double *inv_diag_E, int km, int kn);
void kernel_dtrsm_nn_ru_inv_4x4_lib4(int k, double *A, double *B, int sdb, double *C, double *D, double *E, double *inv_diag_E);
void kernel_dtrsm_nn_ru_inv_4x4_vs_lib4(int k, double *A, double *B, int sdb, double *C, double *D, double *E, double *inv_diag_E, int km, int kn);
void kernel_dtrsm_nn_ll_one_4x4_lib4(int k, double *A, double *B, int sdb, double *C, double *D, double *E);
void kernel_dtrsm_nn_ll_one_4x4_vs_lib4(int k, double *A, double *B, int sdb, double *C, double *D, double *E, int km, int kn);
void kernel_dtrsm_nn_lu_inv_4x4_lib4(int kmax, double *A, double *B, int sdb, double *C, double *D, double *E, double *inv_diag_E);
void kernel_dtrsm_nn_lu_inv_4x4_vs_lib4(int kmax, double *A, double *B, int sdb, double *C, double *D, double *E, double *inv_diag_E, int km, int kn);
// diag
void kernel_dgemm_diag_right_4_a0_lib4(int kmax, double *alpha, double *A, int sda, double *B, double *D, int sdd);
void kernel_dgemm_diag_right_4_lib4(int kmax, double *alpha, double *A, int sda, double *B, double *beta, double *C, int sdc, double *D, int sdd);
void kernel_dgemm_diag_right_3_lib4(int kmax, double *alpha, double *A, int sda, double *B, double *beta, double *C, int sdc, double *D, int sdd);
void kernel_dgemm_diag_right_2_lib4(int kmax, double *alpha, double *A, int sda, double *B, double *beta, double *C, int sdc, double *D, int sdd);
void kernel_dgemm_diag_right_1_lib4(int kmax, double *alpha, double *A, int sda, double *B, double *beta, double *C, int sdc, double *D, int sdd);
void kernel_dgemm_diag_left_4_a0_lib4(int kmax, double *alpha, double *A, double *B, double *D);
void kernel_dgemm_diag_left_4_lib4(int kmax, double *alpha, double *A, double *B, double *beta, double *C, double *D);
void kernel_dgemm_diag_left_3_lib4(int kmax, double *alpha, double *A, double *B, double *beta, double *C, double *D);
void kernel_dgemm_diag_left_2_lib4(int kmax, double *alpha, double *A, double *B, double *beta, double *C, double *D);
void kernel_dgemm_diag_left_1_lib4(int kmax, double *alpha, double *A, double *B, double *beta, double *C, double *D);



// LAPACK
// 12x4
void kernel_dpotrf_nt_l_12x4_lib4(int k, double *A, int sda, double *B, double *C, int sdc, double *D, int sdd, double *inv_diag_D);
void kernel_dpotrf_nt_l_12x4_vs_lib4(int k, double *A, int sda, double *B, double *C, int sdc, double *D, int sdd, double *inv_diag_D, int km, int kn);
void kernel_dgetrf_nn_l_12x4_lib4(int k, double *A, int sda, double *B, int sdb, double *C, int sdc, double *D, int sdd, double *inv_diag_D);
void kernel_dgetrf_nn_m_12x4_lib4(int k, double *A, int sda, double *B, int sdb, double *C, int sdc, double *D, int sdd, double *inv_diag_D);
void kernel_dgetrf_nn_r_12x4_lib4(int k, double *A, int sda, double *B, int sdb, double *C, int sdc, double *D, int sdd, double *inv_diag_D);
void kernel_dgetrf_nn_l_12x4_vs_lib4(int k, double *A, int sda, double *B, int sdb, double *C, int sdc, double *D, int sdd, double *inv_diag_D, int km, int kn);
void kernel_dgetrf_nn_m_12x4_vs_lib4(int k, double *A, int sda, double *B, int sdb, double *C, int sdc, double *D, int sdd, double *inv_diag_D, int km, int kn);
void kernel_dgetrf_nn_r_12x4_vs_lib4(int k, double *A, int sda, double *B, int sdb, double *C, int sdc, double *D, int sdd, double *inv_diag_D, int km, int kn);
// 8x4
void kernel_dpotrf_nt_l_8x4_lib4(int k, double *A, int sda, double *B, double *C, int sdc, double *D, int sdd, double *inv_diag_D);
void kernel_dpotrf_nt_l_8x4_vs_lib4(int k, double *A, int sda, double *B, double *C, int sdc, double *D, int sdd, double *inv_diag_D, int km, int kn);
void kernel_dgetrf_nn_l_8x4_lib4(int k, double *A, int sda, double *B, int sdb, double *C, int sdc, double *D, int sdd, double *inv_diag_D);
void kernel_dgetrf_nn_r_8x4_lib4(int k, double *A, int sda, double *B, int sdb, double *C, int sdc, double *D, int sdd, double *inv_diag_D);
void kernel_dgetrf_nn_l_8x4_vs_lib4(int k, double *A, int sda, double *B, int sdb, double *C, int sdc, double *D, int sdd, double *inv_diag_D, int km, int kn);
void kernel_dgetrf_nn_r_8x4_vs_lib4(int k, double *A, int sda, double *B, int sdb, double *C, int sdc, double *D, int sdd, double *inv_diag_D, int km, int kn);
// 4x4
void kernel_dpotrf_nt_l_4x4_lib4(int k, double *A, double *B, double *C, double *D, double *inv_diag_D);
void kernel_dpotrf_nt_l_4x4_vs_lib4(int k, double *A, double *B, double *C, double *D, double *inv_diag_D, int km, int kn);
#if defined(TARGET_X64_INTEL_SANDY_BRIDGE)
void kernel_dlauum_nt_4x4_lib4(int k, double *alpha, double *A, double *B, double *beta, double *C, double *D); //
void kernel_dlauum_nt_4x4_vs_lib4(int k, double *alpha, double *A, double *B, double *beta, double *C, double *D, int km, int kn); //
#endif
void kernel_dgetrf_nn_4x4_lib4(int k, double *A, double *B, int sdb, double *C, double *D, double *inv_diag_D);
void kernel_dgetrf_nn_4x4_vs_lib4(int k, double *A, double *B, int sdb, double *C, double *D, double *inv_diag_D, int km, int kn);
void kernel_dgetrf_pivot_4_lib4(int m, double *pA, int sda, double *inv_diag_A, int* ipiv);
void kernel_dgetrf_pivot_4_vs_lib4(int m, int n, double *pA, int sda, double *inv_diag_A, int* ipiv);



// merged routines
// 12x4
void kernel_dgemm_dtrsm_nt_rl_inv_12x4_vs_lib4(int kp, double *Ap, int sdap, double *Bp, int km_, double *Am, int sdam, double *Bm, double *C, int sdc, double *D, int sdd, double *E, double *inv_diag_E, int km, int kn);
void kernel_dgemm_dtrsm_nt_rl_inv_12x4_lib4(int kp, double *Ap, int sdap, double *Bp, int km_, double *Am, int sdam, double *Bm, double *C, int sdc, double *D, int sdd, double *E, double *inv_diag_E);
void kernel_dsyrk_dpotrf_nt_l_12x4_vs_lib4(int kp, double *Ap, int sdap, double *Bp, int km_, double *Am, int sdam, double *Bm, double *C, int sdc, double *D, int sdd, double *inv_diag_D, int km, int kn);
void kernel_dsyrk_dpotrf_nt_l_12x4_lib4(int kp, double *Ap, int sdap, double *Bp, int km_, double *Am, int sdam, double *Bm, double *C, int sdc, double *D, int sdd, double *inv_diag_D);
// 8x4
void kernel_dgemm_dtrsm_nt_rl_inv_8x4_vs_lib4(int kp, double *Ap, int sdap, double *Bp, int km_, double *Am, int sdam, double *Bm, double *C, int sdc, double *D, int sdd, double *E, double *inv_diag_E, int km, int kn);
void kernel_dgemm_dtrsm_nt_rl_inv_8x4_lib4(int kp, double *Ap, int sdap, double *Bp, int km_, double *Am, int sdam, double *Bm, double *C, int sdc, double *D, int sdd, double *E, double *inv_diag_E);
void kernel_dsyrk_dpotrf_nt_l_8x4_lib4(int kp, double *Ap, int sdap, double *Bp, int km_, double *Am, int sdam, double *Bm, double *C, int sdc, double *D, int sdd, double *inv_diag_D);
void kernel_dsyrk_dpotrf_nt_l_8x4_vs_lib4(int kp, double *Ap, int sdap, double *Bp, int km_, double *Am, int sdam, double *Bm, double *C, int sdc, double *D, int sdd, double *inv_diag_D, int km, int kn);
// 4x4
void kernel_dgemm_dtrsm_nt_rl_inv_4x4_lib4(int kp, double *Ap, double *Bp, int km_, double *Am, double *Bm, double *C, double *D, double *E, double *inv_diag_E);
void kernel_dgemm_dtrsm_nt_rl_inv_4x4_vs_lib4(int kp, double *Ap, double *Bp, int km_, double *Am, double *Bm, double *C, double *D, double *E, double *inv_diag_E, int km, int kn);
void kernel_dsyrk_dpotrf_nt_l_4x4_vs_lib4(int kp, double *Ap, double *Bp, int km_, double *Am, double *Bm, double *C, double *D, double *inv_diag_D, int km, int kn);
void kernel_dsyrk_dpotrf_nt_l_4x4_lib4(int kp, double *Ap, double *Bp, int km_, double *Am, double *Bm, double *C, double *D, double *inv_diag_D);



// auxiliary routines
void kernel_dgecp_8_0_lib4(int tri, int kmax, double alpha, double *A0, int sda,  double *B0, int sdb);
void kernel_dgecp_8_1_lib4(int tri, int kmax, double alpha, double *A0, int sda, double *B0, int sdb);
void kernel_dgecp_8_2_lib4(int tri, int kmax, double alpha, double *A0, int sda, double *B0, int sdb);
void kernel_dgecp_8_3_lib4(int tri, int kmax, double alpha, double *A0, int sda, double *B0, int sdb);
void kernel_dgecp_4_0_lib4(int tri, int kmax, double alpha, double *A, double *B);
void kernel_dgecp_4_1_lib4(int tri, int kmax, double alpha, double *A0, int sda, double *B);
void kernel_dgecp_4_2_lib4(int tri, int kmax, double alpha, double *A0, int sda, double *B);
void kernel_dgecp_4_3_lib4(int tri, int kmax, double alpha, double *A0, int sda, double *B);
void kernel_dgecp_3_0_lib4(int tri, int kmax, double alpha, double *A, double *B);
void kernel_dgecp_3_2_lib4(int tri, int kmax, double alpha, double *A0, int sda, double *B);
void kernel_dgecp_3_3_lib4(int tri, int kmax, double alpha, double *A0, int sda, double *B);
void kernel_dgecp_2_0_lib4(int tri, int kmax, double alpha, double *A, double *B);
void kernel_dgecp_2_3_lib4(int tri, int kmax, double alpha, double *A0, int sda, double *B);
void kernel_dgecp_1_0_lib4(int tri, int kmax, double alpha, double *A, double *B);
void kernel_dgead_8_0_lib4(int kmax, double alpha, double *A0, int sda,  double *B0, int sdb);
void kernel_dgead_8_1_lib4(int kmax, double alpha, double *A0, int sda, double *B0, int sdb);
void kernel_dgead_8_2_lib4(int kmax, double alpha, double *A0, int sda, double *B0, int sdb);
void kernel_dgead_8_3_lib4(int kmax, double alpha, double *A0, int sda, double *B0, int sdb);
void kernel_dgead_4_0_lib4(int kmax, double alpha, double *A, double *B);
void kernel_dgead_4_1_lib4(int kmax, double alpha, double *A0, int sda, double *B);
void kernel_dgead_4_2_lib4(int kmax, double alpha, double *A0, int sda, double *B);
void kernel_dgead_4_3_lib4(int kmax, double alpha, double *A0, int sda, double *B);
void kernel_dgead_3_0_lib4(int kmax, double alpha, double *A, double *B);
void kernel_dgead_3_2_lib4(int kmax, double alpha, double *A0, int sda, double *B);
void kernel_dgead_3_3_lib4(int kmax, double alpha, double *A0, int sda, double *B);
void kernel_dgead_2_0_lib4(int kmax, double alpha, double *A, double *B);
void kernel_dgead_2_3_lib4(int kmax, double alpha, double *A0, int sda, double *B);
void kernel_dgead_1_0_lib4(int kmax, double alpha, double *A, double *B);
void kernel_dgeset_4_lib4(int kmax, double alpha, double *A);
void kernel_dtrset_4_lib4(int kmax, double alpha, double *A);
void kernel_dgetr_8_lib4(int tri, int kmax, int kna, double alpha, double *A, int sda, double *C, int sdc);
void kernel_dgetr_4_lib4(int tri, int kmax, int kna, double alpha, double *A, double *C, int sdc);
void kernel_dgetr_3_lib4(int tri, int kmax, int kna, double alpha, double *A, double *C, int sdc);
void kernel_dgetr_2_lib4(int tri, int kmax, int kna, double alpha, double *A, double *C, int sdc);
void kernel_dgetr_1_lib4(int tri, int kmax, int kna, double alpha, double *A, double *C, int sdc);



#ifdef __cplusplus
}
#endif
