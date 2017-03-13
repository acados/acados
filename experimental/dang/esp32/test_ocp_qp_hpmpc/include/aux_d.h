/**************************************************************************************************
* /home/dang/acados/external/hpmpc/include/aux_d.h                                                                                                *
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

#ifdef __cplusplus
extern "C" {
#endif

// 2017.03.13 Dang add target.h
#include "target.h"


#if ! defined(BLASFEO)
void d_zeros(double **pA, int row, int col);
void d_zeros_align(double **pA, int row, int col);
void v_zeros_align(void **pA, int size_in_bytes);
void d_free(double *pA);
void d_free_align(double *pA);
void v_free_align(void *pA);
void d_ones(double **pA, int row, int col);
void d_ones_align(double **pA, int row, int col);
#endif
void d_eye(double **pA, int row);
void d_rep_mat(int reps, int row, int col, double *A, int lda, double *B, int ldb);
void dadd_mat(int row, int col, double alpha, double *A, int lda, double *B, int ldb);
void dax_mat(int row, int col, double alpha, double *A, int lda, double *B, int ldb);
float d_max_mat(int row, int col, double *A, int lda);
float d_min_mat(int row, int col, double *A, int lda);
void d_set_mat(int row, int col, double alpha, double *A, int lda);
//void d_set_pmat(int row, int col, double alpha, int offset, double *pA, int sda);
void d_scale_mat(int row, int col, double alpha, double *A, int lda);
void d_scale_pmat(int row, int col, double alpha, int offset, double *pA, int sda);
void d_copy_mat(int row, int col, double *A, int lda, double *B, int ldb);
void d_tran_mat(int row, int col, double *A, int lda, double *B, int ldb);
void d_copy_pmat(int row, int col, int bs, double *A, int sda, double *B, int sdb);
void d_copy_pmat_l(int row, int bs, double *A, int sda, double *B, int sdb);
void d_copy_pmat_panel(int row, int col, int offset, double *A, double *B, int sdb);
void d_align_pmat_panel(int row, int col, int offset, double *A, int sda, double *B);
//void d_transpose_pmat_lo(int row, int offset, double *A, int sda, double *B, int sdb);
void d_align_pmat(int row, int col, int offset, int bs, double *A, int sda, double *B, int sdb);
#if ! defined(BLASFEO)
void d_cvt_mat2pmat(int row, int col, double *A, int lda, int offset, double *pA, int sda);
void d_cvt_tran_mat2pmat(int row, int col, double *A, int lda, int offset, double *pA, int sda);
void d_cvt_pmat2mat(int row, int col, int offset, double *pA, int sda, double *A, int lda);
void d_cvt_tran_pmat2mat(int row, int col, int offset, double *pA, int sda, double *A, int lda);
void d_print_mat(int row, int col, double *A, int lda);
void d_print_pmat(int row, int col, int bs, double *A, int sda);
#endif


// (new) routines for pmat (they are in _lib form)
void dgeset_lib(int row, int col, double alpha, int offset, double *pA, int sda);
void dtrset_lib(int row, double alpha, int offset, double *pA, int sda);




// integer stuff
void int_zeros(int **pA, int row, int col);
void int_free(int *pA);
void int_print_mat(int row, int col, int *A, int lda);



#ifdef __cplusplus
}
#endif
